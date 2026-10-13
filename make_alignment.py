#!/usr/bin/env python
"""
Step 2 - align the template with the full-length target and find the gaps.

  * aligns <code>.seq (template) with the full-length FASTA using Clustal
    Omega, or reads an alignment you already made (--aligned)
  * every gap in the TEMPLATE row = missing residues to be modelled
  * writes the Modeller .ali file and a small JSON file listing the gaps
    (read later by model-single.py for loop refinement)

Usage
  python make_alignment.py --seq 6I07.seq --target full_length.fasta
  python make_alignment.py --seq 6I07.seq --target full_length.fasta \
                           --aligned aligned.fasta
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import textwrap


def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(seq)))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        recs.append((name, ''.join(seq)))
    return recs


def parse_seq(path):
    """Modeller .seq file -> (code, sequence, first, chain, last)."""
    with open(path) as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]
    head = [f.strip() for f in lines[1].split(':')]
    seq = ''.join(lines[2:]).replace('*', '').replace(' ', '')
    return head[1], seq.upper(), head[2], head[3], head[4]


def clustalo_align(code, template_seq, target_seq):
    exe = shutil.which('clustalo')
    if not exe:
        sys.exit('clustalo not found on PATH. Install it with\n'
                 '  conda install -c bioconda clustalo\n'
                 'or align elsewhere and pass the result with --aligned.')
    infile, outfile = code + '_two_seqs.fasta', code + '_aligned.fasta'
    with open(infile, 'w') as fh:
        fh.write('>template\n%s\n>target\n%s\n' % (template_seq, target_seq))
    subprocess.run([exe, '-i', infile, '-o', outfile, '--outfmt=fasta',
                    '--force'], check=True)
    return outfile


def load_alignment(path):
    recs = read_fasta(path)
    if len(recs) != 2:
        sys.exit('%s must contain exactly 2 sequences (template first, '
                 'target second).' % path)
    return recs[0][1].upper(), recs[1][1].upper()


def find_gaps(template_aln, target_aln):
    """Runs where the template has a gap and the target has residues,
    in TARGET numbering (1-based). Returns (gaps, target_length)."""
    gaps, run, idx = [], [], 0
    for tc, gc in zip(template_aln, target_aln):
        if gc != '-':
            idx += 1
        if tc == '-' and gc != '-':
            run.append(idx)
        elif run:
            gaps.append([run[0], run[-1]])
            run = []
    if run:
        gaps.append([run[0], run[-1]])
    return gaps, idx


def identity(template_aln, target_aln):
    same = total = 0
    for a, b in zip(template_aln, target_aln):
        if a != '-' and b != '-':
            total += 1
            same += (a == b)
    return (100.0 * same / total) if total else 0.0


def write_ali(path, code, target_name, tmpl_aln, target_aln,
              first, chain, last):
    def wrap(s):
        return '\n'.join(textwrap.wrap(s + '*', 75))
    with open(path, 'w') as fh:
        fh.write('>P1;%s\n' % code)
        fh.write('structureX:%s:%s:%s:%s:%s:::-1.00:-1.00\n'
                 % (code, first, chain, last, chain))
        fh.write(wrap(tmpl_aln) + '\n\n')
        fh.write('>P1;%s\n' % target_name)
        fh.write('sequence:%s:::::::0.00:0.00\n' % target_name)
        fh.write(wrap(target_aln) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--seq', required=True, help='<code>.seq from getseq.py')
    ap.add_argument('--target', required=True,
                    help='FASTA with the full-length target sequence')
    ap.add_argument('--aligned', help='aligned FASTA (template, then target); '
                                      'skips Clustal Omega')
    ap.add_argument('--target-name', help='default <code>_modelled')
    ap.add_argument('--ali', help='output .ali (default <code>_alignment.ali)')
    ap.add_argument('--gaps', help='output JSON (default <code>_gaps.json)')
    ap.add_argument('--min-identity', type=float, default=90.0,
                    help='warn if template/target identity is below this %%')
    args = ap.parse_args()

    code, tmpl_seq, first, chain, last = parse_seq(args.seq)
    target_name = args.target_name or code + '_modelled'
    ali_path = args.ali or code + '_alignment.ali'
    gaps_path = args.gaps or code + '_gaps.json'

    recs = read_fasta(args.target)
    if not recs:
        sys.exit('No sequence found in %s' % args.target)
    target_seq = recs[0][1].upper()

    aligned = args.aligned or clustalo_align(code, tmpl_seq, target_seq)
    tmpl_aln, target_aln = load_alignment(aligned)

    if tmpl_aln.replace('-', '') != tmpl_seq:
        sys.exit('Template row of %s does not match %s. Align the sequence '
                 'from the .seq file (residues with coordinates), not the '
                 'full-length one.' % (aligned, args.seq))
    if target_aln.replace('-', '') != target_seq:
        sys.exit('Target row of %s does not match %s.'
                 % (aligned, args.target))
    if len(tmpl_aln) != len(target_aln):
        sys.exit('Aligned rows have different lengths.')

    gaps, n_res = find_gaps(tmpl_aln, target_aln)
    ident = identity(tmpl_aln, target_aln)
    extra = sum(1 for a, b in zip(tmpl_aln, target_aln) if b == '-' and a != '-')

    print('Template residues with coordinates: %d' % len(tmpl_seq))
    print('Target length: %d' % n_res)
    print('Identity over aligned positions: %.1f%%' % ident)
    print('Missing regions (target numbering): %s'
          % (', '.join('%d-%d' % tuple(g) for g in gaps) or 'none'))
    print('Missing residues in total: %d'
          % sum(e - s + 1 for s, e in gaps))
    if ident < args.min_identity:
        print('WARNING: identity is below %.0f%%. The PDB and the target '
              'differ (mutations, another construct or species?). '
              'Modelling will continue, but check the alignment.'
              % args.min_identity)
    if extra:
        print('WARNING: %d template residues have no counterpart in the '
              'target (tags or insertions in the PDB).' % extra)

    write_ali(ali_path, code, target_name, tmpl_aln, target_aln,
              first, chain, last)
    with open(gaps_path, 'w') as fh:
        json.dump({'code': code, 'target_name': target_name,
                   'n_res': n_res, 'identity': round(ident, 2),
                   'gaps': gaps}, fh, indent=2)
    print('Wrote %s and %s' % (ali_path, gaps_path))


if __name__ == '__main__':
    main()
