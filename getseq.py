#!/usr/bin/env python
"""
Step 1 - extract the template sequence from a PDB file.

Writes <code>.seq containing only the residues that have coordinates in the
chosen chain. Gaps in the alignment against the full-length sequence (step 2)
are what Modeller will build later.

Usage
  python getseq.py --pdb 6I07.pdb --chain A
  python getseq.py --pdb 6I07.pdb            # first chain in the file
"""
import argparse
import os
import sys


def pdb_chains(path):
    """Chain IDs that have ATOM records, in file order (pure Python)."""
    chains = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('ATOM') and len(line) > 21:
                if line[21] not in chains:
                    chains.append(line[21])
    return chains


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--pdb', required=True, help='template PDB file')
    ap.add_argument('--chain', help='chain ID (default: first chain)')
    ap.add_argument('--outdir', default='.', help='where to write <code>.seq')
    args = ap.parse_args()

    if not os.path.exists(args.pdb):
        sys.exit('PDB file not found: %s' % args.pdb)
    code = os.path.splitext(os.path.basename(args.pdb))[0]

    chains = pdb_chains(args.pdb)
    if not chains:
        sys.exit('No ATOM records found in %s (mmCIF is not supported; '
                 'convert to PDB format).' % args.pdb)
    chain = args.chain if args.chain else chains[0]
    if chain not in chains:
        sys.exit('Chain %r not in %s. Available chains: %s'
                 % (chain, args.pdb, ', '.join(repr(c) for c in chains)))

    from modeller import Environ, Model, Alignment
    env = Environ()
    env.io.atom_files_directory = ['.', os.path.dirname(os.path.abspath(args.pdb))]
    seg = chain.strip()
    mdl = Model(env, file=code, model_segment=('FIRST:' + seg, 'LAST:' + seg))
    aln = Alignment(env)
    aln.append_model(mdl, align_codes=code)

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, code + '.seq')
    aln.write(file=out)
    print('Wrote %s (chain %s, %d residues with coordinates)'
          % (out, chain, len(aln[0].residues)))


if __name__ == '__main__':
    main()
