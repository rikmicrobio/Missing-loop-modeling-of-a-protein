#!/usr/bin/env python
"""
Step 3 - build the models (and optionally refine the missing loops).

  plain      AutoModel builds N models; missing residues are built de novo
  --loop     LoopModel additionally refines the missing regions listed in the
             gaps JSON from make_alignment.py (+ flanking residues)

Models are ranked by DOPE (lower = better). Outputs:
  <sequence>_ranking.tsv     all models with scores
  <sequence>_best.pdb        copy of the best model

Usage
  python model-single.py --ali 6I07_alignment.ali --knowns 6I07 \
         --sequence 6I07_modelled --gaps 6I07_gaps.json --models 5 --loop
"""
import argparse
import json
import os
import shutil
import sys

from modeller import Environ, Selection
from modeller.automodel import AutoModel, LoopModel, assess, refine

LOOPS = []     # filled in main(); read by MissingLoopModel


class MissingLoopModel(LoopModel):
    def select_loop_atoms(self):
        ch = self.chains[0].name
        return Selection(*[
            self.residue_range('%d:%s' % (s, ch), '%d:%s' % (e, ch))
            for s, e in LOOPS])


def pad_ranges(gaps, n_res, flank, min_len, include_termini):
    """Filter gaps, add flanking residues, clip to the chain, merge overlaps."""
    kept = []
    for s, e in gaps:
        if e - s + 1 < min_len:
            continue
        if not include_termini and (s == 1 or e == n_res):
            continue
        kept.append((max(1, s - flank), min(n_res, e + flank)))
    kept.sort()
    merged = []
    for s, e in kept:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def ga341(o):
    v = o.get('GA341 score')
    if isinstance(v, (list, tuple)):
        v = v[0] if v else None
    return v


def rank_and_save(outputs, prefix):
    ok = [o for o in outputs
          if o.get('failure') is None and 'DOPE score' in o]
    if not ok:
        sys.exit('No model finished successfully (or no DOPE scores). '
                 'Check the log.')
    ok.sort(key=lambda o: o['DOPE score'])
    tsv = prefix + '_ranking.tsv'
    with open(tsv, 'w') as fh:
        fh.write('model\tDOPE\tGA341\tmolpdf\n')
        for o in ok:
            ga = ga341(o)
            fh.write('%s\t%.2f\t%s\t%s\n' % (
                o['name'], o['DOPE score'],
                '%.4f' % ga if ga is not None else 'NA',
                '%.2f' % o['molpdf'] if 'molpdf' in o else 'NA'))
    best = prefix + '_best.pdb'
    shutil.copy(ok[0]['name'], best)
    print('\n=== Models ranked by DOPE (lower is better) ===')
    for o in ok:
        print('%-42s %12.2f' % (o['name'], o['DOPE score']))
    print('\nBest model: %s  ->  %s' % (ok[0]['name'], best))
    print('Ranking table: %s' % tsv)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--ali', required=True, help='alignment (.ali) file')
    ap.add_argument('--knowns', required=True, help='template code')
    ap.add_argument('--sequence', required=True, help='target code in the .ali')
    ap.add_argument('--models', type=int, default=5,
                    help='number of models (default 5)')
    ap.add_argument('--atom-dir', default='.', help='folder with the PDB file')
    ap.add_argument('--soap', action='store_true',
                    help='also compute SOAP-Protein-OD (needs its library)')
    ap.add_argument('--loop', action='store_true',
                    help='refine missing regions with LoopModel')
    ap.add_argument('--gaps', help='gaps JSON from make_alignment.py '
                                   '(needed with --loop)')
    ap.add_argument('--loop-models', type=int, default=10,
                    help='loop refinements per model (default 10)')
    ap.add_argument('--md-level', default='very_fast',
                    choices=['very_fast', 'fast', 'slow', 'very_slow'])
    ap.add_argument('--flank', type=int, default=2,
                    help='residues refined on each side of a gap (default 2)')
    ap.add_argument('--min-loop', type=int, default=1,
                    help='ignore gaps shorter than this (default 1)')
    ap.add_argument('--include-termini', action='store_true',
                    help='also refine missing N-/C-terminal tails')
    args = ap.parse_args()

    use_loop = False
    if args.loop:
        if not args.gaps:
            sys.exit('--loop needs --gaps <code>_gaps.json')
        with open(args.gaps) as fh:
            info = json.load(fh)
        LOOPS[:] = pad_ranges(info['gaps'], info['n_res'], args.flank,
                              args.min_loop, args.include_termini)
        if LOOPS:
            use_loop = True
            print('Regions to refine: %s'
                  % ', '.join('%d-%d' % l for l in LOOPS))
        else:
            print('No internal gaps to refine (use --include-termini for '
                  'terminal tails). Running plain AutoModel.')

    env = Environ()
    env.io.atom_files_directory = ['.', args.atom_dir]

    scorers = [assess.DOPE, assess.GA341]
    if args.soap:
        from modeller import soap_protein_od
        scorers.append(soap_protein_od.Scorer())
    common = dict(alnfile=args.ali, knowns=args.knowns,
                  sequence=args.sequence, assess_methods=tuple(scorers))

    if use_loop:
        a = MissingLoopModel(env, loop_assess_methods=tuple(scorers), **common)
        a.loop.starting_model = 1
        a.loop.ending_model = args.loop_models
        a.loop.md_level = getattr(refine, args.md_level)
    else:
        a = AutoModel(env, **common)
    a.starting_model = 1
    a.ending_model = args.models
    a.make()

    rank_and_save(a.loop.outputs if use_loop else a.outputs, args.sequence)


if __name__ == '__main__':
    main()
