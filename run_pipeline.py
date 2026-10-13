#!/usr/bin/env python
"""
MODELLER pipeline: template PDB + full-length sequence -> ranked models.

  getseq.py          template sequence from the PDB          -> <code>.seq
  pdb_to_fasta.py    full-length sequence from the PDB SEQRES -> <code>_seqres.fasta
                     (only when --target is not given)
  make_alignment.py  Clustal Omega alignment, find the gaps   -> .ali + gaps.json
  model-single.py    build models (+ optional loop refinement) -> ranked PDBs

Everything happens inside one work folder (default <code>_run), with a log
per step in <workdir>/logs/. Your input files are never modified.

Usage
  # full-length sequence taken from the PDB's own SEQRES records:
  python run_pipeline.py --pdb 6I07.pdb --chain A

  # or give your own full-length sequence (e.g. from UniProt):
  python run_pipeline.py --pdb 6I07.pdb --chain A --target full_length.fasta
  python run_pipeline.py --pdb 6I07.pdb --chain A --target full_length.fasta \
                         --loop --models 5 --loop-models 20

  # redo only the modelling after editing options (earlier outputs reused)
  python run_pipeline.py --pdb 6I07.pdb --chain A --target full_length.fasta \
                         --start-at model --loop

Run with the conda Python where Modeller is installed.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from getseq import pdb_chains          # noqa: E402  (pure Python helper)

STEPS = ['getseq', 'align', 'model']
RCSB_URL = 'https://files.rcsb.org/download/{}.pdb'
ID_RE = re.compile(r'^[0-9][A-Za-z0-9]{3}$')


def download_pdb(pdb_id, dest):
    url = RCSB_URL.format(pdb_id.upper())
    print('  downloading %s' % url)
    try:
        with urllib.request.urlopen(url, timeout=60) as r, open(dest, 'wb') as fh:
            fh.write(r.read())
    except Exception as exc:
        fail('could not download PDB entry %s (%s). Check the ID and your '
             'internet connection, or download it manually.' % (pdb_id, exc))
    if os.path.getsize(dest) < 200:
        fail('downloaded file for %s is empty. Check the ID.' % pdb_id)


def write_first_model(src, dst):
    """Copy a PDB file keeping only the first model (NMR ensembles).
    Returns the number of models found in the source."""
    n_models, done = 0, False
    with open(src) as fin, open(dst, 'w') as fout:
        for line in fin:
            if line.startswith('MODEL'):
                n_models += 1
            if not done:
                fout.write(line)
                if line.startswith('ENDMDL'):
                    done = True
    return n_models


def ensure_clustalo(install):
    if shutil.which('clustalo'):
        return
    if install and shutil.which('conda'):
        print('  clustalo not found, installing with conda ...')
        subprocess.run(['conda', 'install', '-y', '-c', 'bioconda',
                        '-c', 'conda-forge', 'clustalo'])
        if shutil.which('clustalo'):
            return
    fail('clustalo not found. Install it with:\n'
         '  conda install -c bioconda clustalo\n'
         'or re-run with --install-missing, or give --aligned.')


def fail(msg):
    sys.exit('\nERROR: ' + msg)


def run_step(name, cmd, workdir, verbose):
    """Run one step, log to logs/<name>.log, stop the pipeline on failure."""
    logdir = os.path.join(workdir, 'logs')
    os.makedirs(logdir, exist_ok=True)
    log = os.path.join(logdir, name + '.log')
    env = dict(os.environ, PYTHONUNBUFFERED='1')
    lines = []
    with open(log, 'w') as lf:
        p = subprocess.Popen(cmd, cwd=workdir, env=env, text=True, bufsize=1,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout:
            lf.write(line)
            lines.append(line)
            if verbose:
                sys.stdout.write(line)
        rc = p.wait()
    if rc != 0:
        if not verbose:
            sys.stdout.write(''.join(lines[-25:]))
        fail('step "%s" failed (exit code %d). Full log: %s' % (name, rc, log))
    if not verbose:
        # show the short, useful lines of the lightweight steps only
        if name != 'model':
            sys.stdout.write(''.join(lines))
    return log


def need(workdir, files, step):
    missing = [f for f in files if not os.path.exists(os.path.join(workdir, f))]
    if missing:
        fail('step "%s" did not produce: %s' % (step, ', '.join(missing)))


def preflight(args):
    if not os.path.exists(args.pdb):
        fail('PDB file not found: %s' % args.pdb)
    if args.target and not os.path.exists(args.target):
        fail('target FASTA not found: %s' % args.target)
    if args.aligned and not os.path.exists(args.aligned):
        fail('aligned FASTA not found: %s' % args.aligned)
    chains = pdb_chains(args.pdb)
    if not chains:
        fail('no ATOM records in %s (PDB format required, not mmCIF).'
             % args.pdb)
    if args.chain not in chains:
        fail('chain %r not in %s. Available: %s'
             % (args.chain, args.pdb, ', '.join(repr(c) for c in chains)))
    if args.start_at in ('getseq', 'align') and not args.aligned:
        ensure_clustalo(args.install_missing)
    if not args.skip_preflight:
        r = subprocess.run([sys.executable, '-c',
                            'from modeller import Environ; Environ()'],
                           capture_output=True, text=True)
        if r.returncode != 0:
            fail('Modeller could not start with %s.\n%s\nCheck that you are '
                 'in the right conda env and the license key is set in '
                 'modeller/config.py.' % (sys.executable,
                                          (r.stderr or r.stdout)[-600:]))


def show_results(workdir, target):
    tsv = os.path.join(workdir, target + '_ranking.tsv')
    print('\n=== Results ===')
    with open(tsv) as fh:
        rows = [l.rstrip('\n').split('\t') for l in fh][1:]
    print('%-42s %12s %8s' % ('model', 'DOPE', 'GA341'))
    for name, dope, ga, _ in rows:
        print('%-42s %12s %8s' % (name, dope, ga))
    print('\nBest model : %s' % os.path.join(workdir, target + '_best.pdb'))
    print('Full table : %s' % tsv)
    print('Logs       : %s' % os.path.join(workdir, 'logs'))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--pdb',
                    help='native/template PDB file; if omitted, a file picker opens')
    ap.add_argument('--chain', default=None, help='template chain (default: first protein chain)')
    ap.add_argument('--target',
                    help='FASTA with the full-length target sequence '
                         '(default: extracted from the SEQRES of the PDB)')
    ap.add_argument('--aligned', help='aligned FASTA (template, target); '
                                      'skips Clustal Omega')
    ap.add_argument('--workdir', help='work folder (default <code>_run)')
    ap.add_argument('--start-at', choices=STEPS, default='getseq',
                    help='resume from this step; earlier outputs must exist')
    ap.add_argument('--models', type=int, default=5,
                    help='models to build (default 5)')
    ap.add_argument('--loop', action='store_true',
                    help='refine the missing regions with LoopModel')
    ap.add_argument('--loop-models', type=int, default=10,
                    help='loop refinements per model (default 10)')
    ap.add_argument('--md-level', default='very_fast',
                    choices=['very_fast', 'fast', 'slow', 'very_slow'])
    ap.add_argument('--flank', type=int, default=2)
    ap.add_argument('--min-loop', type=int, default=1)
    ap.add_argument('--include-termini', action='store_true',
                    help='also refine missing N-/C-terminal tails')
    ap.add_argument('--soap', action='store_true',
                    help='also compute SOAP-Protein-OD (needs its library)')
    ap.add_argument('--verbose', action='store_true',
                    help='stream every step\'s full output to the screen')
    ap.add_argument('--install-missing', action='store_true',
                    help='install clustalo with conda if it is missing')
    ap.add_argument('--skip-preflight', action='store_true',
                    help='do not test-start Modeller first')
    args = ap.parse_args()

    # Branded entry point: custom green dot-matrix ASCII logo for Modus-Prime-RGLAB.
    # ANSI color is supported by modern terminals; fall back gracefully if needed.
    try:
        os.system('')  # Enables ANSI processing in many Windows terminals.
    except Exception:
        pass
    green = '\033[92m'
    reset = '\033[0m'
    banner = r'''•   •  •••  ••••  •   •  ••••      ••••  ••••  ••• •   • •••••      ••••   •••  •      •••  ••••
•• •• •   • •   • •   • •          •   • •   •  •  •• •• •          •   • •     •     •   • •   •
• • • •   • •   • •   •  •••  •••• ••••  ••••   •  • • • ••••  •••• ••••  •  •• •     ••••• ••••
•   • •   • •   • •   •     •      •     •  •   •  •   • •          •  •  •   • •     •   • •   •
•   •  •••  ••••   •••  ••••       •     •   • ••• •   • •••••      •   •  •••  ••••• •   • ••••'''
    print('\n' + green + banner + reset)
    print('             Native PDB Modelling Pipeline')
    print('=' * 78)
    if not args.pdb:
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            root.attributes('-topmost', True)
            selected = filedialog.askopenfilename(
                title='Modus-Prime-RGLAB | Select a native PDB file',
                filetypes=[('PDB structure files', '*.pdb *.ent'),
                           ('All files', '*.*')])
            root.destroy()
            if not selected:
                fail('No PDB file selected. Run again and select a native PDB file.')
            args.pdb = selected
        except Exception:
            args.pdb = input('Enter the path to your native PDB file: ').strip().strip('\"')
            if not args.pdb:
                fail('No PDB file provided.')

    is_id = (not os.path.exists(args.pdb)) and bool(ID_RE.match(args.pdb))
    code = args.pdb.upper() if is_id else \
        os.path.splitext(os.path.basename(args.pdb))[0]
    target = code + '_modelled'
    workdir = os.path.abspath(args.workdir or code + '_run')

    print('\nModus-Prime-RGLAB | MODELLER pipeline')
    print('  native PDB : %s' % args.pdb)
    print('  chain      : %s' % (args.chain or 'first available protein chain'))
    print('  target   : %s' % (args.target or 'SEQRES records of the PDB'))
    print('  workdir  : %s' % workdir)
    print('  mode     : %s, %d models%s' % (
        'AutoModel + LoopModel' if args.loop else 'AutoModel', args.models,
        ', %d loop models each' % args.loop_models if args.loop else ''))

    # stage inputs into the work folder (originals untouched)
    os.makedirs(workdir, exist_ok=True)
    pdb_local = code + '.pdb'
    pdb_path = os.path.join(workdir, pdb_local)
    if is_id:
        raw = os.path.join(workdir, code + '_download.pdb')
        download_pdb(args.pdb, raw)
        src = raw
    elif os.path.exists(args.pdb):
        src = args.pdb
    else:
        fail('PDB file not found: %s (use a file path or a 4-letter PDB ID)'
             % args.pdb)
    n_models = write_first_model(src, pdb_path)
    if n_models > 1:
        print('  %d models in the file: keeping model 1 only' % n_models)
    args.pdb = pdb_path
    available_chains = pdb_chains(args.pdb)
    if not args.chain:
        args.chain = available_chains[0] if available_chains else 'A'
        print('  automatically selected chain: %s' % args.chain)
    preflight(args)
    if args.target:
        target_local = os.path.basename(args.target)
        shutil.copy(args.target, os.path.join(workdir, target_local))
    else:
        target_local = code + '_seqres.fasta'
    aligned_local = None
    if args.aligned:
        aligned_local = os.path.basename(args.aligned)
        shutil.copy(args.aligned, os.path.join(workdir, aligned_local))

    py = sys.executable
    ali, gaps, seq = code + '_alignment.ali', code + '_gaps.json', code + '.seq'
    start = STEPS.index(args.start_at)
    t0 = time.time()

    # ---- step 1 ----
    if start <= 0:
        print('\n[1/3] getseq: template sequence')
        run_step('getseq', [py, os.path.join(HERE, 'getseq.py'),
                            '--pdb', pdb_local, '--chain', args.chain],
                 workdir, args.verbose)
    need(workdir, [seq], 'getseq')

    # ---- step 2 ----
    if start <= 1:
        if not args.target:
            print('\n[2/3] pdb_to_fasta: full-length sequence from SEQRES')
            run_step('fasta', [py, os.path.join(HERE, 'pdb_to_fasta.py'),
                               '--pdb', pdb_local, '--chain', args.chain,
                               '--out', target_local],
                     workdir, args.verbose)
            need(workdir, [target_local], 'fasta')
        print('\n[2/3] make_alignment: align and find missing regions')
        cmd = [py, os.path.join(HERE, 'make_alignment.py'),
               '--seq', seq, '--target', target_local,
               '--target-name', target]
        if aligned_local:
            cmd += ['--aligned', aligned_local]
        run_step('align', cmd, workdir, args.verbose)
    need(workdir, [ali, gaps], 'align')

    # ---- step 3 ----
    print('\n[3/3] model-single: building models (this can take a while; '
          'see %s)' % os.path.join(workdir, 'logs', 'model.log'))
    cmd = [py, os.path.join(HERE, 'model-single.py'),
           '--ali', ali, '--knowns', code, '--sequence', target,
           '--models', str(args.models)]
    if args.loop:
        cmd += ['--loop', '--gaps', gaps,
                '--loop-models', str(args.loop_models),
                '--md-level', args.md_level, '--flank', str(args.flank),
                '--min-loop', str(args.min_loop)]
        if args.include_termini:
            cmd.append('--include-termini')
    if args.soap:
        cmd.append('--soap')
    run_step('model', cmd, workdir, args.verbose)
    need(workdir, [target + '_ranking.tsv', target + '_best.pdb'], 'model')

    show_results(workdir, target)
    print('\nDone in %.1f min.' % ((time.time() - t0) / 60.0))


if __name__ == '__main__':
    main()
