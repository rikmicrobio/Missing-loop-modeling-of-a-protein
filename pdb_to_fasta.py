#!/usr/bin/env python
"""
Extract the full sequence of a chain from the SEQRES records of a PDB file.

SEQRES lists every residue of the construct, including those that were not
resolved (no coordinates). Compared with the residues that have coordinates,
the difference is exactly the missing loops/tails to be modelled.

Pure Python, no Modeller needed.

Usage
  python pdb_to_fasta.py --pdb 2LGV.pdb --chain A --out full_length.fasta
"""
import argparse
import os
import sys

AA3 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # common modified residues -> parent residue
    'MSE': 'M', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'CSO': 'C',
    'HYP': 'P', 'MLY': 'K', 'KCX': 'K', 'CME': 'C', 'SEC': 'C',
    'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'CYX': 'C', 'ASH': 'D',
    'GLH': 'E',
}


def seqres(path, chain):
    """List of 3-letter residue names for one chain, from SEQRES lines."""
    res = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('SEQRES') and len(line) > 19 and line[11] == chain:
                res.extend(line[19:].split())
    return res


def coordinate_sequence(path, chain):
    """Extract the observed protein sequence when a native PDB lacks SEQRES."""
    residues = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith('ATOM') or len(line) < 27 or line[21] != chain:
                continue
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            key = (line[22:27], line[17:20].strip())
            if key in seen:
                continue
            seen.add(key)
            residues.append(line[17:20].strip())
    return residues


def seqres_chains(path):
    chains = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('SEQRES') and len(line) > 11 \
                    and line[11] not in chains:
                chains.append(line[11])
    return chains


def to_one_letter(res3):
    seq, unknown = [], []
    for r in res3:
        if r in AA3:
            seq.append(AA3[r])
        else:
            seq.append('X')
            unknown.append(r)
    return ''.join(seq), unknown


def extract(pdb, chain):
    """Return (sequence, unknown residue names). Exits on problems."""
    if not os.path.exists(pdb):
        sys.exit('PDB file not found: %s' % pdb)
    chains = seqres_chains(pdb)
    if chains and chain in chains:
        residues = seqres(pdb, chain)
        source = 'SEQRES'
    else:
        # Many experimentally solved/native PDBs omit SEQRES records. In that
        # case use the observed C-alpha residue sequence so the pipeline can
        # still run without requiring an external FASTA file.
        residues = coordinate_sequence(pdb, chain)
        source = 'observed ATOM coordinates (no SEQRES found)'
        if not residues:
            available = []
            with open(pdb) as fh:
                for line in fh:
                    if line.startswith('ATOM') and len(line) > 21 and line[21] not in available:
                        available.append(line[21])
            sys.exit('Chain %r has no usable protein residues in %s. Available chains: %s'
                     % (chain, pdb, ', '.join(repr(c) for c in available)))
    seq, unknown = to_one_letter(residues)
    if 'X' in seq:
        sys.exit('Chain %s contains non-standard or non-protein residues '
                 '(%s), so a clean protein sequence cannot be extracted. '
                 'Provide the sequence with --target.'
                 % (chain, ', '.join(sorted(set(unknown)))))
    return seq, source


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--pdb', required=True)
    ap.add_argument('--chain', default='A')
    ap.add_argument('--out', help='output FASTA (default <code>_seqres.fasta)')
    args = ap.parse_args()

    code = os.path.splitext(os.path.basename(args.pdb))[0]
    out = args.out or code + '_seqres.fasta'
    seq, source = extract(args.pdb, args.chain)
    with open(out, 'w') as fh:
        fh.write('>%s_%s extracted_sequence\n' % (code, args.chain))
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + '\n')
    print('Wrote %s (chain %s, %d residues; source: %s)'
          % (out, args.chain, len(seq), source))


if __name__ == '__main__':
    main()
