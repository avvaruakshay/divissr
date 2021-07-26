#! /usr/bin/env python
# pylint: disable=C0111, C0301

# Future
from __future__ import print_function, division

# Generic/Built-in
import gzip, sys
from itertools import takewhile, repeat

# list storing kmer names for each motif size in sorted order
kmer_names = [
    'Monomer', 'Dimer', 'Trimer', 'Tetramer', 'Pentamer',
    'Hexamer', 'Heptamer', 'Octamer', 'Nonamer', 'Decamer',
    'Undecamer', 'Dodecamer', 'Tridecamer', 'Tetradecamer',
    'Pentadecamer', 'Hexadecamer', 'Heptadecamer', 'Octadecamer',
    'Nonadecamer', 'Icosamer', 'Uncosamer', 'Docosamer',
    'Tricosamer', 'Tetracosamer', 'Pentacosamer', 'Hexacosamer',
    'Heptacosamer', 'Octacosamer', 'Nonacosamer', 'Triacontamer',
    'Untriacontamer', 'Dotriacontamer', 'Tritriacontamer',
    'Tetratriacontamer', 'Pentatriacontamer', 'Hexatriacontamer',
    'Heptatriacontamer', 'Octatriacontamer', 'Nonatriacontamer',
    'Tetracontamer', 'Untetracontamer', 'Dotetracontamer',
    'Tritetracontamer', 'Tetratetracontamer', 'Pentatetracontamer',
    'Hexatetracontamer', 'Heptatetracontamer', 'Octatetracontamer', 
    'Nonatetracontamer', 'Pentacontamer',
]


def rawcharCount(filename, char):
    """
    Counts the occurrences of a character in the file.

    Parameters
    ----------
    filename : str, input file name
    char : str, character to be counted

    Returns
    -------
    INT, count of the character
    """
    if filename.endswith('gz'):
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(char.encode('ASCII')) for buf in bufgen if buf )


def get_cycles(motif):
    """
    Generates cyclical variations of a motif
    
    Parameters
    ----------
    motif : str, nucleotide motif

    Returns
    -------
    LIST, List of cyclical variations of the motif sorted in alphabetical order.
    """
    cycles = set()
    for i in range(len(motif)):
        cycles.add(motif[i:] + motif[:i])
    cycles = sorted(list(cycles))
    return cycles


def build_cycVariations(motif):
    """
    Description
    -----------
    Builds the set of cyclical variations of the motif and cyclical variations
    of the reverse complement of the motif.
    
    Parameters
    ----------
    motif : str, nucleotide motif

    Returns
    -------
    LIST,  
    Cyclical variations of itself and the reverse complement.
    The list has the cyclical variations of the motif sorted in 
    alphabetical order followed by alphabetically sorted cyclical 
    variations of reverse complement.
    """
    cycles = get_cycles(motif)
    rev_cycles = get_cycles(rev_comp(motif))
    for r in rev_cycles:
        if r not in cycles: cycles.append(r)
    return cycles


def getGC(basesCounter):
    """
    Description
    -----------
    Calculate GC percentage of the genome.

    Parameters
    ----------
    baseCounter : collections Counter object which has counts of each of the 
    nucleotide.

    Returns
    -------
    FLOAT, GC percentage value
    """
    totalBases = sum(basesCounter.values())
    try:
        GC = (float(basesCounter['G'] + basesCounter['C'])/(totalBases-basesCounter['N']))*100
    except KeyError:
        GC = (float(basesCounter['G'] + basesCounter['C'])/totalBases)*100
    return GC


def getGenomeInfo(filename):
    """
    Description
    -----------
    Calculate GC percentage of the genome.

    Parameters
    ----------
    filename : Name of the input file.

    Returns
    -------
    LIST, [Genome size, GC percentage value]
    """
    bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    if filename.endswith('gz'): fh = gzip.open(filename, 'rt')
    else: fh = open(filename, 'r')
    for line in fh:
        if not line.strip().startswith('>'):
            for nuc in bases: bases[nuc] += line.upper().count(nuc)
    gsize = sum(bases.values())
    try:
        GC = (float(bases['G'] + bases['C'])/(gsize-bases['N']))*100
    except KeyError:
        GC = (float(bases['G'] + bases['C'])/gsize)*100
    return [gsize, GC]


def rev_comp(motif):
    """
    Description
    -----------
    Outputs reverse complement of a nucleotide sequence

    Parameters
    ----------
    motif : str, input nucleotide sequence

    Returns
    -------
    STR, reverse complement sequence    
    """
    if sys.version_info.major == 2:
        import motif as st
        complement = motif.translate(st.maketrans('ACGT', 'TGCA'))
    else:
        complement = motif.translate(str.maketrans('ACGT', 'TGCA'))
    return complement[::-1]
