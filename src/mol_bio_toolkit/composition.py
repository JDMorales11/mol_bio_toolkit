"""Analisis de composicion de nucleotidos en secuencias de DNA."""

from Bio.SeqUtils import gc_fraction


def nucleotide_composition(sequence):
    """Calcula la composicion completa de una secuencia de DNA."""
    seq_str = str(sequence).upper()
    total = len(seq_str)
    return {
        "A": seq_str.count("A") / total * 100,
        "T": seq_str.count("T") / total * 100,
        "G": seq_str.count("G") / total * 100,
        "C": seq_str.count("C") / total * 100,
        "GC%": gc_fraction(sequence) * 100,
    }