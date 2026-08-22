"""
mol_bio_toolkit
===============

Tools for the analysis and design of DNA sequences for synthetic
biology: ORF detection, nucleotide composition, primer design, codon
optimization (deterministic and probabilistic via order-1 Markov),
back-translation, and homology search against NCBI BLAST.
"""

from .orf_finder import find_orfs
from .composition import nucleotide_composition
from .primer_design import (
    calculate_tm,
    gc_clamp,
    has_hairpin,
    find_optimal_primer,
    design_primers,
)
from .codon_optimizer import (
    codon_optimize,
    back_translate,
    CODON_TABLES,
    CODON_FREQUENCIES,
)
from .sequence_search import (
    blast_protein,
    blast_nucleotide,
    identify_and_backtranslate,
)

__version__ = "0.2.0"

__all__ = [
    "find_orfs",
    "nucleotide_composition",
    "calculate_tm",
    "gc_clamp",
    "has_hairpin",
    "find_optimal_primer",
    "design_primers",
    "codon_optimize",
    "back_translate",
    "CODON_TABLES",
    "CODON_FREQUENCIES",
    "blast_protein",
    "blast_nucleotide",
    "identify_and_backtranslate",
]