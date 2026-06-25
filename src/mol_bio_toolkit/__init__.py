"""
mol_bio_toolkit
===============

Herramientas para analisis y diseno de secuencias de DNA orientadas a
biologia sintetica: deteccion de ORFs, composicion nucleotidica,
diseno de primers y optimizacion de codones.
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
from .codon_optimizer import codon_optimize, CODON_TABLES

__version__ = "0.1.0"

__all__ = [
    "find_orfs",
    "nucleotide_composition",
    "calculate_tm",
    "gc_clamp",
    "has_hairpin",
    "find_optimal_primer",
    "design_primers",
    "codon_optimize",
    "CODON_TABLES",
]