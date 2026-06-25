"""Tests para mol_bio_toolkit.composition"""

from Bio.Seq import Seq
from mol_bio_toolkit import nucleotide_composition


def test_equal_composition():
    comp = nucleotide_composition(Seq("ATGC"))
    assert comp["A"] == 25.0
    assert comp["T"] == 25.0
    assert comp["G"] == 25.0
    assert comp["C"] == 25.0
    assert comp["GC%"] == 50.0


def test_all_at_zero_gc():
    comp = nucleotide_composition(Seq("ATATAT"))
    assert comp["GC%"] == 0.0
    assert comp["A"] == 50.0


def test_percentages_sum_to_100():
    comp = nucleotide_composition(Seq("ATGCATGCATGC"))
    total = comp["A"] + comp["T"] + comp["G"] + comp["C"]
    assert abs(total - 100.0) < 1e-6