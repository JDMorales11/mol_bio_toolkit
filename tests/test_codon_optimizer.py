"""Tests para mol_bio_toolkit.codon_optimizer"""

import pytest
from Bio.Seq import Seq
from mol_bio_toolkit import codon_optimize, CODON_TABLES


def test_protein_identity_preserved():
    """La propiedad mas importante: optimizar nunca debe cambiar la proteina."""
    seq = "ATGAAACCCGGGTAA"
    result = codon_optimize(seq, host="e_coli_k12")
    original_protein = str(Seq(seq[:-3]).translate())
    optimized_protein = str(Seq(result["optimized_seq"][:-3]).translate())
    assert original_protein == optimized_protein


def test_invalid_host_raises_keyerror():
    with pytest.raises(KeyError):
        codon_optimize("ATGAAATAA", host="martian_bacteria")


def test_all_registered_hosts_work():
    seq = "ATGAAACCCGGGTAA"
    for host in CODON_TABLES:
        result = codon_optimize(seq, host=host)
        assert result["host"] == host
        assert len(result["optimized_seq"]) == len(seq)


def test_gc_content_keys_present():
    result = codon_optimize("ATGAAACCCGGGTAA", host="e_coli_k12")
    assert "original_gc%" in result
    assert "optimized_gc%" in result
    assert 0 <= result["original_gc%"] <= 100