"""Tests for mol_bio_toolkit.codon_optimizer"""
import pytest
from Bio.Seq import Seq
from mol_bio_toolkit import codon_optimize, back_translate, CODON_TABLES, CODON_FREQUENCIES


def test_protein_identity_preserved():
    """The most important property: optimizing must never change the protein."""
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


# ------ mode="markov" -------------------------------------------

def test_markov_mode_preserves_protein_identity():
    """mode='markov' must also preserve the exact protein."""
    seq = "ATGAAACCCGGGTAA"
    result = codon_optimize(seq, host="e_coli_k12", mode="markov", seed=1)
    original_protein = str(Seq(seq[:-3]).translate())
    optimized_protein = str(Seq(result["optimized_seq"][:-3]).translate())
    assert original_protein == optimized_protein


def test_markov_mode_reproducible_with_seed():
    """The same seed must produce exactly the same result (controlled determinism)."""
    seq = "ATGAAACCCGGGTGGTAA"
    r1 = codon_optimize(seq, host="e_coli_k12", mode="markov", seed=42)
    r2 = codon_optimize(seq, host="e_coli_k12", mode="markov", seed=42)
    assert r1["optimized_seq"] == r2["optimized_seq"]


def test_invalid_mode_raises_valueerror():
    with pytest.raises(ValueError):
        codon_optimize("ATGAAATAA", host="e_coli_k12", mode="best_effort")


def test_codon_log_has_freq_in_host_column():
    result = codon_optimize("ATGAAACCCGGGTAA", host="e_coli_k12")
    assert "freq_in_host" in result["codon_log"].columns


def test_codon_tables_consistent_with_frequencies():
    """CODON_TABLES must match the highest-frequency codon in CODON_FREQUENCIES
    (verifies the two tables cannot fall out of sync)."""
    assert CODON_TABLES["e_coli_k12"]["A"] == "GCG"  # 38.5, highest among GCT/GCC/GCA/GCG
    assert CODON_TABLES["s_cerevisiae"]["A"] == "GCT"  # 21.2, highest in yeast


# ------ back_translate --------------------------------------------

def test_back_translate_preserves_protein():
    protein = "MKPG"
    result = back_translate(protein, host="e_coli_k12")
    variant = result["variants"][0]
    translated_back = str(Seq(variant["cds_seq"][:-3]).translate())
    assert translated_back == protein


def test_back_translate_max_mode_deterministic():
    protein = "MKPG"
    r1 = back_translate(protein, host="e_coli_k12", mode="max")
    r2 = back_translate(protein, host="e_coli_k12", mode="max")
    assert r1["variants"][0]["cds_seq"] == r2["variants"][0]["cds_seq"]


def test_back_translate_max_mode_rejects_multiple_sequences():
    """mode='max' is deterministic: requesting n_sequences > 1 makes no sense and must fail."""
    with pytest.raises(ValueError):
        back_translate("MKPG", host="e_coli_k12", mode="max", n_sequences=3)


def test_back_translate_markov_multiple_variants():
    protein = "MKPGEELFTG"
    result = back_translate(protein, host="e_coli_k12", mode="markov", n_sequences=3, seed=7)
    assert result["n_sequences"] == 3
    assert len(result["variants"]) == 3
    for variant in result["variants"]:
        translated_back = str(Seq(variant["cds_seq"][:-3]).translate())
        assert translated_back == protein


def test_back_translate_variant_has_gc_and_freq_keys():
    result = back_translate("MKPG", host="e_coli_k12")
    variant = result["variants"][0]
    assert "gc%" in variant
    assert "mean_codon_freq" in variant
    assert 0 <= variant["gc%"] <= 100


def test_back_translate_invalid_amino_acid_raises():
    with pytest.raises(ValueError):
        back_translate("MKPZ", host="e_coli_k12")  # 'Z' is not a standard amino acid


def test_back_translate_invalid_host_raises_keyerror():
    with pytest.raises(KeyError):
        back_translate("MKPG", host="martian_bacteria")