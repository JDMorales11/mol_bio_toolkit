"""Tests para mol_bio_toolkit.primer_design"""

import pytest
from mol_bio_toolkit import (
    calculate_tm,
    gc_clamp,
    has_hairpin,
    find_optimal_primer,
    design_primers,
)


def test_tm_wallace_short_primer():
    # 4 GC, 0 AT, < 14 nt -> Wallace: 4*4 + 2*0 = 16
    tm = calculate_tm("GCGC")
    assert tm == 16.0


def test_tm_long_primer_returns_float():
    tm = calculate_tm("ATGCATGCATGCATGCATGC")  # 20 nt
    assert isinstance(tm, float)
    assert tm > 0


def test_gc_clamp_true_when_gc_in_tail():
    assert gc_clamp("AAAAATTTG") is True


def test_gc_clamp_false_when_only_at_in_tail():
    assert gc_clamp("GGGGGCCCAAT") is False


def test_hairpin_detected():
    assert has_hairpin("AAAAGGGGCCCCTTTT", min_stem=4) is True


def test_no_hairpin_in_simple_seq():
    assert has_hairpin("ATATATATATATATAT", min_stem=8) is False


def test_find_optimal_primer_returns_expected_keys():
    seq = "ATGC" * 30
    result = find_optimal_primer(seq, anchor=0, direction="forward")
    assert "Tm_C" in result
    assert "sequence" in result
    assert "all_candidates" in result


def test_design_primers_returns_pair():
    seq = "ATGC" * 50
    result = design_primers(seq, amplicon_start=0, amplicon_end=100)
    assert "forward" in result
    assert "reverse" in result
    assert result["amplicon_size"] == 100


def test_invalid_direction_raises():
    with pytest.raises(ValueError):
        find_optimal_primer("ATGC" * 10, anchor=0, direction="sideways")