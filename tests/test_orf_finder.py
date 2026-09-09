"""Tests para mol_bio_toolkit.orf_finder"""

from Bio.Seq import Seq
from mol_bio_toolkit import find_orfs


def test_finds_single_simple_orf():
    seq = Seq("ATGAAACCCGGGTAA")
    orfs = find_orfs(seq, min_length=10)
    assert len(orfs) == 1
    assert orfs[0]["frame"] == 1
    assert orfs[0]["protein"] == "MKPG"


def test_respects_min_length_filter():
    seq = Seq("ATGAAATAA")  
    orfs = find_orfs(seq, min_length=100)
    assert orfs == []


def test_no_atg_returns_empty():
    seq = Seq("CCCGGGTTTAAA")
    orfs = find_orfs(seq, min_length=5)
    assert orfs == []


def test_stops_at_first_stop_codon():
    seq = Seq("ATGAAATAACCCTAG")
    orfs = find_orfs(seq, min_length=1)
    assert len(orfs) == 1
    assert orfs[0]["sequence"] == "ATGAAATAA"


def test_finds_orfs_in_multiple_frames():
    seq = Seq("A" + "ATGAAACCCGGGTAA")
    orfs = find_orfs(seq, min_length=5)
    frames_found = {o["frame"] for o in orfs}
    assert 2 in frames_found



def test_forward_orf_has_plus_strand():
    seq = Seq("ATGAAACCCGGGTAA")
    orfs = find_orfs(seq, min_length=10)
    assert orfs[0]["strand"] == "+"
    assert orfs[0]["start"] == 0
    assert orfs[0]["end"] == 15


def test_finds_orf_on_reverse_strand():
    fwd_orf = "ATGAAACCCGGGTAA"
    rc = str(Seq(fwd_orf).reverse_complement())
    orfs = find_orfs(rc, min_length=10)
    assert len(orfs) == 1
    assert orfs[0]["strand"] == "-"
    assert orfs[0]["frame"] == -1
    assert orfs[0]["sequence"] == fwd_orf
    assert orfs[0]["protein"] == "MKPG"


def test_reverse_orf_coordinates_map_to_original_sequence():
    fwd_orf = "ATGAAACCCGGGTAA"
    rc = str(Seq(fwd_orf).reverse_complement())
    padded = "CCC" + rc + "AAA"
    orfs = find_orfs(padded, min_length=10)
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf["strand"] == "-"
    assert orf["start"] == 3
    assert orf["end"] == 18
    assert str(Seq(padded[orf["start"]:orf["end"]]).reverse_complement()) == fwd_orf


def test_both_strands_detected_in_same_sequence():
    fwd_orf = "ATGAAACCCGGGTAA"
    rc = str(Seq(fwd_orf).reverse_complement())
    combined = fwd_orf + rc
    orfs = find_orfs(combined, min_length=10)
    strands = {o["strand"] for o in orfs}
    assert "+" in strands
    assert "-" in strands


def test_reverse_min_length_filter_applies():
    fwd_orf = "ATGAAATAA"  
    rc = str(Seq(fwd_orf).reverse_complement())
    orfs = find_orfs(rc, min_length=100)
    assert orfs == []