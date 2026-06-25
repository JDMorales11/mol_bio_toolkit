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