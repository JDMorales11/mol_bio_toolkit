"""Tests for mol_bio_toolkit.sequence_search (NCBI mocked, no real network calls)."""
import io

import pytest

from mol_bio_toolkit import blast_protein, blast_nucleotide, identify_and_backtranslate

# Minimal but valid BLAST XML, with 2 hits:
# - Hit 1: low e_value (1e-05) -> should pass the default filter (0.001)
# - Hit 2: high e_value (0.5)  -> should be filtered out by default
# Query length is 10 so percentages come out as round numbers.
_SAMPLE_BLAST_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-len>10</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-len>10</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gb|ABC123.1|</Hit_id>
          <Hit_def>green fluorescent protein [Aequorea victoria]</Hit_def>
          <Hit_accession>ABC123</Hit_accession>
          <Hit_len>10</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>20.0</Hsp_bit-score>
              <Hsp_score>45</Hsp_score>
              <Hsp_evalue>1e-05</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>10</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>10</Hsp_hit-to>
              <Hsp_identity>9</Hsp_identity>
              <Hsp_positive>10</Hsp_positive>
              <Hsp_gaps>1</Hsp_gaps>
              <Hsp_align-len>10</Hsp_align-len>
              <Hsp_qseq>MKPGTEELFT</Hsp_qseq>
              <Hsp_hseq>MKPGTEELFT</Hsp_hseq>
              <Hsp_midline>MKPGTEELFT</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gb|XYZ999.1|</Hit_id>
          <Hit_def>hypothetical protein [Escherichia coli]</Hit_def>
          <Hit_accession>XYZ999</Hit_accession>
          <Hit_len>10</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>8.0</Hsp_bit-score>
              <Hsp_score>15</Hsp_score>
              <Hsp_evalue>0.5</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>5</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>5</Hsp_hit-to>
              <Hsp_identity>3</Hsp_identity>
              <Hsp_positive>4</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>5</Hsp_align-len>
              <Hsp_qseq>MKPGT</Hsp_qseq>
              <Hsp_hseq>MKAGT</Hsp_hseq>
              <Hsp_midline>MK GT</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def _mock_qblast(*args, **kwargs):
    """Replaces NCBIWWW.qblast: returns the fabricated XML instead of calling NCBI."""
    return io.StringIO(_SAMPLE_BLAST_XML)


def test_blast_protein_parses_hits_correctly():
    result = blast_protein("MKPGTEELFT", _qblast_fn=_mock_qblast)
    assert result["n_hits"] == 1  # the second hit (e_value=0.5) is filtered out by default
    hit = result["hits"].iloc[0]
    assert hit["accession"] == "ABC123"
    assert hit["organism"] == "Aequorea victoria"
    assert hit["pct_identity"] == 90.0   # 9/10 * 100
    assert hit["pct_similarity"] == 100.0  # 10/10 * 100
    assert hit["pct_coverage"] == 100.0  # (10-1+1)/10 * 100
    assert hit["pct_gaps"] == 10.0  # 1/10 * 100


def test_blast_protein_e_value_threshold_filters_hits():
    result = blast_protein(
        "MKPGTEELFT", _qblast_fn=_mock_qblast, e_value_threshold=1.0
    )
    assert result["n_hits"] == 2  # with a lax threshold, both hits pass


def test_blast_protein_empty_sequence_raises():
    with pytest.raises(ValueError):
        blast_protein("", _qblast_fn=_mock_qblast)


def test_blast_protein_hits_sorted_by_e_value():
    result = blast_protein(
        "MKPGTEELFT", _qblast_fn=_mock_qblast, e_value_threshold=1.0
    )
    e_values = result["hits"]["e_value"].tolist()
    assert e_values == sorted(e_values)


def test_blast_nucleotide_parses_hits_correctly():
    result = blast_nucleotide("ATGAAACCCGGGACTGAATTATTTACT", _qblast_fn=_mock_qblast)
    assert result["n_hits"] == 1
    assert result["database"] == "nt"


def test_blast_nucleotide_empty_sequence_raises():
    with pytest.raises(ValueError):
        blast_nucleotide("", _qblast_fn=_mock_qblast)


def test_organism_extraction_handles_missing_brackets():
    """If the title does not carry an organism in brackets, it must not raise."""
    xml_without_brackets = _SAMPLE_BLAST_XML.replace(
        "green fluorescent protein [Aequorea victoria]",
        "green fluorescent protein",
    )

    def mock_no_organism(*args, **kwargs):
        return io.StringIO(xml_without_brackets)

    result = blast_protein("MKPGTEELFT", _qblast_fn=mock_no_organism)
    assert result["hits"].iloc[0]["organism"] == "unknown"


# --- identify_and_backtranslate ------------------------------------------

def test_identify_and_backtranslate_returns_both_parts():
    protein = "MKPGTEELFT"
    result = identify_and_backtranslate(
        protein, host="e_coli_k12", _qblast_fn=_mock_qblast
    )
    assert result["protein"] == protein
    assert "backtranslation" in result
    assert "identification" in result
    assert result["backtranslation"]["variants"][0]["cds_seq"]
    assert result["identification"]["n_hits"] == 1


def test_identify_and_backtranslate_respects_host_and_mode():
    result = identify_and_backtranslate(
        "MKPG", host="s_cerevisiae", mode="markov", n_sequences=2, seed=3,
        _qblast_fn=_mock_qblast,
    )
    assert result["backtranslation"]["host"] == "s_cerevisiae"
    assert len(result["backtranslation"]["variants"]) == 2