"""Homology search via NCBI BLAST and a combined identification +
back-translation workflow.

blast_protein and blast_nucleotide wrap Bio.Blast.NCBIWWW for remote
searches against NCBI databases. These calls are slow (typically
30-120 seconds) because they run on NCBI's servers, not locally, and
NCBI rate-limits requests from the same IP -- avoid looping searches
without a pause.
"""

import re
from typing import Callable, IO

from Bio.Blast import NCBIWWW, NCBIXML

import pandas as pd

from .codon_optimizer import back_translate

_ORGANISM_PATTERN = re.compile(r"\[([^\[\]]+)\]\s*$")


def _extract_organism(hit_description: str) -> str:
    """
    Extracts the organism name from an NCBI hit title.

    The typical NCBI format is "gene description [Organism]", with
    the organism in brackets at the end of the title.
    """
    match = _ORGANISM_PATTERN.search(hit_description)
    return match.group(1) if match else "unknown"


def _parse_blast_xml(xml_handle: IO[str], query_length: int) -> pd.DataFrame:
    """Internal helper: converts BLAST XML into a DataFrame of hits."""
    records = NCBIXML.parse(xml_handle)
    rows = []
    for record in records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                alignment_length = hsp.align_length
                coverage = (hsp.query_end - hsp.query_start + 1) / query_length * 100
                rows.append({
                    "accession": alignment.accession,
                    "description": alignment.hit_def,
                    "organism": _extract_organism(alignment.hit_def),
                    "pct_identity": round(hsp.identities / alignment_length * 100, 1),
                    "pct_similarity": round(hsp.positives / alignment_length * 100, 1),
                    "pct_coverage": round(coverage, 1),
                    "e_value": hsp.expect,
                    "bit_score": hsp.bits,
                    "alignment_length": alignment_length,
                    "pct_gaps": round(hsp.gaps / alignment_length * 100, 1) if hsp.gaps else 0.0,
                })
    return pd.DataFrame(rows)


def blast_protein(
    protein_sequence: str,
    database: str = "nr",
    hitlist_size: int = 10,
    e_value_threshold: float = 0.001,
    _qblast_fn: Callable = NCBIWWW.qblast,
) -> dict:
    """
    Searches for homologous protein sequences against NCBI databases.

    Sends the sequence to NCBI's servers with blastp and downloads the
    resulting alignments. Useful for identifying the likely gene and
    organism a protein sequence comes from (e.g. a designed protein,
    or one translated from an ORF detected with orf_finder whose
    identity is unknown).

    Parameters
    ----------
    protein_sequence : str
        Amino acid sequence (one-letter code).
    database : str, default "nr"
        NCBI database to query (e.g. "nr", "swissprot").
    hitlist_size : int, default 10
        Maximum number of hits to return.
    e_value_threshold : float, default 0.001
        E-value threshold; hits with a higher E-value are discarded.
    _qblast_fn : callable, default NCBIWWW.qblast
        Injection point for testing (allows replacing the real NCBI
        call with a mock in tests).

    Returns
    -------
    dict
        query, database, n_hits, hits (DataFrame with columns:
        accession, description, organism, pct_identity,
        pct_similarity, pct_coverage, e_value, bit_score,
        alignment_length, pct_gaps), filtered and sorted by E-value.

    Raises
    ------
    ValueError
        If protein_sequence is empty.

    Notes
    -----
    This function depends on the availability of NCBI's remote
    service and can take 30-120 seconds. It has no automatic retries:
    if NCBI is down or there are network issues, the exception
    propagates as raised by Biopython/urllib.
    """
    if not protein_sequence.strip():
        raise ValueError("protein_sequence cannot be empty")

    seq = protein_sequence.strip().upper()
    result_handle = _qblast_fn("blastp", database, seq, hitlist_size=hitlist_size)
    hits = _parse_blast_xml(result_handle, query_length=len(seq))

    if not hits.empty:
        hits = hits[hits["e_value"] <= e_value_threshold].sort_values("e_value")
        hits = hits.reset_index(drop=True)

    return {
        "query": seq,
        "database": database,
        "n_hits": len(hits),
        "hits": hits,
    }


def blast_nucleotide(
    nucleotide_sequence: str,
    database: str = "nt",
    hitlist_size: int = 10,
    e_value_threshold: float = 0.001,
    _qblast_fn: Callable = NCBIWWW.qblast,
) -> dict:
    """
    Searches for homologous DNA sequences against NCBI databases.

    Same logic as blast_protein but uses blastn on nucleotides.
    Useful, for example, to verify whether an optimized CDS is still
    recognizable as the original gene, or to confirm the identity of
    a sequence of unknown origin.

    Parameters
    ----------
    nucleotide_sequence : str
        DNA sequence (A/T/G/C).
    database : str, default "nt"
    hitlist_size : int, default 10
    e_value_threshold : float, default 0.001
    _qblast_fn : callable, default NCBIWWW.qblast
        Injection point for testing.

    Returns
    -------
    dict
        query, database, n_hits, hits (same format as blast_protein).

    Raises
    ------
    ValueError
        If nucleotide_sequence is empty.
    """
    if not nucleotide_sequence.strip():
        raise ValueError("nucleotide_sequence cannot be empty")

    seq = nucleotide_sequence.strip().upper()
    result_handle = _qblast_fn("blastn", database, seq, hitlist_size=hitlist_size)
    hits = _parse_blast_xml(result_handle, query_length=len(seq))

    if not hits.empty:
        hits = hits[hits["e_value"] <= e_value_threshold].sort_values("e_value")
        hits = hits.reset_index(drop=True)

    return {
        "query": seq,
        "database": database,
        "n_hits": len(hits),
        "hits": hits,
    }


def identify_and_backtranslate(
    protein_sequence: str,
    host: str = "e_coli_k12",
    mode: str = "max",
    n_sequences: int = 1,
    hitlist_size: int = 5,
    e_value_threshold: float = 0.001,
    seed: int | None = None,
    _qblast_fn: Callable = NCBIWWW.qblast,
) -> dict:
    """
    Combined workflow: from a protein, generates DNA candidates
    (back-translation) and identifies its likely origin (BLAST).

    This is the "reverse" functionality originally requested:
    starting only from an amino acid sequence, two things are
    obtained in a single call -- (1) one or more candidate nucleotide
    sequences optimized for host, and (2) the organism and gene this
    protein most likely comes from, via comparison against NCBI
    databases.

    It is literally back_translate() + blast_protein() over the same
    input, bundled together because they conceptually answer the same
    question: "I have this protein, give me usable nucleotides and
    tell me what it is."

    Parameters
    ----------
    protein_sequence : str
        Amino acid sequence (one-letter code).
    host : str, default "e_coli_k12"
        Target host for the back-translation.
    mode : {"max", "markov"}, default "max"
        Same codon selection logic as back_translate.
    n_sequences : int, default 1
        Number of nucleotide variants to generate.
    hitlist_size : int, default 5
        Maximum number of BLAST hits to return.
    e_value_threshold : float, default 0.001
    seed : int, optional
        Seed for the back-translation (see back_translate).
    _qblast_fn : callable, default NCBIWWW.qblast
        Injection point for testing.

    Returns
    -------
    dict
        protein, backtranslation (full back_translate output),
        identification (full blast_protein output).

    Raises
    ------
    ValueError
        If protein_sequence is empty or contains invalid characters.

    Notes
    -----
    BLAST identification can take 30-120 seconds since it is a remote
    NCBI call; back-translation is instantaneous since it runs
    locally. If you only need one of the two, use back_translate() or
    blast_protein() directly instead of this combined function.
    """
    backtranslation = back_translate(
        protein_sequence, host=host, mode=mode, n_sequences=n_sequences, seed=seed
    )
    identification = blast_protein(
        protein_sequence,
        hitlist_size=hitlist_size,
        e_value_threshold=e_value_threshold,
        _qblast_fn=_qblast_fn,
    )

    return {
        "protein": backtranslation["protein"],
        "backtranslation": backtranslation,
        "identification": identification,
    }