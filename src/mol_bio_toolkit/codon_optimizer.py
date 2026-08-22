"""Codon optimization for heterologous expression.

Includes two codon selection strategies (mode="max" and
mode="markov"), probabilistic back-translation of a protein into
nucleotides (back_translate), and penalization of CpG/TpA
dinucleotides at codon junctions.

Frequency table source: Codon Usage Database (Kazusa DNA Research
Institute), based on Nakamura, Y., Gojobori, T. and Ikemura, T.
(2000) "Codon usage tabulated from the international DNA sequence
databases: status for the year 2000." Nucleic Acids Research 28, 292.
Tables downloaded from kazusa.or.jp/codon/ for:
  - Escherichia coli K12 (taxid 83333): 14 CDS's, 5122 codons
  - Saccharomyces cerevisiae (taxid 4932): 14411 CDS's, 6534504 codons
"""

import random

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

# ---------------------------------------------------------------------
# Standard genetic code: amino acid -> list of synonymous codons
# ---------------------------------------------------------------------
AMINO_ACID_CODONS = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "*": ["TAA", "TAG", "TGA"],
}

# codon -> amino acid (reverse lookup, useful for reading an existing CDS)
CODON_TO_AA = {
    codon: aa for aa, codons in AMINO_ACID_CODONS.items() for codon in codons
}

# ---------------------------------------------------------------------
# Kazusa frequency tables (codons per 1000 codons)
# Single source of truth: everything else ("max" codon, CODON_TABLES,
# sampling weights) is derived from here, so the two tables can never
# fall out of sync.
# ---------------------------------------------------------------------
CODON_FREQUENCIES = {
    "e_coli_k12": {
        "TTT": 19.7, "TTC": 15.0, "TTA": 15.2, "TTG": 11.9,
        "CTT": 11.9, "CTC": 10.5, "CTA": 5.3, "CTG": 46.9,
        "ATT": 30.5, "ATC": 18.2, "ATA": 3.7,
        "ATG": 24.8,
        "GTT": 16.8, "GTC": 11.7, "GTA": 11.5, "GTG": 26.4,
        "TCT": 5.7, "TCC": 5.5, "TCA": 7.8, "TCG": 8.0,
        "CCT": 8.4, "CCC": 6.4, "CCA": 6.6, "CCG": 26.7,
        "ACT": 8.0, "ACC": 22.8, "ACA": 6.4, "ACG": 11.5,
        "GCT": 10.7, "GCC": 31.6, "GCA": 21.1, "GCG": 38.5,
        "TAT": 16.8, "TAC": 14.6, "TAA": 1.8, "TAG": 0.0,
        "CAT": 15.8, "CAC": 13.1, "CAA": 12.1, "CAG": 27.7,
        "AAT": 21.9, "AAC": 24.4, "AAA": 33.2, "AAG": 12.1,
        "GAT": 37.9, "GAC": 20.5, "GAA": 43.7, "GAG": 18.4,
        "TGT": 5.9, "TGC": 8.0, "TGA": 1.0, "TGG": 10.7,
        "CGT": 21.1, "CGC": 26.0, "CGA": 4.3, "CGG": 4.1,
        "AGT": 7.2, "AGC": 16.6, "AGA": 1.4, "AGG": 1.6,
        "GGT": 21.3, "GGC": 33.4, "GGA": 9.2, "GGG": 8.6,
    },
    "s_cerevisiae": {
        "TTT": 26.1, "TTC": 18.4, "TTA": 26.2, "TTG": 27.2,
        "CTT": 12.3, "CTC": 5.4, "CTA": 13.4, "CTG": 10.5,
        "ATT": 30.1, "ATC": 17.2, "ATA": 17.8,
        "ATG": 20.9,
        "GTT": 22.1, "GTC": 11.8, "GTA": 11.8, "GTG": 10.8,
        "TCT": 23.5, "TCC": 14.2, "TCA": 18.7, "TCG": 8.6,
        "CCT": 13.5, "CCC": 6.8, "CCA": 18.3, "CCG": 5.3,
        "ACT": 20.3, "ACC": 12.7, "ACA": 17.8, "ACG": 8.0,
        "GCT": 21.2, "GCC": 12.6, "GCA": 16.2, "GCG": 6.2,
        "TAT": 18.8, "TAC": 14.8, "TAA": 1.1, "TAG": 0.5,
        "CAT": 13.6, "CAC": 7.8, "CAA": 27.3, "CAG": 12.1,
        "AAT": 35.7, "AAC": 24.8, "AAA": 41.9, "AAG": 30.8,
        "GAT": 37.6, "GAC": 20.2, "GAA": 45.6, "GAG": 19.2,
        "TGT": 8.1, "TGC": 4.8, "TGA": 0.7, "TGG": 10.4,
        "CGT": 6.4, "CGC": 2.6, "CGA": 3.0, "CGG": 1.7,
        "AGT": 14.2, "AGC": 9.8, "AGA": 21.3, "AGG": 9.2,
        "GGT": 23.9, "GGC": 9.8, "GGA": 10.9, "GGG": 6.0,
    },
}

_DINUCLEOTIDE_PENALTY = 0.15  # multiplicative factor, not a binary exclusion
_PENALIZED_DINUCLEOTIDES = {"CG", "TA"}


def _codon_weights(host: str, aa: str) -> dict[str, float]:
    """Kazusa frequencies (unnormalized) of the synonymous codons for `aa`."""
    table = CODON_FREQUENCIES[host]
    return {c: table[c] for c in AMINO_ACID_CODONS[aa]}


def _best_codon(host: str, aa: str) -> str:
    """Highest-frequency synonymous codon for `aa` in `host` per Kazusa."""
    weights = _codon_weights(host, aa)
    return max(weights, key=lambda c: weights[c])


# CODON_TABLES: mapping host -> {amino_acid: most frequent codon}.
# Kept for backward compatibility with v0.1.x (notebook, tests, and
# README already import it), but it is now DERIVED from
# CODON_FREQUENCIES instead of being maintained by hand, so the two
# tables can never fall out of sync.
CODON_TABLES = {
    host: {aa: _best_codon(host, aa) for aa in AMINO_ACID_CODONS}
    for host in CODON_FREQUENCIES
}


def _sample_codon(
    host: str, aa: str, prev_codon: str | None, rng: random.Random
) -> str:
    """
    Samples a synonymous codon for `aa`, weighted by its Kazusa
    frequency and penalized if it forms a CpG or TpA dinucleotide
    with the last base of `prev_codon` (order-1 dependence: the
    probability of each state depends on the previous state).
    """
    weights = dict(_codon_weights(host, aa))
    if prev_codon is not None:
        last_base = prev_codon[-1]
        for codon in weights:
            junction = last_base + codon[0]
            if junction in _PENALIZED_DINUCLEOTIDES:
                weights[codon] *= _DINUCLEOTIDE_PENALTY
    codons = list(weights.keys())
    scores = list(weights.values())
    return rng.choices(codons, weights=scores, k=1)[0]


def codon_optimize(
    cds_sequence: str | Seq,
    host: str = "e_coli_k12",
    mode: str = "max",
    seed: int | None = None,
) -> dict:
    """
    Recodes a CDS using the preferred codons of the host organism.

    The genetic code is degenerate: 61 codons encode only 20 amino
    acids, so most amino acids have multiple possible synonymous
    codons. Each organism uses those synonyms with different
    frequency (codon usage bias), reflecting the relative abundance
    of tRNAs in that cell. This function preserves the resulting
    protein exactly: it only changes codons, never the amino acid
    sequence.

    Two selection strategies:

    - mode="max" (deterministic, v0.1.x-compatible): at each position
      always uses the synonymous codon with the highest Kazusa
      frequency for that amino acid. Maximizes codon usage bias, but
      a gene that is 100% "max codon" is an unnatural pattern (no
      real gene uses only the top codon at every position) and can
      generate repetitive stretches that favor mRNA secondary
      structure or synthesis issues (homopolymers, repeated
      restriction sites).

    - mode="markov" (stochastic, order 1): at each position samples a
      synonymous codon with probability proportional to its Kazusa
      frequency, but the distribution at that position depends on the
      immediately preceding codon: if the last nucleotide of the
      previous codon plus the first nucleotide of the candidate form
      a CpG or TpA dinucleotide, its weight is reduced. Both
      dinucleotides are underrepresented in bacterial and yeast
      genomes due to mutational pressure (CpG is additionally a
      target for methylation and restriction-modification systems
      such as Dcm/Dam in E. coli). This mode produces sequences that
      are more natural for the host. Closes Issue #5.

    Methodological note: this is a first-order approximation based on
    marginal codon frequencies (Kazusa) plus a junction penalty, not a
    true codon-pair bias table (which would require frequencies of
    consecutive codon pairs measured directly from the genome; Kazusa
    does not provide these). Documented as a known limitation.

    Parameters
    ----------
    cds_sequence : Bio.Seq or str
        Full CDS in reading frame (divisible by 3, includes the stop
        codon).
    host : str, default "e_coli_k12"
        Target organism. Options: "e_coli_k12", "s_cerevisiae".
    mode : {"max", "markov"}, default "max"
    seed : int, optional
        Random number generator seed (only applies to mode="markov").

    Returns
    -------
    dict
        host, mode, original_seq, optimized_seq, aa_length,
        codon_changes, identity_%, original_gc%, optimized_gc%,
        codon_log (DataFrame with columns: position, aa, original,
        optimized, changed, freq_in_host).

    Raises
    ------
    KeyError
        If host is not in CODON_FREQUENCIES.
    ValueError
        If mode is neither "max" nor "markov".
    """
    if host not in CODON_FREQUENCIES:
        raise KeyError(
            f"Host '{host}' not recognized. Options: {list(CODON_FREQUENCIES)}"
        )
    if mode not in ("max", "markov"):
        raise ValueError("mode must be 'max' or 'markov'")

    freq_table = CODON_FREQUENCIES[host]
    seq_str = str(cds_sequence).upper()
    protein = str(Seq(seq_str[:-3]).translate())
    rng = random.Random(seed)

    optimized = ""
    changes = 0
    codon_log = []
    prev_codon = None

    for i, aa in enumerate(protein):
        original_codon = seq_str[i * 3: i * 3 + 3]
        if mode == "max":
            new_codon = _best_codon(host, aa)
        else:
            new_codon = _sample_codon(host, aa, prev_codon, rng)
        optimized += new_codon
        changed = original_codon != new_codon
        if changed:
            changes += 1
        codon_log.append({
            "position": i + 1,
            "aa": aa,
            "original": original_codon,
            "optimized": new_codon,
            "changed": changed,
            "freq_in_host": freq_table[new_codon],
        })
        prev_codon = new_codon

    stop_codon = (
        _best_codon(host, "*") if mode == "max"
        else _sample_codon(host, "*", prev_codon, rng)
    )
    optimized += stop_codon
    identity = (len(protein) + 1 - changes) / (len(protein) + 1) * 100

    return {
        "host": host,
        "mode": mode,
        "original_seq": seq_str,
        "optimized_seq": optimized,
        "aa_length": len(protein),
        "codon_changes": changes,
        "identity_%": round(identity, 1),
        "original_gc%": round(gc_fraction(Seq(seq_str)) * 100, 1),
        "optimized_gc%": round(gc_fraction(Seq(optimized)) * 100, 1),
        "codon_log": pd.DataFrame(codon_log),
    }


def _generate_one_variant(
    protein: str, host: str, mode: str, include_stop: bool, rng: random.Random
) -> dict:
    """Internal helper: generates a single back-translation variant."""
    freq_table = CODON_FREQUENCIES[host]
    cds = ""
    codon_log = []
    prev_codon = None

    for i, aa in enumerate(protein):
        codon = _best_codon(host, aa) if mode == "max" else _sample_codon(host, aa, prev_codon, rng)
        cds += codon
        codon_log.append({"position": i + 1, "aa": aa, "codon": codon, "freq_in_host": freq_table[codon]})
        prev_codon = codon

    if include_stop:
        stop_codon = (
            _best_codon(host, "*") if mode == "max"
            else _sample_codon(host, "*", prev_codon, rng)
        )
        cds += stop_codon
        codon_log.append({
            "position": len(protein) + 1, "aa": "*", "codon": stop_codon,
            "freq_in_host": freq_table[stop_codon],
        })

    log_df = pd.DataFrame(codon_log)
    return {
        "cds_seq": cds,
        "gc%": round(gc_fraction(Seq(cds)) * 100, 1),
        "mean_codon_freq": round(log_df["freq_in_host"].mean(), 2),
        "codon_log": log_df,
    }


def back_translate(
    protein_sequence: str,
    host: str = "e_coli_k12",
    mode: str = "max",
    n_sequences: int = 1,
    include_stop: bool = True,
    seed: int | None = None,
) -> dict:
    """
    Translates a protein sequence back into one or more DNA CDSs.

    Unlike codon_optimize (which recodes an already-existing DNA
    CDS), back_translate starts only from the amino acid sequence:
    there are no "original" codons to preserve, each position is
    generated from scratch according to the host's frequency table.
    This is the relevant tool when a protein has been designed (e.g.
    a synthetic sequence, or one taken from another species) and a
    candidate gene is needed for synthesis. Closes Issue #4 partially.

    With mode="markov" and n_sequences > 1, each variant is an
    independent probabilistic sample (same order-1 logic as
    codon_optimize): the function does not pick "the best one" for
    you, it returns the n variants with their GC% and mean codon
    frequency so you can compare and choose based on your own
    criteria (e.g. GC% within a synthesis-friendly range, or the
    variant with the highest mean frequency if expression is the
    priority).

    Parameters
    ----------
    protein_sequence : str
        Amino acid sequence (one-letter code, without a trailing "*").
    host : str, default "e_coli_k12"
    mode : {"max", "markov"}, default "max"
        With mode="max" the result is deterministic, so
        n_sequences > 1 would not produce distinct variants and is
        rejected.
    n_sequences : int, default 1
        Number of variants to generate. Only meaningful > 1 with
        mode="markov".
    include_stop : bool, default True
    seed : int, optional
        Base seed; each variant uses seed + index so it is
        individually reproducible and distinct from the others.

    Returns
    -------
    dict
        host, mode, protein, aa_length, n_sequences,
        variants: list of dicts, each with cds_seq, gc%,
        mean_codon_freq, codon_log (DataFrame).

    Raises
    ------
    KeyError
        If host is not in CODON_FREQUENCIES.
    ValueError
        If mode is invalid, if n_sequences < 1, if mode="max" and
        n_sequences > 1, or if the protein contains a character that
        does not correspond to a standard amino acid.
    """
    if host not in CODON_FREQUENCIES:
        raise KeyError(
            f"Host '{host}' not recognized. Options: {list(CODON_FREQUENCIES)}"
        )
    if mode not in ("max", "markov"):
        raise ValueError("mode must be 'max' or 'markov'")
    if n_sequences < 1:
        raise ValueError("n_sequences must be >= 1")
    if mode == "max" and n_sequences > 1:
        raise ValueError(
            "mode='max' is deterministic: n_sequences > 1 would "
            "generate identical copies. Use mode='markov' for variants."
        )

    protein = protein_sequence.upper()
    unknown = set(protein) - set(AMINO_ACID_CODONS) - {"*"}
    if unknown:
        raise ValueError(f"Unrecognized amino acids: {sorted(unknown)}")

    variants = []
    for k in range(n_sequences):
        variant_seed = None if seed is None else seed + k
        rng = random.Random(variant_seed)
        variants.append(_generate_one_variant(protein, host, mode, include_stop, rng))

    return {
        "host": host,
        "mode": mode,
        "protein": protein,
        "aa_length": len(protein),
        "n_sequences": n_sequences,
        "variants": variants,
    }