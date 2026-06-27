# mol_bio_toolkit

[![Tests](https://github.com/JDMorales11/mol_bio_toolkit/actions/workflows/tests.yml/badge.svg)](https://github.com/JDMorales11/mol_bio_toolkit/actions/workflows/tests.yml)

> A Python toolkit for core molecular biology sequence analysis, designed for researchers working at the intersection of wet and dry lab workflows.

---

## 1. Motivation

Sequence-level analysis tasks which include identifying open reading frames, calculating primer melting temperatures, assessing nucleotide composition, and optimizing codon usage are foundational operations in molecular biology and synthetic biology workflows. Despite their ubiquity, these tasks are often performed using disconnected online tools or large framework dependencies that obscure the underlying logic and are difficult to integrate into reproducible pipelines. `mol_bio_toolkit` provides a Python implementation of these core operations built on BioPython, designed to be readable, testable, and composable within larger bioinformatics workflows. It is particularly suited for synthetic biology applications where sequence design and optimization are iterative processes tightly coupled to experimental decisions.

---

## 2. Biological Background

Gene expression in heterologous systems depends critically on sequence-level properties that span multiple biological layers. Open reading frames define the translatable units of an mRNA, but their efficiency is modulated by codon usage bias (the non-uniform frequency with which synonymous codons are used across organisms). When a gene from one organism is expressed in another (e.g., a human gene expressed in *E. coli*), rare codons in the host can stall ribosomes and reduce protein yield. Codon optimization replaces rare codons with host-preferred synonyms while preserving the amino acid sequence, a strategy widely applied in recombinant protein production and synthetic gene design.

Primer design is equally critical: the melting temperature (Tm) of a primer determines the annealing conditions for PCR, and small sequence changes can shift Tm by several degrees, affecting specificity and yield. Nucleotide composition, particularly GC content, underlies both Tm and secondary structure propensity, making it a first-pass quality metric for any synthetic or amplified sequence.

This toolkit implements these analyses using standard thermodynamic models and published codon usage tables (Nakamura et al. 2000), providing transparent, auditable implementations that can be validated against experimental outcomes.

---

## 3. Methods & Implementation

### Installation

Requires Python 3.8+. Dependencies: BioPython, pandas.

```bash
git clone https://github.com/JDMorales11/mol_bio_toolkit.git
cd mol_bio_toolkit
pip install -e .
```

### Module Overview

| Module | Description |
|---|---|
| `orf_finder` | ORF detection across 3 forward reading frames |
| `composition` | Nucleotide frequency and GC% calculation |
| `primer_design` | Iterative Tm-based primer length search with GC-clamp and hairpin detection |
| `codon_optimizer` | Codon harmonization for *E. coli* K-12 and *S. cerevisiae* |

### Usage

**ORF Finder**
```python
from mol_bio_toolkit.orf_finder import find_orfs

orfs = find_orfs("ATGAAACCCGGGTTTTAA", min_length=100)
for orf in orfs:
    print(orf)
```

**Sequence Composition**
```python
from mol_bio_toolkit.composition import nucleotide_composition

print(nucleotide_composition("ATGCGATAGCGCGC"))
# Returns: {"A": ..., "T": ..., "G": ..., "C": ..., "GC%": ...}
```

**Primer Design**
```python
from mol_bio_toolkit.primer_design import design_primers

result = design_primers(sequence, amplicon_start=0, amplicon_end=717)
print(result["forward"])
print(result["reverse"])
```

**Codon Optimizer**
```python
from mol_bio_toolkit.codon_optimizer import codon_optimize

result = codon_optimize(cds_sequence, host="e_coli_k12")
print(result["optimized_seq"])
print(result["codon_changes"])
```

### Testing

```bash
pytest tests/
```

The test suite (21 tests) uses the sfGFP sequence (superfolder GFP, Pédelacq et al. 2006, 717 bp) as a reference — a widely used synthetic biology reporter with well-characterized codon usage in *E. coli*.

---

## 4. Results on Public Data

Analysis of the sfGFP coding sequence (iGEM Registry: [BBa_I746916](http://parts.igem.org/Part:BBa_I746916), 717 bp, *E. coli* codon-optimized, license: CC-BY-4.0) demonstrates the toolkit's core functionality on a biologically meaningful reference:

- **ORF detection** identified 11 ORFs (≥100 bp, forward frames only) across the 717 bp sequence. The primary sfGFP ORF was correctly recovered at frame +1, positions 0–717, encoding 238 amino acids, consistent with the annotated 239-codon CDS (238 aa + stop codon).
- **Nucleotide composition** returned a GC content of 38.77% (A: 34.59%, T: 26.64%, G: 19.67%, C: 19.11%), characteristic of *E. coli*-optimized synthetic genes with A/T-biased codon usage.
- **Primer design** targeting the full sfGFP CDS (717 bp amplicon) produced a 30 nt forward primer (5′-ATGCGTAAAGGAGAAGAACTTTTCACTGGA-3′, Tm = 58.9°C, GC = 40.0%) and a 30 nt reverse primer (5′-TTATTTGTATAGTTCATCCATGCCATGTGT-3′, Tm = 56.2°C, GC = 33.3%), both within standard PCR annealing range (55–65°C).
- **Codon optimization** of the sfGFP protein sequence (239 codons) using the *E. coli* K-12 table changed 128/239 codons and shifted GC content from 38.8% to 42.1%, reflecting substitution toward high-frequency *E. coli* codons. Optimization for *S. cerevisiae* changed 111/239 codons and shifted GC content to 35.8%, consistent with the characteristically lower GC bias of yeast-preferred codons.

> Pédelacq, J.D., Cabantous, S., Tran, T., Terwilliger, T.C., & Waldo, G.S. (2006). Engineering and characterization of a superfolder green fluorescent protein. *Nature Biotechnology*, 24(1), 79–88. https://doi.org/10.1038/nbt1172

> Nakamura, Y., Gojobori, T., & Ikemura, T. (2000). Codon usage tabulated from international DNA sequence databases: status for the year 2000. *Nucleic Acids Research*, 28(1), 292. https://doi.org/10.1093/nar/28.1.292

---

## 5. Limitations & Future Work

- **Forward frames only:** ORF detection currently scans only the 3 forward reading frames. Reverse complement ORF finding is not yet implemented and would be required for complete genome-scale annotation.
- **Simplified Tm model:** Tm calculation uses the Wallace rule for short oligos (<14 nt) and a GC-length correction for longer primers. Salt correction, DMSO adjustment, and full nearest-neighbor thermodynamics (SantaLucia 1998) are not yet implemented.
- **Heuristic hairpin detection:** The hairpin check identifies sequence-level self-complementarity but does not compute free energy (ΔG). Primer3 or equivalent tools should be used for final experimental validation.
- **Codon table coverage:** Optimization is currently available for *E. coli* K-12 and *S. cerevisiae*. Expanding to additional expression hosts is a planned extension.
- **No secondary structure prediction:** Codon optimization does not account for mRNA secondary structure, which can significantly affect translation efficiency in practice.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Citation

If you use this toolkit in your research, please cite:

> Morales-Cortes, J.D. (2026). mol_bio_toolkit: A Python toolkit for molecular biology sequence analysis. GitHub. https://github.com/JDMorales11/mol_bio_toolkit
