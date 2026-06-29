# Changelog

All notable changes to this project will be documented in this file.

# Changelog

All notable changes to `mol_bio_toolkit` will be documented in this file.

This project adheres to [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.1.0] - 2026-06-29

### Added

#### Package infrastructure
- `src`-layout package structure (`src/mol_bio_toolkit/`) with editable install support (`pip install -e .`)
- `pyproject.toml` with full project metadata, dependency declarations (`biopython`, `pandas`, `matplotlib`), and `requires-python = ">=3.10"`
- MIT License (2026, JD Morales)
- `.gitignore` for Python projects
- GitHub Actions CI pipeline (`.github/workflows/tests.yml`) running `pytest` on Python 3.10, 3.11, and 3.12 on every push and pull request

#### Core modules
- **`orf_finder`** — ORF detection across 3 forward reading frames with configurable minimum length (`min_length`, default 100 bp). Returns ORF start position, end position, reading frame, and translated amino acid sequence. Validated on sfGFP (BBa_I746916, 717 bp): correctly identifies 11 ORFs ≥100 bp; primary ORF recovered at frame +1, positions 0–717 (238 aa).
- **`composition`** — Per-nucleotide frequency (A, T, G, C) and GC% calculation. Validated on sfGFP: GC% = 38.77%.
- **`primer_design`** — Iterative Tm-based primer length search using the Wallace rule (Tm = 2·[A+T] + 4·[G+C] with GC-length correction for primers >13 nt), with GC-clamp check and heuristic hairpin detection via self-complementarity scoring. Returns forward and reverse primers with Tm, GC%, and length. Validated on sfGFP 717 bp amplicon: forward 30 nt / Tm 58.9°C / GC 40.0%; reverse 30 nt / Tm 56.2°C / GC 33.3%.
- **`codon_optimizer`** — Codon harmonization using published frequency tables (Nakamura et al. 2000). Supports *E. coli* K-12 (`e_coli_k12`) and *S. cerevisiae* (`s_cerevisiae`). Preserves amino acid sequence. Returns optimized nucleotide sequence, number of codon changes, and GC% shift. Validated on sfGFP 239 codons: *E. coli* optimization changes 128/239 codons, GC 38.8% → 42.1%; *S. cerevisiae* changes 111/239 codons, GC 38.8% → 35.8%.

#### Tests
- 21 pytest tests covering all four modules
- sfGFP (BBa_I746916, iGEM Registry, CC-BY-4.0; Pédelacq et al. 2006) used as biological reference sequence throughout

#### Documentation
- `README.md` with motivation, biological background, module overview, usage examples with verified numerical outputs, limitations, and citation

### Known limitations (tracked as open issues)
- ORF detection scans forward frames only; reverse complement scanning not yet implemented
- Tm calculation uses Wallace rule / GC-length correction; nearest-neighbor thermodynamics (SantaLucia 1998) and salt correction not implemented
- Hairpin detection is heuristic (sequence-level self-complementarity); ΔG estimation not implemented
- Codon optimization available for *E. coli* K-12 and *S. cerevisiae* only

---

## References

- Pédelacq, J.D., Cabantous, S., Tran, T., Terwilliger, T.C., & Waldo, G.S. (2006). Engineering and characterization of a superfolder green fluorescent protein. *Nature Biotechnology*, 24(1), 79–88. https://doi.org/10.1038/nbt1172
- Nakamura, Y., Gojobori, T., & Ikemura, T. (2000). Codon usage tabulated from international DNA sequence databases: status for the year 2000. *Nucleic Acids Research*, 28(1), 292. https://doi.org/10.1093/nar/28.1.292

[Unreleased]: https://github.com/JDMorales11/mol_bio_toolkit/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/JDMorales11/mol_bio_toolkit/releases/tag/v0.1.0
