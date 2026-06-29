# Contributing to mol_bio_toolkit

Thank you for your interest in contributing to `mol_bio_toolkit`. This document explains how to report issues, propose changes, and submit code.

---

## Table of contents

- [Code of conduct](#code-of-conduct)
- [Reporting bugs](#reporting-bugs)
- [Suggesting enhancements](#suggesting-enhancements)
- [Development setup](#development-setup)
- [Submitting a pull request](#submitting-a-pull-request)
- [Coding conventions](#coding-conventions)
- [Running tests](#running-tests)
- [Commit message style](#commit-message-style)

---

## Code of conduct

This project follows the [Contributor Covenant](https://www.contributor-covenant.org/) Code of Conduct (v2.1). By participating, you agree to uphold a welcoming and respectful environment for everyone.

---

## Reporting bugs

Before opening a new issue, please search existing issues to avoid duplicates.

When filing a bug report, include:

- A short, descriptive title
- The Python version and operating system you are using
- The exact command or code that triggers the bug
- The full error message or unexpected output
- The expected behavior

Use the **Bug report** issue template on GitHub.

---

## Suggesting enhancements

Feature requests are welcome. Open an issue using the **Feature request** template and describe:

- The biological or computational problem the feature would address
- A proposed interface or API sketch (if you have one)
- Any relevant references (papers, existing tools)

Known planned enhancements are already tracked as open issues. Check those before opening a duplicate.

---

## Development setup

```bash
# 1. Fork and clone the repository
git clone https://github.com/<your-username>/mol_bio_toolkit.git
cd mol_bio_toolkit

# 2. Create and activate a conda environment
conda create -n mol_bio_toolkit python=3.11
conda activate mol_bio_toolkit

# 3. Install the package in editable mode with dev dependencies
pip install -e ".[dev]"
```

Dependencies are declared in `pyproject.toml`. Runtime dependencies are `biopython`, `pandas`, and `matplotlib`; the `[dev]` extra adds `pytest`, `pytest-cov`, `black`, and `ruff`.

---

## Submitting a pull request

1. Open or find an existing issue that describes the problem or feature.
2. Comment on the issue to let others know you are working on it.
3. Create a feature branch from `main`:

   ```bash
   git checkout -b feat/reverse-complement-orfs
   ```

4. Make your changes. Keep commits small and focused (one logical change per commit).
5. Add or update tests in `tests/` so that all 21+ existing tests continue to pass and new behavior is covered.
6. Run the full test suite locally before pushing (see [Running tests](#running-tests)).
7. Push your branch and open a pull request against `main`.
8. Reference the related issue in the PR description using `Closes #<issue-number>`.

Pull requests are reviewed by the maintainer. Feedback will be provided within a reasonable time. Please keep the PR focused — large, unrelated changes are harder to review.

---

## Coding conventions

- **Style:** Follow [PEP 8](https://peps.python.org/pep-0008/). Line length limit is 88 characters (compatible with `black`).
- **Type hints:** Add type annotations to all public functions. Example:

  ```python
  def find_orfs(sequence: str, min_length: int = 100) -> list[dict]:
      ...
  ```

- **Docstrings:** Use Google-style docstrings for all public functions and classes. Include `Args`, `Returns`, and `Raises` sections where applicable.
- **Biological accuracy:** Any change to a core algorithm (Tm calculation, codon tables, ORF detection logic) must include a citation to the primary literature in the docstring.
- **No unnecessary dependencies:** Do not introduce new runtime dependencies without discussion in an issue first.

---

## Running tests

```bash
# Run the full test suite
pytest tests/

# Run with verbose output
pytest tests/ -v

# Run a single test file
pytest tests/test_orf_finder.py
```

The test suite uses the sfGFP sequence (BBa_I746916, iGEM Registry, CC-BY-4.0) as a reference. All numerical outputs in tests must be verified against actual code execution — do not hard-code expected values without running the function first.

CI runs automatically on every push and pull request via GitHub Actions (`.github/workflows/tests.yml`), targeting Python 3.10, 3.11, and 3.12.

---

## Commit message style

This project uses [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <short description>

[optional body]

[optional footer]
```

Common types:

| Type       | Use for                                              |
|------------|------------------------------------------------------|
| `feat`     | New feature or biological capability                 |
| `fix`      | Bug fix                                              |
| `test`     | Adding or updating tests                             |
| `docs`     | Documentation changes only                           |
| `refactor` | Code restructuring with no functional change         |
| `ci`       | Changes to GitHub Actions or CI configuration        |
| `chore`    | Maintenance tasks (dependency updates, etc.)         |

Examples:

```
feat(orf_finder): add reverse complement frame scanning
fix(primer_design): correct Tm formula for primers >13 nt
test(codon_optimizer): add invariant amino acid sequence test for yeast host
docs(readme): update primer Tm output values after Wallace rule fix
```

---

## Questions

If you are unsure whether a change is in scope, open an issue and ask before investing time in an implementation. This keeps contributions aligned with the project's goals and avoids wasted effort.
