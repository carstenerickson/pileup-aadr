# Contributing to pileup-aadr

Thanks for your interest. This document covers the local dev loop, the
integration-test setup, lint/type expectations, and the release process.
For user-facing install + usage docs see [README.md](README.md). For a
tour of the codebase — architecture, key contracts, common-task recipes —
see [DEVELOPMENT.md](DEVELOPMENT.md); read it before your first PR.

## Dev install

```bash
git clone https://github.com/carstenerickson/pileup-aadr.git
cd pileup-aadr
pip install -e ".[dev]"
```

The `[dev]` extra pulls in `pytest`, `pytest-cov`, `pytest-mock`, `mypy`,
and `ruff`. No additional system packages are required for the unit suite.

## Running the test suite

```bash
pytest                          # unit suite (~230 tests, ~5s)
pytest tests/integration/ -v    # integration tests (skip silently w/o binaries)
pytest -m slow                  # the gated >10s tests (opt-in)
```

The unit suite stubs every external binary via `pytest-mock`. CI runs it on
a 6-cell matrix (Python 3.11/3.12/3.13 × Linux/macOS); locally any single
combo is fine.

### Integration tests against real binaries

`tests/integration/` runs against real `samtools` / `pileupCaller` /
`Picard` / `mosdepth` / `pgen-samplebind` / AT2 and skips silently when
they're missing. To opt in locally, point the env at your installs:

```bash
# Picard (homebrew on macOS): the wrapper at /opt/homebrew/bin/picard exists
# but the version probe needs the JAR + java on PATH directly.
export PICARD_JAR=/opt/homebrew/Cellar/picard-tools/3.4.0/libexec/picard.jar
export PATH=/opt/homebrew/opt/openjdk/bin:$PATH

# pgen-samplebind (sister tool from the same author; install separately and
# add its venv's bin to PATH). The integration tests need it on PATH.
export PATH=/path/to/pgen-samplebind/.venv/bin:$PATH

# AT2 (R package): install the production fork — upstream 2.0.10 doesn't
# read pgen-samplebind's PFILE output cleanly:
#   Rscript -e 'remotes::install_github("carstenerickson/admixtools", ref="production/v1.0")'

# LLD #19 full-chain f2 invariant test (Layer B, gated on a real BAM):
export PILEUP_AADR_LLD19_BAM=/path/to/sample.bam
pytest tests/integration/test_lld19_full_chain.py -v
```

Without these, `pytest tests/integration/` exits cleanly with skips. CI's
`bio-tools` job does this setup via bioconda automatically — local-only
setup is for fast iteration without round-tripping through CI.

## Lint and type-checking

```bash
ruff check .
ruff format --check .
mypy pileup_aadr
```

All three must pass on `main`. CI enforces them in the `lint` job. Ruff
config lives in `pyproject.toml`; mypy is configured in `[tool.mypy]`
with `strict = true`. Per-file ignores are documented inline in
`pyproject.toml` with rationale.

## Benchmarks

`benchmarks/` uses `bench_*.py` (not `test_*.py`) so the default `pytest`
run doesn't collect them. To run:

```bash
pytest benchmarks/ -v
```

## Commit + PR conventions

- Short subject (≤72 chars), blank line, detailed body when warranted
- Reference HLD/LLD section numbers when the change implements a spec item
- One logical change per commit; rebase before opening a PR
- CI must be green before merge

## Release process

Releases are cut by tagging `vX.Y.Z` on `main`:

1. Bump `version` in `pyproject.toml`
2. Add a `## [X.Y.Z] — YYYY-MM-DD` heading + entries to `CHANGELOG.md`
3. Commit on `main`, then `git tag -a vX.Y.Z -m "..."` and push the tag
4. `.github/workflows/release.yml` picks up the tag, runs the wheel-smoke
   matrix (6 cells), and publishes to PyPI via OIDC trusted publishing
5. Verify with `pip install --upgrade pileup-aadr`

No PyPI token is configured; the workflow uses GitHub's OIDC integration.
The PyPI project is locked to publishing only from this repo's
`release.yml` workflow.

## Where things live

| Path                       | Contents                                              |
|----------------------------|-------------------------------------------------------|
| `pileup_aadr/`             | Package source                                        |
| `pileup_aadr/cli.py`       | Click root + subcommand dispatch                      |
| `pileup_aadr/extract_*.py` | The 4-stage extract pipeline                          |
| `pileup_aadr/data/`        | Bundled UCSC chain (SHA-verified at startup)          |
| `tests/`                   | Unit tests (mocked binaries)                          |
| `tests/integration/`       | Real-binary tests; skip cleanly without setup         |
| `benchmarks/`              | Perf harness, opt-in via explicit invocation          |
| `scripts/`                 | One-off ops scripts (e.g., LLD #19 baseline regen)    |

The HLD and LLD design docs live outside this repo. CHANGELOG entries
reference HLD/LLD section IDs where applicable.
