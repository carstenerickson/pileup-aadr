"""Exception hierarchy for pileup-aadr.

Per HLD §"Error class taxonomy": every exception inherits from `PileupAadrError`,
carries an `exit_code` class attribute in {1, 2, 3, 4}, and is constructed with three
keyword-only arguments (`what`, `why`, `fix`) for the standard message-format template.
The orchestrator's outermost handler catches `PileupAadrError`, formats the message
per the template, writes to stderr, and returns the class's `exit_code`.

Exit-code mapping (consumed by ancestry-pipeline-tool's `AncestryPipelineError`):
    0  success
    1  soft-validation failure (LiftoverYieldError, CoverageGateFailure)
    2  I/O failure (chain/FASTA/BAM not found, subprocess crashed, lock held)
    3  invariant violation (build mismatch, AADR malformed, defensive sanity check)
    4  usage error (bad CLI args, missing/wrong-version external binary)
"""


class PileupAadrError(Exception):
    """Base exception for all pileup-aadr errors. Carries exit_code + What/Why/Fix message format.

    All subclasses set `exit_code` as a class attribute (one of {1, 2, 3, 4}).
    Constructor takes three keyword-only str args (`what`, `why`, `fix`) to enforce the
    HLD's `What:/Why:/Fix:` template — no positional message strings allowed.

    Subclass discipline: every subclass MUST set exit_code (validated via __init_subclass__).
    """

    exit_code: int = -1  # subclass must override

    def __init_subclass__(cls, **kwargs: object) -> None:
        super().__init_subclass__(**kwargs)
        # Skip validation for explicit intermediate base classes that don't set their own
        # exit_code (e.g., ToolNotFoundError); they inherit from a class that already passed.
        if cls.exit_code not in (1, 2, 3, 4):
            raise TypeError(
                f"{cls.__name__} must define exit_code class attr in {{1, 2, 3, 4}} "
                f"(got {cls.exit_code})"
            )

    def __init__(self, *, what: str, why: str, fix: str) -> None:
        self.what = what
        self.why = why
        self.fix = fix
        super().__init__(f"{type(self).__name__}: {why}")

    def format(self) -> str:
        """Render the standard What/Why/Fix message block.

        Format:
            pileup-aadr: <ClassName>: <one-line description (uses why)>

              What:  <self.what>
              Why:   <self.why>
              Fix:   <self.fix>
        """
        return (
            f"pileup-aadr: {type(self).__name__}: {self.why}\n"
            f"\n"
            f"  What:  {self.what}\n"
            f"  Why:   {self.why}\n"
            f"  Fix:   {self.fix}\n"
        )


def format_error(e: PileupAadrError) -> str:
    """Render a PileupAadrError as the standard What/Why/Fix message block."""
    return e.format()


# --- Exit code 1: soft-validation failures ---


class LiftoverYieldError(PileupAadrError):
    """Picard liftover yield below --liftover-yield-fail-pct (default 70%)."""

    exit_code = 1


class CoverageGateFailure(PileupAadrError):
    """Non-missing autosomal calls below --min-coverage (default 500K)."""

    exit_code = 1


# --- Exit code 2: I/O failures ---


class ChainFileNotFound(PileupAadrError):
    """Explicit --chain PATH does not exist.

    The bundled chain is always present after install; this only fires for --chain override.
    """

    exit_code = 2


class ChainFileSHAError(PileupAadrError):
    """Bundled chain SHA mismatches expected, OR --strict-chain-sha set on user-supplied
    chain that mismatches the package-pinned SHA."""

    exit_code = 2


class ReferenceFastaNotFound(PileupAadrError):
    """Target FASTA missing from --ref-fasta + env + BAM @PG."""

    exit_code = 2


class ReferenceFastaBuildMismatchError(PileupAadrError):
    """Target FASTA's chr1 length doesn't match BAM @SQ chr1 length (LLD v2.1 H11)."""

    exit_code = 2


class BAMSampleNameAmbiguous(PileupAadrError):
    """BAM has multiple @RG SM: values that disagree AND --sample-name not given."""

    exit_code = 2


class OutputLockHeldError(PileupAadrError):
    """Advisory lock on <prefix>.lock held by another process."""

    exit_code = 2


class ToolSubprocessError(PileupAadrError):
    """A subprocess (Picard / samtools / pileupCaller / mosdepth) exited non-zero
    or was killed by signal."""

    exit_code = 2


class BAMBaiMissingError(PileupAadrError):
    """BAM file present but .bai index not found.

    v0.1: raised only by `validate` subcommand; the extract orchestrator lets samtools
    mpileup fail later if the index is missing. v0.2 may promote this to extract-startup.
    """

    exit_code = 2


class BAMParseError(PileupAadrError):
    """BAM/CRAM file unparseable (wrong format, truncated, gzip-but-not-BAM)."""

    exit_code = 2


# --- Exit code 3: invariant violations ---


class UnsupportedReferenceBuild(PileupAadrError):
    """BAM @SQ chr1 length doesn't match hg19 or hg38 ± 1Mb (e.g., T2T-CHM13).

    User can override with --bam-build to skip detection.
    """

    exit_code = 3


class UnsupportedAADRBuild(PileupAadrError):
    """AADR .snp chr1 max position doesn't match hg19 or hg38."""

    exit_code = 3


class BAMBuildMismatchError(PileupAadrError):
    """Reserved for v0.2: a strict-build-check mode that compares --bam-build override
    against detected length. v0.1 honors the override silently (no `--no-lift` flag in v0.1).
    """

    exit_code = 3


class AADRDuplicateRsidError(PileupAadrError):
    """AADR .snp contains duplicate SNP-name (Stage 4 rsID join requires unique IDs)."""

    exit_code = 3


class AADRParseError(PileupAadrError):
    """AADR .snp row has wrong column count, non-ACGT alleles, or unparseable position."""

    exit_code = 3


class PileupAadrInternalError(PileupAadrError):
    """Defensive sanity-check failure (e.g., Stage 4 allele mismatch beyond tolerance,
    Picard stderr unparseable suggesting tool format change). Should never fire under
    correct usage; surface a bug-report path."""

    exit_code = 3


# --- Exit code 4: usage errors / missing dependencies ---


class OutputExistsError(PileupAadrError):
    """Output file(s) exist at startup and --overwrite not set."""

    exit_code = 4


class ToolNotFoundError(PileupAadrError):
    """Generic tool-not-found base. Subclasses below for specific tools."""

    exit_code = 4


class SamtoolsNotFoundError(ToolNotFoundError):
    """samtools binary not found on PATH."""


class PileupCallerNotFoundError(ToolNotFoundError):
    """pileupCaller binary not found on PATH."""


class PicardNotFoundError(ToolNotFoundError):
    """picard.jar not found via $PICARD_JAR, conda paths, wrapper script, or known install paths."""


class JavaNotFoundError(ToolNotFoundError):
    """java binary not found on PATH (transitive Picard dependency)."""


class MosdepthNotFoundError(ToolNotFoundError):
    """mosdepth binary not found on PATH (used by `coverage` subcommand)."""


class ToolVersionError(PileupAadrError):
    """A required binary is present but version too old."""

    exit_code = 4


__all__ = [
    "PileupAadrError",
    "format_error",
    # Exit 1
    "LiftoverYieldError",
    "CoverageGateFailure",
    # Exit 2
    "ChainFileNotFound",
    "ChainFileSHAError",
    "ReferenceFastaNotFound",
    "ReferenceFastaBuildMismatchError",
    "BAMSampleNameAmbiguous",
    "OutputLockHeldError",
    "ToolSubprocessError",
    "BAMBaiMissingError",
    "BAMParseError",
    # Exit 3
    "UnsupportedReferenceBuild",
    "UnsupportedAADRBuild",
    "BAMBuildMismatchError",
    "AADRDuplicateRsidError",
    "AADRParseError",
    "PileupAadrInternalError",
    # Exit 4
    "OutputExistsError",
    "ToolNotFoundError",
    "SamtoolsNotFoundError",
    "PileupCallerNotFoundError",
    "PicardNotFoundError",
    "JavaNotFoundError",
    "MosdepthNotFoundError",
    "ToolVersionError",
]
