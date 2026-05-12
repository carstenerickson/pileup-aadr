"""Error class hierarchy + format_error compliance tests.

Maps to LLD test #40 (error message format compliance): every named error class
must follow the What/Why/Fix template.
"""
from __future__ import annotations

import inspect

import pytest

from pileup_aadr import errors
from pileup_aadr.errors import (
    AADRDuplicateRsidError,
    AADRParseError,
    BAMBaiMissingError,
    BAMBuildMismatchError,
    BAMParseError,
    BAMSampleNameAmbiguous,
    ChainFileNotFound,
    ChainFileSHAError,
    CoverageGateFailure,
    JavaNotFoundError,
    LiftoverYieldError,
    MosdepthNotFoundError,
    OutputExistsError,
    OutputLockHeldError,
    PicardNotFoundError,
    PileupAadrError,
    PileupAadrInternalError,
    PileupCallerNotFoundError,
    ReferenceFastaBuildMismatchError,
    ReferenceFastaNotFound,
    SamtoolsNotFoundError,
    ToolNotFoundError,
    ToolSubprocessError,
    ToolVersionError,
    UnsupportedAADRBuild,
    UnsupportedReferenceBuild,
    format_error,
)


def _all_concrete_error_classes() -> list[type[PileupAadrError]]:
    """Iterate every concrete (raisable) PileupAadrError subclass declared in errors.py."""
    return [
        cls
        for _, cls in inspect.getmembers(errors, inspect.isclass)
        if issubclass(cls, PileupAadrError)
        and cls is not PileupAadrError
        and cls.exit_code in (1, 2, 3, 4)
    ]


def test_count_of_named_classes() -> None:
    """25 concrete classes total (HLD's "20 named" was a logical-category count;
    counting ToolNotFoundError + 5 subclasses individually gives 25 actually-raisable
    distinct classes).

    Distribution per `test_exit_code_distribution` below: 2 + 9 + 6 + 8 = 25.
    """
    assert len(_all_concrete_error_classes()) == 25


@pytest.mark.parametrize(
    "cls",
    [
        LiftoverYieldError,
        CoverageGateFailure,
        ChainFileNotFound,
        ChainFileSHAError,
        ReferenceFastaNotFound,
        ReferenceFastaBuildMismatchError,
        BAMSampleNameAmbiguous,
        OutputLockHeldError,
        ToolSubprocessError,
        BAMBaiMissingError,
        BAMParseError,
        UnsupportedReferenceBuild,
        UnsupportedAADRBuild,
        BAMBuildMismatchError,
        AADRDuplicateRsidError,
        AADRParseError,
        PileupAadrInternalError,
        OutputExistsError,
        ToolNotFoundError,
        SamtoolsNotFoundError,
        PileupCallerNotFoundError,
        PicardNotFoundError,
        JavaNotFoundError,
        MosdepthNotFoundError,
        ToolVersionError,
    ],
)
def test_what_why_fix_compliance(cls: type[PileupAadrError]) -> None:
    """Every error formats to a message containing What:, Why:, Fix:."""
    e = cls(what="<test what>", why="<test why>", fix="<test fix>")
    msg = format_error(e)
    assert "What:" in msg
    assert "Why:" in msg
    assert "Fix:" in msg
    assert "<test what>" in msg
    assert "<test why>" in msg
    assert "<test fix>" in msg
    assert cls.__name__ in msg


def test_exit_codes_all_in_valid_range() -> None:
    """Every concrete class declares exit_code in {1, 2, 3, 4}."""
    for cls in _all_concrete_error_classes():
        assert cls.exit_code in (1, 2, 3, 4), f"{cls.__name__}.exit_code = {cls.exit_code}"


def test_exit_code_distribution() -> None:
    """Counts per exit code: 2 / 9 / 6 / 8 = 25 total."""
    counts = {1: 0, 2: 0, 3: 0, 4: 0}
    for cls in _all_concrete_error_classes():
        counts[cls.exit_code] += 1
    # 2 exit-1: LiftoverYieldError + CoverageGateFailure
    assert counts[1] == 2, f"exit-1 count: {counts[1]}"
    # 9 exit-2: ChainFileNotFound + ChainFileSHAError + ReferenceFastaNotFound +
    # ReferenceFastaBuildMismatchError + BAMSampleNameAmbiguous + OutputLockHeldError +
    # ToolSubprocessError + BAMBaiMissingError + BAMParseError
    assert counts[2] == 9, f"exit-2 count: {counts[2]}"
    # 6 exit-3: UnsupportedReferenceBuild + UnsupportedAADRBuild + BAMBuildMismatchError +
    # AADRDuplicateRsidError + AADRParseError + PileupAadrInternalError
    assert counts[3] == 6, f"exit-3 count: {counts[3]}"
    # 8 exit-4: OutputExistsError + ToolNotFoundError (base) + 5 ToolNotFoundError
    # subclasses + ToolVersionError
    assert counts[4] == 8, f"exit-4 count: {counts[4]}"


def test_subclass_without_exit_code_raises() -> None:
    """The __init_subclass__ guard rejects subclasses with bad exit_code."""
    with pytest.raises(TypeError, match=r"exit_code"):

        class BadError(PileupAadrError):
            exit_code = 99  # type: ignore[misc]


def test_init_requires_keyword_only_args() -> None:
    """All three of what/why/fix must be passed as kwargs."""
    with pytest.raises(TypeError):
        LiftoverYieldError("positional", "args", "fail")  # type: ignore[misc]


def test_format_error_starts_with_pileup_aadr_prefix() -> None:
    e = LiftoverYieldError(what="x", why="y", fix="z")
    msg = format_error(e)
    assert msg.startswith("pileup-aadr:")
