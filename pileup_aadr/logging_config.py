"""Logging configuration for pileup-aadr.

Per HLD §"Logging architecture":
- Default: human-readable format to stderr
- $PILEUP_AADR_JSON_LOGS=1 → JSON-lines format to stderr (for ancestry-pipeline-tool's
  stderr-streams-to-disk discipline)
- Per-module loggers via `log = logging.getLogger(__name__)` pattern
- ruff T20 rule (in pyproject.toml) forbids `print()` in module code
"""
import json
import logging
import os
import sys
from datetime import UTC, datetime
from typing import Any, Final

_HUMAN_FORMAT: Final = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_HUMAN_DATEFMT: Final = "%Y-%m-%d %H:%M:%S"


class JsonLinesFormatter(logging.Formatter):
    """One JSON object per record, written as a single line.

    For framework integration (ancestry-pipeline-tool's stderr-streams-to-disk discipline
    can parse without regex). Extra fields attached via `logger.info(..., extra={...})`
    are passed through; reserved LogRecord attributes are filtered automatically (M16 fix:
    derived from a fresh LogRecord rather than a hardcoded list, so Python version
    differences like 3.12's `taskName` are handled transparently).
    """

    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, Any] = {
            "ts": datetime.fromtimestamp(record.created, tz=UTC).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "msg": record.getMessage(),
        }
        if record.exc_info:
            payload["exc_info"] = self.formatException(record.exc_info)
        # Include any extra fields attached to the record via logger.*(..., extra={...})
        for key, value in record.__dict__.items():
            if key in _LOGRECORD_RESERVED_KEYS or key.startswith("_"):
                continue
            try:
                json.dumps(value)  # check serializability
                payload[key] = value
            except (TypeError, ValueError):
                payload[key] = str(value)
        return json.dumps(payload, separators=(",", ":"))


def _compute_reserved_keys() -> frozenset[str]:
    """Derive the reserved-keys set from a fresh LogRecord rather than a hardcoded list.

    Catches Python version differences automatically (e.g., Python 3.12 added `taskName`).
    Adds two keys not in __dict__ at construction but populated at format time:
    `message` (set by Formatter.format) and `asctime` (set if format uses it).
    """
    sample = logging.makeLogRecord({})
    return frozenset(sample.__dict__.keys()) | {"message", "asctime"}


_LOGRECORD_RESERVED_KEYS: Final = _compute_reserved_keys()


def configure_logging(*, level: int = logging.INFO) -> None:
    """Initialize root logger. Called once at startup by `cli.py`.

    Args:
        level: log level for the root logger (DEBUG / INFO / WARNING / etc.)

    Behavior:
        - Default: human-readable format to stderr
        - $PILEUP_AADR_JSON_LOGS=1 → JSON-lines format to stderr (for framework ingestion)
        - Re-callable safely (clears existing handlers first; useful for tests)
    """
    root = logging.getLogger()
    # Clear any pre-existing handlers (idempotent across test runs)
    for handler in list(root.handlers):
        root.removeHandler(handler)

    handler = logging.StreamHandler(sys.stderr)
    if os.environ.get("PILEUP_AADR_JSON_LOGS") == "1":
        handler.setFormatter(JsonLinesFormatter())
    else:
        handler.setFormatter(logging.Formatter(fmt=_HUMAN_FORMAT, datefmt=_HUMAN_DATEFMT))

    root.addHandler(handler)
    root.setLevel(level)

    # Quiet noisy third-party loggers
    logging.getLogger("urllib3").setLevel(logging.WARNING)


__all__ = ["JsonLinesFormatter", "configure_logging"]
