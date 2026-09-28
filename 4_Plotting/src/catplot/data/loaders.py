"""Explicit tabular file loading with conservative format inference."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from catplot.errors import DataError

_FORMAT_ALIASES = {
    "csv": "csv",
    "tsv": "tsv",
    "txt": "whitespace",
    "text": "whitespace",
    "dat": "whitespace",
    "whitespace": "whitespace",
    "excel": "excel",
    "xls": "excel",
    "xlsx": "excel",
    "parquet": "parquet",
    "feather": "feather",
    "json": "json",
}

_SUFFIX_FORMATS = {
    ".csv": "csv",
    ".tsv": "tsv",
    ".txt": "whitespace",
    ".dat": "whitespace",
    ".xls": "excel",
    ".xlsx": "excel",
    ".parquet": "parquet",
    ".feather": "feather",
    ".json": "json",
}


def infer_format(path: str | Path, file_format: str | None = None) -> str:
    """Resolve an explicit format or infer one from a known file extension."""
    if file_format is not None:
        normalized = file_format.lower().lstrip(".")
        try:
            return _FORMAT_ALIASES[normalized]
        except KeyError as exc:
            choices = ", ".join(sorted(set(_FORMAT_ALIASES.values())))
            raise DataError(f"Unknown data format {file_format!r}. Choose from: {choices}") from exc

    suffix = Path(path).suffix.lower()
    try:
        return _SUFFIX_FORMATS[suffix]
    except KeyError as exc:
        raise DataError(
            f"Cannot infer a data format from extension {suffix or '<none>'!r}; "
            "set file_format explicitly"
        ) from exc


def _normalise_headers(frame: pd.DataFrame, *, strip_headers: bool) -> pd.DataFrame:
    headers = [str(column).strip() if strip_headers else str(column) for column in frame.columns]
    duplicates = sorted({name for name in headers if headers.count(name) > 1})
    if duplicates:
        raise DataError(f"Column names collide after header normalization: {duplicates}")
    frame.columns = headers
    return frame


def load_table(
    path: str | Path,
    *,
    file_format: str | None = None,
    delimiter: str | None = None,
    sheet_name: str | int = 0,
    header: int | None = 0,
    names: Sequence[str] | None = None,
    comment: str | None = "#",
    decimal: str = ".",
    encoding: str = "utf-8",
    na_values: Sequence[Any] | None = None,
    columns: Sequence[str] | None = None,
    strip_headers: bool = True,
) -> pd.DataFrame:
    """Load one table and attach source metadata to ``DataFrame.attrs``.

    Format-specific engines are invoked only for their matching formats. No
    delimiter or header guessing occurs beyond the documented extension map.
    """
    source = Path(path).expanduser()
    if not source.is_file():
        raise DataError(f"Data file does not exist: {source}")
    resolved = infer_format(source, file_format)

    common = {
        "header": header,
        "names": list(names) if names is not None else None,
        "na_values": list(na_values) if na_values is not None else None,
    }
    try:
        if resolved == "csv":
            frame = pd.read_csv(
                source,
                sep=delimiter or ",",
                comment=comment,
                decimal=decimal,
                encoding=encoding,
                **common,
            )
        elif resolved == "tsv":
            frame = pd.read_csv(
                source,
                sep=delimiter or "\t",
                comment=comment,
                decimal=decimal,
                encoding=encoding,
                **common,
            )
        elif resolved == "whitespace":
            frame = pd.read_csv(
                source,
                sep=delimiter or r"\s+",
                comment=comment,
                decimal=decimal,
                encoding=encoding,
                **common,
            )
        elif resolved == "excel":
            frame = pd.read_excel(
                source,
                sheet_name=sheet_name,
                header=header,
                names=list(names) if names is not None else None,
                na_values=list(na_values) if na_values is not None else None,
            )
        elif resolved == "parquet":
            frame = pd.read_parquet(source)
        elif resolved == "feather":
            frame = pd.read_feather(source)
        elif resolved == "json":
            frame = pd.read_json(source, encoding=encoding)
        else:  # pragma: no cover - infer_format guarantees this branch is unreachable.
            raise DataError(f"Unsupported data format: {resolved}")
    except ImportError as exc:
        raise DataError(
            f"Loading {resolved!r} data needs an optional IO dependency: {exc}"
        ) from exc
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise DataError(f"Could not load {source}: {exc}") from exc

    if not isinstance(frame, pd.DataFrame):
        raise DataError(f"Loading {source} did not produce one table")
    if frame.empty and len(frame.columns) == 0:
        raise DataError(f"Data file contains no table: {source}")

    frame = _normalise_headers(frame, strip_headers=strip_headers)
    if columns is not None:
        from .columns import require_columns

        requested = [str(column) for column in columns]
        require_columns(frame, requested, source=str(source))
        frame = frame.loc[:, requested].copy()

    frame.attrs["source_path"] = str(source.resolve())
    frame.attrs["source_format"] = resolved
    return frame


def load_many(
    sources: Mapping[str, str | Path],
    *,
    options: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, pd.DataFrame]:
    """Load several explicitly named sources in input order."""
    if not sources:
        raise DataError("At least one named data source is required")
    unknown_options = set(options or {}) - set(sources)
    if unknown_options:
        raise DataError(f"Load options reference unknown sources: {sorted(unknown_options)}")

    loaded: dict[str, pd.DataFrame] = {}
    for source_id, source_path in sources.items():
        if not source_id or not isinstance(source_id, str):
            raise DataError("Every data source needs a non-empty string identifier")
        loaded[source_id] = load_table(source_path, **dict((options or {}).get(source_id, {})))
    return loaded

