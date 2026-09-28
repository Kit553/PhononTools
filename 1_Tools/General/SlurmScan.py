#!/usr/bin/env python3

"""
slurm_scan

Scan Slurm output files and summarize current job state, diagnostics and
resource usage.

Default mode
------------
Without a subcommand, scan immediate subdirectories of the current directory.

For each directory, only the highest numeric Slurm job ID is considered the
current run. Older runs can be displayed with --history.

Examples
--------
    slurm_scan
    slurm_scan --all
    slurm_scan -p
    slurm_scan --history
    slurm_scan --details
    slurm_scan --recursive

    slurm_scan --state FAILED,TIMEOUT
    slurm_scan --diag OOM
    slurm_scan --diag HFI/PSM2,EAGER_MSG
    slurm_scan --diag OOM --paths

    slurm_scan --failed
    slurm_scan --failed restart_dirs.txt

    slurm_scan array 123456
    slurm_scan array 123456 123457 -p

Diagnostics
-----------
    OOM
    TIME_LIMIT
    HFI/PSM2
    EAGER_MSG
    MPI_ABORT
    SEGFAULT
    KILLED
    CANCELLED
    NODE_FAIL
    VASP_ERROR
    ERROR

Exit status is non-zero only for invocation, discovery, or parsing errors.
Detected job failures themselves do not make slurm_scan fail.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable, NoReturn, Optional


RUNNING_DIR_RE = re.compile(
    r"^Running directory:\s*(?P<running_dir>\S+)\s+in DirList line\s+"
    r"(?P<dirlist_line>\d+),\s+Array Task ID:\s*(?P<array_task_id>\d+),"
    r"\s+time\s+(?P<start_ts>.+?)\s*$",
    re.MULTILINE,
)

STEP_START_RE = re.compile(
    r"^STEP TEST START\s+(?P<start_ts>.+?)\s*$",
    re.MULTILINE,
)

FINISHED_AT_RE = re.compile(
    r"^Finished at\s+(?P<finish_ts>.+?)\s*$",
    re.MULTILINE,
)

JOB_ID_RE = re.compile(
    r"^Job ID:\s*(?P<job_id>\S+)\s*$",
    re.MULTILINE,
)

ARRAY_JOB_ID_RE = re.compile(
    r"^Array Job ID:\s*(?P<array_job_id>\S+)\s*$",
    re.MULTILINE,
)

STATE_RE = re.compile(
    r"^State:\s*(?P<state>[A-Z_]+)"
    r"(?:\s*\((?P<detail>[^)]*)\))?\s*$",
    re.MULTILINE,
)

NODES_RE = re.compile(
    r"^Nodes:\s*(?P<nodes>\S+)\s*$",
    re.MULTILINE,
)

CORES_PER_NODE_RE = re.compile(
    r"^Cores per node:\s*(?P<cpn>\S+)\s*$",
    re.MULTILINE,
)

CPU_UTIL_RE = re.compile(
    r"^CPU Utilized:\s*(?P<cpu_util>.+?)\s*$",
    re.MULTILINE,
)

CPU_EFF_RE = re.compile(
    r"^CPU Efficiency:\s*(?P<cpu_eff>[0-9.]+%)\s+of\s+"
    r"(?P<cpu_core_walltime>.+?)\s*$",
    re.MULTILINE,
)

WALL_RE = re.compile(
    r"^Job Wall-clock time:\s*(?P<wall>.+?)\s*$",
    re.MULTILINE,
)

MEM_UTIL_RE = re.compile(
    r"^Memory Utilized:\s*(?P<mem_util>.+?)"
    r"(?:\s+\(estimated maximum\))?\s*$",
    re.MULTILINE,
)

MEM_EFF_RE = re.compile(
    r"^Memory Efficiency:\s*(?P<mem_eff>[0-9.]+%)\s+of\s+"
    r"(?P<mem_req>.+?)\s*$",
    re.MULTILINE,
)

ARRAY_FILENAME_RE = re.compile(
    r"^slurm-(?P<array_id>\d+)_(?P<task_id>\d+)\.out$"
)

JOB_FILENAME_RE = re.compile(
    r"^slurm-(?P<job_id>\d+)\.out$"
)

MEM_RE = re.compile(
    r"^\s*(?P<value>[0-9.]+)\s*(?P<unit>[KMGTP]i?B)\s*$",
    re.IGNORECASE,
)


DIAGNOSTIC_PATTERNS: list[
    tuple[str, re.Pattern[str]]
] = [
    (
        "OOM",
        re.compile(
            r"(out[\s-]*of[\s-]*memory|"
            r"oom[_ -]?kill|"
            r"oom-kill|"
            r"exceeded.*memory|"
            r"memory limit|"
            r"cannot allocate memory)",
            re.IGNORECASE,
        ),
    ),
    (
        "TIME_LIMIT",
        re.compile(
            r"(due to time limit|"
            r"time limit exceeded|"
            r"timelimit|"
            r"\bTIMEOUT\b)",
            re.IGNORECASE,
        ),
    ),
    (
        "HFI/PSM2",
        re.compile(
            r"(can't open hfi|"
            r"hfi_userinit|"
            r"hfp_gen1_context_open|"
            r"assign_context command failed|"
            r"PSM2.*(?:error|can't)|"
            r"OFI endpoint open failed)",
            re.IGNORECASE,
        ),
    ),
    (
        "EAGER_MSG",
        re.compile(
            r"received eager message",
            re.IGNORECASE,
        ),
    ),
    (
        "MPI_ABORT",
        re.compile(
            r"(PMI_Abort|"
            r"MPI_Abort|"
            r"PMPI_Init.*(?:fatal|error)|"
            r"Fatal error in PMPI|"
            r"BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES)",
            re.IGNORECASE,
        ),
    ),
    (
        "SEGFAULT",
        re.compile(
            r"(segmentation fault|"
            r"segfault|"
            r"signal\s+11\b|"
            r"SIGSEGV)",
            re.IGNORECASE,
        ),
    ),
    (
        "KILLED",
        re.compile(
            r"(\bSIGKILL\b|"
            r"\bkilled\b|"
            r"signal\s+9\b)",
            re.IGNORECASE,
        ),
    ),
    (
        "CANCELLED",
        re.compile(
            r"\bCANCELLED\b",
            re.IGNORECASE,
        ),
    ),
    (
        "NODE_FAIL",
        re.compile(
            r"(NODE_FAIL|node failure)",
            re.IGNORECASE,
        ),
    ),
    (
        "VASP_ERROR",
        re.compile(
            r"(VERY BAD NEWS|"
            r"internal error in subroutine|"
            r"ERROR FEXCP|"
            r"ZBRENT: fatal error|"
            r"BRMIX: very serious problems|"
            r"EDDDAV: Call to ZHEGV failed)",
            re.IGNORECASE,
        ),
    ),
]

GENERIC_ERROR_RE = re.compile(
    r"\b(?:fatal|error|failed)\b",
    re.IGNORECASE,
)

KNOWN_DIAGNOSTICS = {
    name
    for name, _ in DIAGNOSTIC_PATTERNS
} | {"ERROR"}


@dataclass
class JobRecord:
    source_file: Path
    mode: str
    root: Path

    array_id: str = ""
    task_id: str = ""

    running_dir: str = ""
    dirlist_line: str = ""
    subdir: str = ""

    job_id: str = ""
    array_job_id: str = ""

    state: str = "NO REPORT"
    state_detail: str = ""

    nodes: str = ""
    cores_per_node: str = ""

    cpu_utilized: str = ""
    cpu_efficiency: str = ""
    cpu_core_walltime: str = ""

    wall_time: str = ""

    memory_utilized: str = ""
    memory_efficiency: str = ""
    requested_memory_total: str = ""
    requested_memory_per_core: str = ""

    start_time_str: str = ""
    finish_time_str: str = ""

    parse_note: str = ""

    diagnostic_counts: Counter[str] = field(
        default_factory=Counter
    )

    diagnostic_examples: dict[
        str,
        list[str],
    ] = field(
        default_factory=dict
    )

    @property
    def is_completed(self) -> bool:
        return self.state == "COMPLETED"

    @property
    def started_at(self) -> Optional[datetime]:
        return parse_log_datetime(
            self.start_time_str
        )

    @property
    def finished_at(self) -> Optional[datetime]:
        return parse_log_datetime(
            self.finish_time_str
        )

    @property
    def runtime_seconds(self) -> Optional[float]:
        start = self.started_at
        finish = self.finished_at

        if (
            start is not None
            and finish is not None
        ):
            runtime = (
                finish - start
            ).total_seconds()

            if runtime >= 0:
                return runtime

        return parse_slurm_duration_to_seconds(
            self.wall_time
        )

    @property
    def numeric_job_id(self) -> int:
        match = JOB_FILENAME_RE.match(
            self.source_file.name
        )

        if match:
            return int(
                match.group("job_id")
            )

        if self.job_id.isdigit():
            return int(
                self.job_id
            )

        return -1

    @property
    def diagnostics(self) -> str:
        if not self.diagnostic_counts:
            return "-"

        parts = []

        for name, _ in DIAGNOSTIC_PATTERNS:
            count = self.diagnostic_counts.get(
                name,
                0,
            )

            if count:
                parts.append(
                    f"{name}x{count}"
                    if count > 1
                    else name
                )

        count = self.diagnostic_counts.get(
            "ERROR",
            0,
        )

        if count:
            parts.append(
                f"ERRORx{count}"
                if count > 1
                else "ERROR"
            )

        return ",".join(parts) or "-"


def die(
    message: str,
    code: int = 1,
) -> "NoReturn":
    print(
        f"Error: {message}",
        file=sys.stderr,
    )

    raise SystemExit(code)


def warn(
    message: str,
) -> None:
    print(
        f"Warning: {message}",
        file=sys.stderr,
    )


def extract_last(
    pattern: re.Pattern[str],
    text: str,
    group: str,
) -> str:
    matches = list(
        pattern.finditer(text)
    )

    return (
        matches[-1]
        .group(group)
        .strip()
        if matches
        else ""
    )


def parse_requested_memory(
    value: str,
) -> tuple[str, str]:
    value = value.strip()

    match = re.match(
        r"^(?P<total>.+?)"
        r"\s*\((?P<per_core>.+?)\)\s*$",
        value,
    )

    if match:
        return (
            match.group("total").strip(),
            match.group("per_core").strip(),
        )

    return value, ""


def parse_log_datetime(
    value: str,
) -> Optional[datetime]:
    value = value.strip()

    if not value:
        return None

    parts = value.split()

    try:
        if len(parts) == 6:
            value_no_tz = " ".join(
                parts[:4]
                + [parts[5]]
            )

            return datetime.strptime(
                value_no_tz,
                "%a %b %d %H:%M:%S %Y",
            )

        if len(parts) == 5:
            return datetime.strptime(
                value,
                "%a %b %d %H:%M:%S %Y",
            )

    except ValueError:
        return None

    return None


def parse_slurm_duration_to_seconds(
    value: str,
) -> Optional[float]:
    value = value.strip()

    if not value:
        return None

    days = 0
    time_part = value

    if "-" in value:
        day_part, time_part = value.split(
            "-",
            1,
        )

        try:
            days = int(
                day_part
            )

        except ValueError:
            return None

    parts = time_part.split(":")

    if len(parts) != 3:
        return None

    try:
        hours = int(
            parts[0]
        )

        minutes = int(
            parts[1]
        )

        seconds = float(
            parts[2]
        )

    except ValueError:
        return None

    return (
        days * 86400
        + hours * 3600
        + minutes * 60
        + seconds
    )


def parse_percent(
    value: str,
) -> Optional[float]:
    value = (
        value.strip()
        .rstrip("%")
    )

    try:
        return float(
            value
        )

    except ValueError:
        return None


def parse_memory_to_bytes(
    value: str,
) -> Optional[float]:
    value = value.strip()

    if not value:
        return None

    match = MEM_RE.match(
        value
    )

    if not match:
        return None

    number = float(
        match.group("value")
    )

    unit = (
        match.group("unit")
        .upper()
    )

    factors = {
        "KB": 1000,
        "MB": 1000**2,
        "GB": 1000**3,
        "TB": 1000**4,
        "PB": 1000**5,

        "KIB": 1024,
        "MIB": 1024**2,
        "GIB": 1024**3,
        "TIB": 1024**4,
        "PIB": 1024**5,
    }

    return (
        number
        * factors[unit]
    )


def format_seconds_as_hms(
    seconds: float,
) -> str:
    total = int(
        round(seconds)
    )

    days, remainder = divmod(
        total,
        86400,
    )

    hours, remainder = divmod(
        remainder,
        3600,
    )

    minutes, seconds = divmod(
        remainder,
        60,
    )

    if days:
        return (
            f"{days}-"
            f"{hours:02d}:"
            f"{minutes:02d}:"
            f"{seconds:02d}"
        )

    return (
        f"{hours:02d}:"
        f"{minutes:02d}:"
        f"{seconds:02d}"
    )


def format_bytes(
    value: float,
) -> str:
    units = [
        "B",
        "KB",
        "MB",
        "GB",
        "TB",
        "PB",
    ]

    size = float(
        value
    )

    for unit in units:
        if (
            size < 1000.0
            or unit == units[-1]
        ):
            return (
                f"{size:.2f} {unit}"
            )

        size /= 1000.0

    return (
        f"{size:.2f} PB"
    )


def natural_key(
    value: str,
) -> list[object]:
    return [
        int(part)
        if part.isdigit()
        else part.lower()

        for part
        in re.split(
            r"(\d+)",
            value,
        )
    ]


def parse_filter(
    value: str | None,
) -> set[str]:
    if not value:
        return set()

    return {
        item.strip().upper()
        for item
        in re.split(
            r"[,;\s]+",
            value,
        )
        if item.strip()
    }


def atomic_write(
    path: Path,
    text: str,
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    descriptor, temp_name = tempfile.mkstemp(
        prefix=f".{path.name}.tmp.",
        dir=path.parent,
        text=True,
    )

    temp_path = Path(
        temp_name
    )

    try:
        with os.fdopen(
            descriptor,
            "w",
        ) as handle:
            handle.write(
                text
            )

            handle.flush()

            os.fsync(
                handle.fileno()
            )

        os.replace(
            temp_path,
            path,
        )

    except Exception:
        temp_path.unlink(
            missing_ok=True
        )

        raise


def scan_diagnostics(
    text: str,
    max_examples: int = 2,
) -> tuple[
    Counter[str],
    dict[str, list[str]],
]:
    counts: Counter[str] = Counter()

    examples: dict[
        str,
        list[str],
    ] = defaultdict(list)

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith(
            "State:"
        ):
            continue

        matched_known = False

        for name, pattern in DIAGNOSTIC_PATTERNS:
            if pattern.search(
                line
            ):
                counts[name] += 1
                matched_known = True

                if (
                    len(
                        examples[name]
                    )
                    < max_examples
                    and line
                    not in examples[name]
                ):
                    examples[
                        name
                    ].append(
                        line
                    )

        if (
            not matched_known
            and GENERIC_ERROR_RE.search(
                line
            )
        ):
            counts["ERROR"] += 1

            if (
                len(
                    examples["ERROR"]
                )
                < max_examples
                and line
                not in examples["ERROR"]
            ):
                examples[
                    "ERROR"
                ].append(
                    line
                )

    return (
        counts,
        dict(examples),
    )


def parse_slurm_out(
    path: Path,
    mode: str,
    root: Path,
) -> JobRecord:
    record = JobRecord(
        source_file=path,
        mode=mode,
        root=root,
    )

    try:
        record.subdir = str(
            path.parent.relative_to(
                root
            )
        )

    except ValueError:
        record.subdir = str(
            path.parent
        )

    if mode == "array":
        match = ARRAY_FILENAME_RE.match(
            path.name
        )

        if match:
            record.array_id = (
                match.group(
                    "array_id"
                )
            )

            record.task_id = (
                match.group(
                    "task_id"
                )
            )

    else:
        match = JOB_FILENAME_RE.match(
            path.name
        )

        if match:
            record.job_id = (
                match.group(
                    "job_id"
                )
            )

    try:
        text = path.read_text(
            encoding="utf-8",
            errors="replace",
        )

    except Exception as exc:
        record.state = "READ_ERROR"
        record.parse_note = str(
            exc
        )

        return record

    (
        record.diagnostic_counts,
        record.diagnostic_examples,
    ) = scan_diagnostics(
        text
    )

    running_matches = list(
        RUNNING_DIR_RE.finditer(
            text
        )
    )

    if running_matches:
        match = running_matches[
            -1
        ]

        record.running_dir = (
            match.group(
                "running_dir"
            ).strip()
        )

        record.dirlist_line = (
            match.group(
                "dirlist_line"
            ).strip()
        )

        record.start_time_str = (
            match.group(
                "start_ts"
            ).strip()
        )

        if not record.task_id:
            record.task_id = (
                match.group(
                    "array_task_id"
                ).strip()
            )

    if not record.start_time_str:
        record.start_time_str = extract_last(
            STEP_START_RE,
            text,
            "start_ts",
        )

    record.finish_time_str = extract_last(
        FINISHED_AT_RE,
        text,
        "finish_ts",
    )

    record.job_id = (
        extract_last(
            JOB_ID_RE,
            text,
            "job_id",
        )
        or record.job_id
    )

    record.array_job_id = extract_last(
        ARRAY_JOB_ID_RE,
        text,
        "array_job_id",
    )

    state_matches = list(
        STATE_RE.finditer(
            text
        )
    )

    if state_matches:
        match = state_matches[
            -1
        ]

        record.state = (
            match.group(
                "state"
            ).strip()
        )

        record.state_detail = (
            match.group(
                "detail"
            )
            or ""
        ).strip()

    record.nodes = extract_last(
        NODES_RE,
        text,
        "nodes",
    )

    record.cores_per_node = extract_last(
        CORES_PER_NODE_RE,
        text,
        "cpn",
    )

    record.cpu_utilized = extract_last(
        CPU_UTIL_RE,
        text,
        "cpu_util",
    )

    cpu_matches = list(
        CPU_EFF_RE.finditer(
            text
        )
    )

    if cpu_matches:
        match = cpu_matches[
            -1
        ]

        record.cpu_efficiency = (
            match.group(
                "cpu_eff"
            ).strip()
        )

        record.cpu_core_walltime = (
            match.group(
                "cpu_core_walltime"
            ).strip()
        )

    record.wall_time = extract_last(
        WALL_RE,
        text,
        "wall",
    )

    record.memory_utilized = extract_last(
        MEM_UTIL_RE,
        text,
        "mem_util",
    )

    memory_matches = list(
        MEM_EFF_RE.finditer(
            text
        )
    )

    if memory_matches:
        match = memory_matches[
            -1
        ]

        record.memory_efficiency = (
            match.group(
                "mem_eff"
            ).strip()
        )

        (
            record.requested_memory_total,
            record.requested_memory_per_core,
        ) = parse_requested_memory(
            match.group(
                "mem_req"
            )
        )

    if record.state == "NO REPORT":
        record.parse_note = (
            "No Slurm job report block found"
        )

    return record


def gather_array_files(
    root: Path,
    array_ids: list[str],
) -> list[Path]:
    files = []

    for array_id in array_ids:
        matches = sorted(
            (
                path
                for path
                in root.glob(
                    f"slurm-{array_id}_*.out"
                )
                if ARRAY_FILENAME_RE.match(
                    path.name
                )
            ),
            key=lambda path: natural_key(
                path.name
            ),
        )

        if not matches:
            warn(
                f"no files found for "
                f"array ID {array_id}"
            )

        files.extend(
            matches
        )

    return files


def gather_subdir_files(
    root: Path,
    dir_regex: str,
    recursive: bool,
) -> list[Path]:
    try:
        directory_pattern = re.compile(
            dir_regex
        )

    except re.error as exc:
        raise ValueError(
            f"invalid --dir-regex: {exc}"
        ) from exc

    if recursive:
        iterator = root.rglob(
            "slurm-*.out"
        )

    else:
        iterator = (
            path
            for child in root.iterdir()
            if child.is_dir()
            for path in child.glob(
                "slurm-*.out"
            )
        )

    files = []

    for path in iterator:
        if not path.is_file():
            continue

        if not JOB_FILENAME_RE.match(
            path.name
        ):
            continue

        if not directory_pattern.fullmatch(
            path.parent.name
        ):
            continue

        files.append(
            path
        )

    return sorted(
        files,
        key=lambda path: natural_key(
            str(
                path.relative_to(
                    root
                )
            )
        ),
    )


def group_latest_subdir_runs(
    records: list[JobRecord],
) -> tuple[
    list[JobRecord],
    dict[str, list[JobRecord]],
]:
    groups: dict[
        str,
        list[JobRecord],
    ] = defaultdict(list)

    for record in records:
        groups[
            record.subdir
        ].append(
            record
        )

    latest = []
    older = {}

    for subdir, group in groups.items():
        group.sort(
            key=lambda record: (
                record.numeric_job_id,
                (
                    record.source_file
                    .stat()
                    .st_mtime
                    if record.source_file.exists()
                    else -1
                ),
            ),
            reverse=True,
        )

        latest.append(
            group[0]
        )

        if len(group) > 1:
            older[subdir] = (
                group[1:]
            )

    latest.sort(
        key=lambda record: natural_key(
            record.subdir
        )
    )

    return latest, older


def print_table(
    headers: list[str],
    rows: list[list[str]],
) -> None:
    if not rows:
        print(
            "No entries found."
        )

        return

    widths = [
        len(header)
        for header in headers
    ]

    for row in rows:
        for index, cell in enumerate(
            row
        ):
            widths[index] = max(
                widths[index],
                len(cell),
            )

    def format_row(
        row: list[str],
    ) -> str:
        return "  ".join(
            cell.ljust(
                widths[index]
            )
            for index, cell
            in enumerate(row)
        )

    print(
        format_row(
            headers
        )
    )

    print(
        format_row(
            [
                "-" * width
                for width
                in widths
            ]
        )
    )

    for row in rows:
        print(
            format_row(
                row
            )
        )


def display_subdir(
    record: JobRecord,
    older: dict[
        str,
        list[JobRecord],
    ],
) -> str:
    return (
        record.subdir
        + (
            "*"
            if record.subdir
            in older
            else ""
        )
    )


def state_detail(
    record: JobRecord,
) -> str:
    return (
        record.state_detail
        or record.parse_note
        or "-"
    )


def matches_filters(
    record: JobRecord,
    states: set[str],
    diagnostics: set[str],
) -> bool:
    if (
        states
        and record.state.upper()
        not in states
    ):
        return False

    if diagnostics:
        if not any(
            record.diagnostic_counts.get(
                diagnostic,
                0,
            )
            for diagnostic
            in diagnostics
        ):
            return False

    return True


def select_current_records(
    records: list[JobRecord],
    show_all: bool,
    states: set[str],
    diagnostics: set[str],
) -> list[JobRecord]:
    filtered = [
        record
        for record in records
        if matches_filters(
            record,
            states,
            diagnostics,
        )
    ]

    if (
        show_all
        or states
        or diagnostics
    ):
        return filtered

    return [
        record
        for record in filtered
        if (
            not record.is_completed
            or bool(
                record.diagnostic_counts
            )
        )
    ]


def print_diagnostic_details(
    records: list[JobRecord],
) -> None:
    interesting = [
        record
        for record in records
        if record.diagnostic_counts
    ]

    if not interesting:
        return

    print(
        "\nDiagnostic excerpts:"
    )

    ordered_names = [
        name
        for name, _
        in DIAGNOSTIC_PATTERNS
    ] + [
        "ERROR"
    ]

    for record in interesting:
        print(
            f"\n[{record.subdir}] "
            f"{record.source_file.name}"
        )

        for name in ordered_names:
            count = (
                record.diagnostic_counts.get(
                    name,
                    0,
                )
            )

            if not count:
                continue

            print(
                f"  {name} "
                f"({count} occurrence"
                f"{'s' if count != 1 else ''}):"
            )

            for line in (
                record.diagnostic_examples.get(
                    name,
                    [],
                )
            ):
                print(
                    f"    {line}"
                )


def report_current_status(
    records: list[JobRecord],
    older: dict[
        str,
        list[JobRecord],
    ],
    show_all: bool,
    details: bool,
    states: set[str],
    diagnostics: set[str],
    paths_only: bool,
) -> None:
    selected = select_current_records(
        records,
        show_all,
        states,
        diagnostics,
    )

    if paths_only:
        for record in selected:
            print(
                record.subdir
            )

        return

    if not selected:
        if (
            states
            or diagnostics
        ):
            print(
                "No latest runs match "
                "the requested filters."
            )

        else:
            print(
                "No failures or diagnostics "
                "in the latest run of any "
                "matched subdirectory."
            )

    else:
        headers = [
            "subdir",
            "job_id",
            "state",
            "detail",
            "diagnostics",
            "file",
        ]

        rows = [
            [
                display_subdir(
                    record,
                    older,
                ),
                record.job_id or "?",
                record.state,
                state_detail(
                    record
                ),
                record.diagnostics,
                record.source_file.name,
            ]

            for record
            in selected
        ]

        print_table(
            headers,
            rows,
        )

    if older:
        print(
            "\n* = older Slurm run(s) "
            "exist for this directory; "
            "use --history to show them."
        )

    if details:
        print_diagnostic_details(
            selected
        )


def summarize_header(
    records: list[JobRecord],
) -> list[str]:
    completed = [
        record
        for record in records
        if record.is_completed
    ]

    if not completed:
        return [
            "No completed latest runs found."
        ]

    def unique(
        values: Iterable[str],
    ) -> list[str]:
        return sorted(
            {
                value
                for value in values
                if value
            }
        )

    nodes = unique(
        record.nodes
        for record in completed
    )

    cores = unique(
        record.cores_per_node
        for record in completed
    )

    memory_total = unique(
        record.requested_memory_total
        for record in completed
    )

    memory_core = unique(
        record.requested_memory_per_core
        for record in completed
    )

    lines = [
        "Latest completed-run resource summary:",
        (
            "  Nodes: "
            + (
                " / ".join(nodes)
                if nodes
                else "unknown"
            )
        ),
        (
            "  Cores per node: "
            + (
                " / ".join(cores)
                if cores
                else "unknown"
            )
        ),
        (
            "  Requested memory total: "
            + (
                " / ".join(
                    memory_total
                )
                if memory_total
                else "unknown"
            )
        ),
        (
            "  Requested memory per core: "
            + (
                " / ".join(
                    memory_core
                )
                if memory_core
                else "unknown"
            )
        ),
    ]

    if (
        len(nodes) > 1
        or len(cores) > 1
        or len(memory_total) > 1
        or len(memory_core) > 1
    ):
        lines.append(
            "  Note: mixed resource "
            "settings detected."
        )

    return lines


def summarize_global_performance(
    records: list[JobRecord],
) -> list[str]:
    completed = [
        record
        for record in records
        if record.is_completed
    ]

    if not completed:
        return [
            "No completed latest runs found."
        ]

    lines = [
        (
            "Global performance summary "
            "(latest run per subdirectory only):"
        ),
        (
            f"  Completed jobs: "
            f"{len(completed)} / {len(records)}"
        ),
    ]

    runtimes = [
        record.runtime_seconds
        for record in completed
        if record.runtime_seconds
        is not None
    ]

    cpu_used = [
        value
        for value
        in (
            parse_slurm_duration_to_seconds(
                record.cpu_utilized
            )
            for record
            in completed
        )
        if value is not None
    ]

    cpu_eff = [
        value
        for value
        in (
            parse_percent(
                record.cpu_efficiency
            )
            for record
            in completed
        )
        if value is not None
    ]

    mem_used = [
        value
        for value
        in (
            parse_memory_to_bytes(
                record.memory_utilized
            )
            for record
            in completed
        )
        if value is not None
    ]

    mem_eff = [
        value
        for value
        in (
            parse_percent(
                record.memory_efficiency
            )
            for record
            in completed
        )
        if value is not None
    ]

    if runtimes:
        mean_runtime = (
            sum(runtimes)
            / len(runtimes)
        )

        lines.append(
            "  Mean job runtime: "
            + format_seconds_as_hms(
                mean_runtime
            )
        )

        if mean_runtime > 0:
            lines.append(
                "  Mean jobs/hour from runtime: "
                f"{3600.0 / mean_runtime:.2f}"
            )

    if cpu_used:
        lines.append(
            "  Mean CPU utilized: "
            + format_seconds_as_hms(
                sum(cpu_used)
                / len(cpu_used)
            )
        )

    if cpu_eff:
        lines.append(
            "  Mean CPU efficiency: "
            f"{sum(cpu_eff) / len(cpu_eff):.2f}%"
        )

    if mem_used:
        lines.append(
            "  Mean memory utilized: "
            + format_bytes(
                sum(mem_used)
                / len(mem_used)
            )
        )

    if mem_eff:
        lines.append(
            "  Mean memory efficiency: "
            f"{sum(mem_eff) / len(mem_eff):.2f}%"
        )

    timed = [
        record
        for record in completed
        if (
            record.started_at
            and record.finished_at
        )
    ]

    if timed:
        earliest = min(
            record.started_at
            for record in timed
            if record.started_at
            is not None
        )

        latest = max(
            record.finished_at
            for record in timed
            if record.finished_at
            is not None
        )

        window = (
            latest - earliest
        ).total_seconds()

        if window > 0:
            lines.append(
                "  Observed throughput over "
                "batch window: "
                f"{len(timed) * 3600.0 / window:.2f} "
                "jobs/hour"
            )

            lines.append(
                "  Batch window: "
                f"{earliest:%Y-%m-%d %H:%M:%S}"
                " -> "
                f"{latest:%Y-%m-%d %H:%M:%S}"
            )

    return lines


def report_performance(
    records: list[JobRecord],
    older: Optional[
        dict[
            str,
            list[JobRecord],
        ]
    ] = None,
    mode: str = "subdir",
) -> None:
    older = older or {}

    completed = [
        record
        for record in records
        if record.is_completed
    ]

    if not completed:
        print(
            "No completed latest runs found."
        )

        return

    for line in summarize_header(
        completed
    ):
        print(
            line
        )

    print()

    for line in summarize_global_performance(
        records
    ):
        print(
            line
        )

    print()

    if mode == "array":
        headers = [
            "array_id",
            "task_id",
            "running_dir",
            "job_id",
            "wall_time",
            "cpu_utilized",
            "cpu_eff",
            "mem_utilized",
            "mem_eff",
            "diagnostics",
        ]

        rows = [
            [
                record.array_id or "?",
                record.task_id or "?",
                record.running_dir or "UNKNOWN",
                record.job_id or "?",
                record.wall_time or "?",
                record.cpu_utilized or "?",
                record.cpu_efficiency or "?",
                record.memory_utilized or "?",
                record.memory_efficiency or "?",
                record.diagnostics,
            ]

            for record
            in completed
        ]

    else:
        headers = [
            "subdir",
            "job_id",
            "wall_time",
            "cpu_utilized",
            "cpu_eff",
            "mem_utilized",
            "mem_eff",
            "diagnostics",
        ]

        rows = [
            [
                display_subdir(
                    record,
                    older,
                ),
                record.job_id or "?",
                record.wall_time or "?",
                record.cpu_utilized or "?",
                record.cpu_efficiency or "?",
                record.memory_utilized or "?",
                record.memory_efficiency or "?",
                record.diagnostics,
            ]

            for record
            in completed
        ]

    print_table(
        headers,
        rows,
    )

    if (
        older
        and mode == "subdir"
    ):
        print(
            "\n* = older Slurm run(s) exist "
            "for this directory; "
            "use --history to show them."
        )


def report_history(
    older: dict[
        str,
        list[JobRecord],
    ],
) -> None:
    if not older:
        print(
            "\nNo older runs found."
        )

        return

    print(
        "\nOlder runs "
        "(not used to determine current status/performance):"
    )

    headers = [
        "subdir",
        "job_id",
        "state",
        "detail",
        "diagnostics",
        "file",
    ]

    rows = []

    for subdir in sorted(
        older,
        key=natural_key,
    ):
        for record in older[
            subdir
        ]:
            rows.append(
                [
                    subdir,
                    record.job_id or "?",
                    record.state,
                    state_detail(
                        record
                    ),
                    record.diagnostics,
                    record.source_file.name,
                ]
            )

    print_table(
        headers,
        rows,
    )


def write_failed_subdirs(
    records: list[JobRecord],
    output: Path,
) -> int:
    failed = [
        record.subdir
        for record in records
        if not record.is_completed
    ]

    text = "".join(
        f"{subdir}\n"
        for subdir in failed
    )

    atomic_write(
        output,
        text,
    )

    return len(
        failed
    )


def report_array(
    records: list[JobRecord],
    details: bool,
    states: set[str],
    diagnostics: set[str],
    paths_only: bool,
) -> None:
    filtered = [
        record
        for record in records
        if matches_filters(
            record,
            states,
            diagnostics,
        )
    ]

    if (
        not states
        and not diagnostics
    ):
        filtered = [
            record
            for record in filtered
            if (
                not record.is_completed
                or bool(
                    record.diagnostic_counts
                )
            )
        ]

    if paths_only:
        for record in filtered:
            print(
                record.running_dir
                or str(
                    record.source_file.parent
                )
            )

        return

    if not filtered:
        print(
            "No failed jobs or diagnostics found."
        )

        return

    headers = [
        "array_id",
        "task_id",
        "running_dir",
        "state",
        "detail",
        "diagnostics",
        "file",
    ]

    rows = [
        [
            record.array_id or "?",
            record.task_id or "?",
            record.running_dir or "UNKNOWN",
            record.state,
            state_detail(
                record
            ),
            record.diagnostics,
            record.source_file.name,
        ]

        for record
        in filtered
    ]

    print_table(
        headers,
        rows,
    )

    if details:
        print_diagnostic_details(
            filtered
        )


def add_filter_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    parser.add_argument(
        "--state",
        metavar="STATE[,STATE...]",
        help=(
            "Only show matching Slurm states, "
            "e.g. FAILED,TIMEOUT"
        ),
    )

    parser.add_argument(
        "--diag",
        metavar="TYPE[,TYPE...]",
        help=(
            "Only show matching diagnostics, "
            "e.g. OOM,EAGER_MSG"
        ),
    )

    parser.add_argument(
        "--paths",
        action="store_true",
        help=(
            "Print only matching directories; "
            "useful for scripting"
        ),
    )


def build_subdir_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="slurm_scan",
        description=(
            "Scan Slurm output files in subdirectories. "
            "Only the highest job ID per directory "
            "determines current status."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  slurm_scan\n"
            "  slurm_scan --all\n"
            "  slurm_scan -p\n"
            "  slurm_scan --diag OOM\n"
            "  slurm_scan --diag OOM --paths\n"
            "  slurm_scan --state FAILED,TIMEOUT\n"
            "  slurm_scan --failed\n"
            "  slurm_scan --recursive --details\n"
            "\n"
            "Array jobs:\n"
            "  slurm_scan array ARRAY_ID [ARRAY_ID ...]\n"
        ),
    )

    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        help=(
            "Show every latest run, "
            "including clean COMPLETED jobs"
        ),
    )

    parser.add_argument(
        "-p",
        "--performance",
        action="store_true",
        help=(
            "Show resource and performance "
            "information for completed latest runs"
        ),
    )

    parser.add_argument(
        "--history",
        action="store_true",
        help="Also show older runs",
    )

    parser.add_argument(
        "--details",
        action="store_true",
        help=(
            "Print representative diagnostic lines"
        ),
    )

    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help=(
            "Search all descendant directories"
        ),
    )

    parser.add_argument(
        "--failed",
        "--failed-file",
        nargs="?",
        const="failed_subdirs.txt",
        default=None,
        metavar="FILE",
        help=(
            "Write latest non-COMPLETED "
            "subdirectories to FILE "
            "(default: failed_subdirs.txt)"
        ),
    )

    parser.add_argument(
        "--root",
        default=".",
        help=(
            "Root directory to scan "
            "(default: current directory)"
        ),
    )

    parser.add_argument(
        "--dir-regex",
        default=r".*",
        help=(
            "Regex applied to parent directory name"
        ),
    )

    add_filter_arguments(
        parser
    )

    return parser


def build_array_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="slurm_scan array",
        description=(
            "Scan slurm-ARRAYID_TASKID.out files."
        ),
    )

    parser.add_argument(
        "array_ids",
        nargs="+",
        help="One or more array job IDs",
    )

    parser.add_argument(
        "-p",
        "--performance",
        action="store_true",
    )

    parser.add_argument(
        "--details",
        action="store_true",
    )

    parser.add_argument(
        "--root",
        default=".",
        help=(
            "Directory containing array output files"
        ),
    )

    add_filter_arguments(
        parser
    )

    return parser


def validate_filters(
    diagnostics: set[str],
) -> None:
    unknown = (
        diagnostics
        - KNOWN_DIAGNOSTICS
    )

    if unknown:
        known = ", ".join(
            sorted(
                KNOWN_DIAGNOSTICS
            )
        )

        die(
            "unknown diagnostic type(s): "
            + ", ".join(
                sorted(unknown)
            )
            + f"\nKnown diagnostics: {known}",
            2,
        )


def main_subdir(
    argv: list[str],
) -> int:
    parser = build_subdir_parser()

    args = parser.parse_args(
        argv
    )

    states = parse_filter(
        args.state
    )

    diagnostics = parse_filter(
        args.diag
    )

    validate_filters(
        diagnostics
    )

    if (
        args.paths
        and args.performance
    ):
        parser.error(
            "--paths and --performance "
            "cannot be combined"
        )

    if (
        args.paths
        and args.details
    ):
        parser.error(
            "--paths and --details "
            "cannot be combined"
        )

    if (
        args.paths
        and args.history
    ):
        parser.error(
            "--paths and --history "
            "cannot be combined"
        )

    root = (
        Path(
            args.root
        )
        .expanduser()
        .resolve()
    )

    if not root.is_dir():
        die(
            f"root directory does not exist: "
            f"{root}"
        )

    try:
        files = gather_subdir_files(
            root,
            args.dir_regex,
            args.recursive,
        )

    except ValueError as exc:
        die(
            str(exc),
            2,
        )

    if not files:
        location = (
            "descendant directories"
            if args.recursive
            else "immediate subdirectories"
        )

        die(
            "no valid slurm-JOBID.out files "
            f"found in {location} of {root}"
        )

    all_records = [
        parse_slurm_out(
            path,
            "subdir",
            root,
        )
        for path in files
    ]

    latest, older = (
        group_latest_subdir_runs(
            all_records
        )
    )

    if args.failed is not None:
        failed_path = Path(
            args.failed
        )

        if not failed_path.is_absolute():
            failed_path = (
                root
                / failed_path
            )

        count = write_failed_subdirs(
            latest,
            failed_path,
        )

        if not args.paths:
            print(
                f"Wrote {count} failed "
                f"subdirector"
                f"{'y' if count == 1 else 'ies'} "
                f"to {failed_path}"
            )

            print()

    if args.performance:
        selected = [
            record
            for record in latest
            if matches_filters(
                record,
                states,
                diagnostics,
            )
        ]

        report_performance(
            selected,
            older=older,
            mode="subdir",
        )

        if args.details:
            print_diagnostic_details(
                [
                    record
                    for record in selected
                    if record.is_completed
                ]
            )

    else:
        report_current_status(
            latest,
            older=older,
            show_all=args.all,
            details=args.details,
            states=states,
            diagnostics=diagnostics,
            paths_only=args.paths,
        )

    if args.history:
        report_history(
            older
        )

    return 0


def main_array(
    argv: list[str],
) -> int:
    parser = build_array_parser()

    args = parser.parse_args(
        argv
    )

    states = parse_filter(
        args.state
    )

    diagnostics = parse_filter(
        args.diag
    )

    validate_filters(
        diagnostics
    )

    if (
        args.paths
        and args.performance
    ):
        parser.error(
            "--paths and --performance "
            "cannot be combined"
        )

    if (
        args.paths
        and args.details
    ):
        parser.error(
            "--paths and --details "
            "cannot be combined"
        )

    root = (
        Path(
            args.root
        )
        .expanduser()
        .resolve()
    )

    if not root.is_dir():
        die(
            f"root directory does not exist: "
            f"{root}"
        )

    files = gather_array_files(
        root,
        args.array_ids,
    )

    if not files:
        die(
            "no matching Slurm array "
            "output files found"
        )

    records = [
        parse_slurm_out(
            path,
            "array",
            root,
        )
        for path in files
    ]

    filtered = [
        record
        for record in records
        if matches_filters(
            record,
            states,
            diagnostics,
        )
    ]

    if args.performance:
        report_performance(
            filtered,
            mode="array",
        )

        return 0

    report_array(
        records,
        details=args.details,
        states=states,
        diagnostics=diagnostics,
        paths_only=args.paths,
    )

    return 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv

    if (
        argv
        and argv[0] == "array"
    ):
        return main_array(
            argv[1:]
        )

    return main_subdir(
        argv
    )


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
