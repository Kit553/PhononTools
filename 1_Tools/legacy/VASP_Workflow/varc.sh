#!/usr/bin/env bash

# Archive the current output of a VASP run and prepare the next restart. Timed-out or otherwise interrupted calculations are accepted as long as a usable CONTCAR exists.

set -Eeuo pipefail
shopt -s nullglob

usage() {
    cat <<'EOF'
Usage: varc.sh [OPTIONS] [DIRECTORY]

Archive the current VASP files as *.sN and copy CONTCAR to POSCAR.
DIRECTORY defaults to the current working directory.
If all current output files match an existing .sN backup, exit successfully
without creating files or changing POSCAR (also after --backup-only).
Stop VASP before using this script; it does not check the scheduler.

Options:
  --finished-only  Refuse to restart unless OUTCAR shows normal termination
  --backup-only    Archive the files but do not replace POSCAR
  --force          Accept a non-empty CONTCAR even if it looks unusually short
  -h, --help       Show this help

Files archived when present:
  OUTCAR, CONTCAR, vasprun.xml, out_VASP, out
EOF
}

finished_only=false
backup_only=false
force=false
run_dir=""

while (($#)); do
    case "$1" in
        --finished-only)
            finished_only=true
            ;;
        --backup-only)
            backup_only=true
            ;;
        --force)
            force=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            if (($# > 1)); then
                echo "Error: only one directory may be specified." >&2
                exit 2
            fi
            run_dir="${1:-}"
            break
            ;;
        -*)
            echo "Error: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *)
            if [[ -n "$run_dir" ]]; then
                echo "Error: only one directory may be specified." >&2
                exit 2
            fi
            run_dir="$1"
            ;;
    esac
    shift
done

run_dir="${run_dir:-$PWD}"
if [[ ! -d "$run_dir" ]]; then
    echo "Error: directory does not exist: $run_dir" >&2
    exit 1
fi

cd -- "$run_dir"
run_dir=$PWD
echo "Inspecting VASP calculation in: $run_dir"

if [[ -f OUTCAR ]] && grep -q 'Voluntary context switches' OUTCAR; then
    echo "Status: VASP terminated normally."
    normally_finished=true
else
    echo "Status: VASP did not terminate normally (for example, timeout/cancel)."
    normally_finished=false
fi

if [[ "$finished_only" == true && "$normally_finished" == false ]]; then
    echo "Error: refusing interrupted run because --finished-only was requested." >&2
    exit 1
fi

if [[ ! -s CONTCAR ]]; then
    echo "Error: CONTCAR is missing or empty; no safe restart is possible." >&2
    exit 1
fi

contcar_lines=$(wc -l < CONTCAR)
if ((contcar_lines < 8)) && [[ "$force" == false ]]; then
    echo "Error: CONTCAR has only $contcar_lines lines and may be truncated." >&2
    echo "Inspect it manually, then use --force if it is genuinely valid." >&2
    exit 1
fi

next_id=1
archive_names=(OUTCAR CONTCAR vasprun.xml out_VASP out)
for name in "${archive_names[@]}"; do
    for archived in "$name".s*; do
        suffix=${archived#"$name.s"}
        if [[ "$suffix" =~ ^[0-9]+$ ]]; then
            suffix_value=$((10#$suffix))
            if ((suffix_value >= next_id)); then
                next_id=$((suffix_value + 1))
            fi
        fi
    done
done

for candidate in CONTCAR.s*; do
    suffix=${candidate#CONTCAR.s}
    [[ "$suffix" =~ ^[0-9]+$ ]] || continue
    identical=true
    for name in "${archive_names[@]}"; do
        if [[ -f "$name" ]]; then
            if [[ ! -f "$name.s$suffix" ]] || ! cmp -s -- "$name" "$name.s$suffix"; then
                identical=false
                break
            fi
        fi
    done
    if [[ "$identical" == true ]]; then
        echo "Already backed up as s$suffix: all current output files are identical."
        echo "Nothing changed; POSCAR left untouched."
        if [[ "$backup_only" == false ]] && ! cmp -s -- CONTCAR POSCAR; then
            echo "Note: POSCAR differs from CONTCAR. To prepare this archived structure explicitly:"
            echo "  cp -- CONTCAR POSCAR"
        fi
        exit 0
    fi
done

if [[ ! -e POSCAR_init ]]; then
    if [[ -f POSCAR ]]; then
        cp -p -- POSCAR POSCAR_init
        echo "Saved original structure as POSCAR_init."
    else
        echo "Warning: POSCAR is missing, so POSCAR_init could not be created." >&2
    fi
fi

echo "Archiving current run as s$next_id:"
archived_count=0
for name in "${archive_names[@]}"; do
    if [[ -f "$name" ]]; then
        cp -p -- "$name" "$name.s$next_id"
        echo "  $name -> $name.s$next_id"
        ((archived_count += 1))
    fi
done

if ((archived_count == 0)); then
    echo "Error: none of the expected VASP files could be archived." >&2
    exit 1
fi

if [[ "$backup_only" == false ]]; then
    temporary_poscar=$(mktemp .POSCAR.restart.XXXXXX)
    trap 'rm -f -- "${temporary_poscar:-}"' EXIT
    cp -p -- CONTCAR "$temporary_poscar"
    mv -f -- "$temporary_poscar" POSCAR
    trap - EXIT
    echo "Restart prepared: CONTCAR -> POSCAR"
else
    echo "Backup complete; POSCAR unchanged (--backup-only)."
fi

echo "Finished."