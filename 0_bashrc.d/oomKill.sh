ohmu() {
    local mode="report"
    local root="."

    case "${1:-}" in
        -c|--cancel)
            mode="cancel"
            shift
            ;;
        -r|--restart)
            mode="restart"
            shift
            ;;
        -h|--help)
            cat <<'EOF'
Usage:
  oomjobs [DIR]              Find active jobs whose slurm output contains OOM
  oomjobs -c [DIR]           Cancel matching jobs
  oomjobs -r [DIR]           Cancel and restart matching jobs
  oomjobs -h                 Show this help

Restarted jobs are submitted from the directory containing their slurm output,
NOT from the directory where oomjobs was called.
EOF
            return 0
            ;;
    esac

    [[ $# -gt 0 ]] && root="$1"

    printf "%-12s | %-10s | %-10s | %s\n" \
        "JOBID" "STATE" "ACTION" "DIRECTORY"
    printf '%*s\n' 110 '' | tr ' ' '-'

    declare -A seen

    while IFS= read -r -d '' file; do

        grep -Eqi \
            'oom[_ -]?kill|out of memory|exceeded.*memory|memory limit|Killed process' \
            "$file" || continue

        local jobid
        jobid=$(basename "$file" |
            sed -nE 's/^slurm-([0-9]+)(_[0-9]+)?\.out$/\1\2/p')

        [[ -z "$jobid" ]] && continue

        [[ -n "${seen[$jobid]:-}" ]] && continue
        seen[$jobid]=1

        local info id state
        info=$(squeue -h -j "$jobid" -o "%i|%T" 2>/dev/null | head -n1)

        [[ -z "$info" ]] && continue

        id=${info%%|*}
        state=${info#*|}

        [[ "$state" != "RUNNING" && "$state" != "PENDING" ]] && continue

        local jobdir
        jobdir=$(dirname "$(realpath "$file")")

        case "$mode" in
            report)
                printf "%-12s | %-10s | %-10s | %s\n" \
                    "$id" "$state" "-" "$jobdir"
                ;;

            cancel)
                if scancel "$id"; then
                    printf "%-12s | %-10s | %-10s | %s\n" \
                        "$id" "$state" "CANCELLED" "$jobdir"
                else
                    printf "%-12s | %-10s | %-10s | %s\n" \
                        "$id" "$state" "FAILED" "$jobdir"
                fi
                ;;

            restart)
                local command script
                command=$(scontrol show job -o "$id" 2>/dev/null |
                    sed -n 's/.* Command=\([^ ]*\).*/\1/p')

                if [[ -z "$command" ]]; then
                    printf "%-12s | %-10s | %-10s | %s\n" \
                        "$id" "$state" "NO SCRIPT" "$jobdir"
                    continue
                fi

                script=$(basename "$command")

                if [[ ! -f "$jobdir/$script" ]]; then
                    printf "%-12s | %-10s | %-10s | %s\n" \
                        "$id" "$state" "NO SCRIPT" "$jobdir/$script"
                    continue
                fi

                if ! scancel "$id"; then
                    printf "%-12s | %-10s | %-10s | %s\n" \
                        "$id" "$state" "CANCEL FAIL" "$jobdir"
                    continue
                fi

                local submit_output
                if submit_output=$(cd "$jobdir" && sbatch "./$script" 2>&1); then
                    printf "%-12s | %-10s | %-10s | %s\n" \
                        "$id" "$state" "RESTARTED" "$jobdir"
                    printf "    -> %s\n" "$submit_output"
                else
                    printf "%-12s | %-10s | %-10s | %s\n" \
                        "$id" "$state" "SUBMIT FAIL" "$jobdir"
                    printf "    -> %s\n" "$submit_output"
                fi
                ;;
        esac

    done < <(find "$root" -type f -name 'slurm*.out' -print0 2>/dev/null)
}
