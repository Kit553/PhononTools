recentjobs() {
    local d="${1:-1h}"
    local n factor seconds
    local number='([0-9]+([.][0-9]*)?|[.][0-9]+)'

    if [[ "$d" =~ ^${number}(m|min|mins|minute|minutes)$ ]]; then
        n="${BASH_REMATCH[1]}"
        factor=60
    elif [[ "$d" =~ ^${number}(h|hr|hrs|hour|hours)$ ]]; then
        n="${BASH_REMATCH[1]}"
        factor=3600
    elif [[ "$d" =~ ^${number}(d|day|days)$ ]]; then
        n="${BASH_REMATCH[1]}"
        factor=86400
    else
        echo "Usage: recentjobs <duration>"
        echo "Examples: recentjobs 10min | recentjobs 0.5h | recentjobs 24h | recentjobs 1.25d"
        return 1
    fi

    # Slurm accepts integer time offsets; round fractional durations to seconds.
    seconds=$(LC_ALL=C awk -v n="$n" -v factor="$factor" 'BEGIN { printf "%.0f", int(n * factor + 0.5) }')

    sacct -X -n -P \
        -S "now-${seconds}seconds" \
        -E now \
        -s COMPLETED,FAILED,CANCELLED,TIMEOUT,OUT_OF_MEMORY,NODE_FAIL,PREEMPTED \
        -o JobID,State,Elapsed,WorkDir%500 |
    awk -F'|' -v user="$USER" '
        BEGIN {
            printf "%-14s %-18s %-12s %s\n", "JOBID", "STATE", "RUNTIME", "WORKDIR"
            printf "%-14s %-18s %-12s %s\n", "-----", "-----", "-------", "-------"
        }
        {
            sub("^/scratch/tmp/" user "/?", "", $4)

            if ($2 ~ /^CANCELLED/)
                $2 = "CANCELLED by user"

            printf "%-14s %-18s %-12s %s\n", $1, $2, $3, $4
        }
    '
}
