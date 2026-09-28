jq() {
    local dir
    dir=$(realpath "$PWD")

    local jobs
    jobs=$(squeue -u "$USER" -h -t R,PD \
        -o "%i|%j|%t|%M|%l|%V|%R|%S|%Q|%Z" \
        | awk -F'|' -v d="$dir" '$10 == d')

    if [[ -z "$jobs" ]]; then
        echo "NO"
        return 1
    fi

    echo "YES"
    printf "%10s | %25s | %2s | %10s | %10s | %19s | %25s | %19s | %8s | %s\n" \
        "JOBID" "NAME" "ST" "TIME" "TIME_LIMIT" "SUBMIT_TIME" \
        "NODELIST(REASON)" "START_TIME" "PRIORITY" "WORK_DIR"

    while IFS='|' read -r jobid name st time limit submit nodes start priority workdir; do
        printf "%10s | %25s | %2s | %10s | %10s | %19s | %25s | %19s | %8s | %s\n" \
            "$jobid" "$name" "$st" "$time" "$limit" "$submit" \
            "$nodes" "$start" "$priority" "$workdir"
    done <<< "$jobs"
}

