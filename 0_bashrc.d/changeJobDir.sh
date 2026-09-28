cdj() {
    local job_id="$1"
    local job_dir=""

    if [[ -z "$job_id" ]]; then
        echo "Usage: cdj <Job-ID>"
        return 2
    fi

    job_dir=$(
        scontrol show job -o "$job_id" 2>/dev/null |
        sed -n 's/.* WorkDir=\([^ ]*\).*/\1/p'
    )

    if [[ -z "$job_dir" ]]; then
        job_dir=$(
            sacct -j "$job_id" -X -n -P -o WorkDir 2>/dev/null |
            awk -F '|' '$1 != "" && $1 != "Unknown" {print $1; exit}'
        )
    fi

    if [[ -z "$job_dir" ]]; then
        echo "cdj: could not find a working directory for job $job_id"
        return 1
    fi

    if [[ ! -d "$job_dir" ]]; then
        echo "cdj: directory no longer exists: $job_dir"
        return 1
    fi

    cd -- "$job_dir" || return
    printf 'Job %s â†’ %s\n' "$job_id" "$PWD"
}
