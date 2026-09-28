jp() {
    local jobid="$1"

    if [[ -z "$jobid" ]]; then
        echo "Usage: jp JOBID"
        return 1
    fi

    sstat -j "$jobid" \
        --format=JobID,AveCPU,AveRSS,MaxRSS,AveVMSize,MaxVMSize

    echo
    scontrol show job "$jobid" \
        | grep -oE 'AllocTRES=[^ ]+' \
        | sed 's/^AllocTRES=/Resources: /'

    echo
    date
}