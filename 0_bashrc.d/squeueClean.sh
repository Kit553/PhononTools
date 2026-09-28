sq() {
squeue -o "%.25i | %.25j | %.2t | %.10M | %.10l | %.16V | %.25R | %.16S | %.8p | %60Z" --me | \
awk -F "|" -v user="$USER" 'NR==1 { gsub(/[ \t]+$/, "", $10); 
		    print $1 "|" $2 "|" $3 "|" $4 "|" $5 "|" $6 "|" $7 "|" $8 "|" $9 "|" $10 }
            NR>1 {
              gsub(/^ *[0-9]{4}-/, "      ", $6);
              gsub("^ */scratch/tmp/" user "/", "", $10);
              gsub(/[ \t]+$/, "", $10);
              print $1 "|" $2 "|" $3 "|" $4 "|" $5 "|" $6 "|" $7 "|" $8 "|" $9 "|" $10
            }'
printf '%*s\n' "$COLUMNS" '' | tr ' ' '-'
job_numR=$(( $(squeue --me --states=RUNNING | wc -l) - 1 ))
job_numP=$(( $(squeue --me --states=PENDING | wc -l) - 1 ))
echo "Number of running jobs is $job_numR, $job_numP jobs are pending"
}
