#!/usr/bin/env bash
set -euo pipefail

PREFIX=""
VERBOSE=false

while (( $# )); do
    case "$1" in
        --prefix) PREFIX="$2"; shift 2 ;;
        -v|--verbose) VERBOSE=true; shift ;;
        -h|--help)
            echo "Usage: ./ConvCheck --prefix PREFIX [-v|--verbose]"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

[[ -n "$PREFIX" ]] || {
    echo "Usage: ./ConvCheck --prefix PREFIX [-v|--verbose]"
    exit 1
}

printf "%-24s %10s %16s %16s %18s %18s %18s %16s\n" \
    "Directory" "ENCUT" "ETOT" "ETOT/Atom" \
    "MaxForceAtom" "MaxCellStress" "PullayStress" "WallTime (s)"

printf '%*s\n' 142 '' | tr ' ' '-'

for dir in ${PREFIX}*/; do
    [[ -d "$dir" ]] || continue

    out="${dir}OUTCAR"

    if [[ ! -f "$out" ]]; then
        $VERBOSE && echo "Skipping $dir: no OUTCAR" >&2
        continue
    fi

    if ! grep -q "General timing and accounting informations" "$out"; then
        $VERBOSE && echo "Warning: $dir may be unfinished" >&2
    fi

    ENCUT=$(awk '/ENCUT[[:space:]]*=/{print $3; exit}' "$out")

    ETOT=$(awk '/free  energy   TOTEN/{e=$5} END{print e}' "$out")

    NIONS=$(awk '/NIONS[[:space:]]*=/{n=$12} END{print n}' "$out")

    if [[ -n "$ETOT" && -n "$NIONS" ]]; then
        EPA=$(awk -v e="$ETOT" -v n="$NIONS" 'BEGIN{printf "%.8f",e/n}')
    else
        EPA="NA"
    fi

    MAXFORCE=$(awk '
        /TOTAL-FORCE/ {
            active=1
            max=0
            next
        }

        active && NF>=6 && $1 ~ /^[-+0-9.]/ {
            f=sqrt($4*$4 + $5*$5 + $6*$6)
            if (f>max) max=f
            last=max
        }

        active && NF==0 {
            active=0
        }

        END {
            if (last!="") printf "%.6f",last
            else print "NA"
        }
    ' "$out")

    MAXCELL=$(awk '
        /in kB/ {
            max=0
            for(i=3;i<=8;i++) {
                x=$i
                if(x<0) x=-x
                if(x>max) max=x
            }
            last=max
        }

        END {
            if(last!="") printf "%.4f",last
            else print "NA"
        }
    ' "$out")

    PULAY=$(awk '
        /Pullay stress/ {
            for(i=1;i<=NF;i++) {
                if($i=="stress") {
                    x=$(i+2)
                    if(x<0) x=-x
                    last=x
                }
            }
        }

        END {
            if(last!="") printf "%.4f",last
            else print "NA"
        }
    ' "$out")

    ELAPSED=$(awk '/Elapsed time \(sec\)/ {t=$4} END{print t}' "$out")
    if [[ -n "$ELAPSED" ]]; then
        WALLTIME=$(awk -v s="$ELAPSED" '
            BEGIN {
                h=int(s/3600)
                m=int((s-h*3600)/60)
                sec=int(s-h*3600-m*60)
                printf "%02d:%02d:%02d", h, m, sec
            }
        ')
    else
        WALLTIME="NA"
    fi

    printf "%-24s %10s %16s %16s %18s %18s %18s %12s\n" \
        "${dir%/}" \
        "${ENCUT:-NA}" \
        "${ETOT:-NA}" \
        "$EPA" \
        "$MAXFORCE" \
        "$MAXCELL" \
        "$PULAY"    \
        "$ELAPSED"

done
