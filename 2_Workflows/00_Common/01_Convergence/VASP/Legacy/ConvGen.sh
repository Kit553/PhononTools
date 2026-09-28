#!/usr/bin/env bash
set -euo pipefail

AXIS=""
MODE=""
CALC="static"
PREFIX=""
SUBMIT=true
VERBOSE=false

VALUES=()
SETS=()

usage() {
    cat <<'TXT'
Usage:
  ./ConvGen.sh AXIS (--man VALUES... | --range LOW HIGH STEP) [OPTIONS]

Axes:
  --encut      Plane-wave cutoff
  --kpoints    Reciprocal k-point density
  --ediff      Electronic SCF threshold
  --ediffg     Ionic force threshold
  --sigma      Smearing width
  --nbands     Number of bands

Value specification:
  --man V1 V2 V3 ...
  --range LOW HIGH STEP

ENCUT convention:
  --encut --range   values are multipliers of max(ENMAX)
  --encut --man     values are absolute ENCUT values in eV

KPOINTS convention:
  values are reciprocal-density parameters (e.g. 50)

Options:
  --prefix PREFIX
      Override the automatic output directory prefix.

      Example:
          --prefix RelaxENCUT_

  --calc static|ions|full|keep

      static   IBRION=-1, NSW=0,   ISIF=2
      ions     IBRION=2,  NSW=100, ISIF=2
      full     IBRION=2,  NSW=100, ISIF=3
      keep     Preserve relaxation tags from the supplied INCAR

  --set TAG=VALUE
      Override an arbitrary INCAR tag. May be repeated.

  --no-submit
      Generate directories but do not submit them.

  -v, --verbose
      Detailed setup/submission information.

  -h, --help
TXT
}

die() {
    echo "Error: $*" >&2
    exit 1
}

is_known_option() {
    case "$1" in
        --encut|--kpoints|--ediff|--ediffg|--sigma|--nbands|\
        --man|--range|--prefix|--calc|--set|--no-submit|\
        -v|--verbose|-h|--help)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

set_tag() {
    local file="$1"
    local tag="$2"
    local value="$3"

    [[ "$tag" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] \
        || die "Invalid INCAR tag '$tag'"

    sed -i "/^[[:space:]]*${tag}[[:space:]]*=/Id" "$file"
    printf '%-10s = %s\n' "$tag" "$value" >> "$file"
}

make_range() {
    awk -v lo="$1" -v hi="$2" -v step="$3" '
    BEGIN {
        if (step == 0)
            exit 1

        if (lo < hi && step < 0)
            exit 1

        if (lo > hi && step > 0)
            exit 1

        for (
            x = lo;
            (step > 0 ? x <= hi + 1e-12 : x >= hi - 1e-12);
            x += step
        )
            printf "%.12g\n", x
    }'
}

get_enmax() {
    awk '
    /ENMAX/ {
        for (i = 1; i <= NF; i++) {
            if ($i == "=") {
                x = $(i+1)
                gsub(/;/, "", x)

                if (x > max)
                    max = x
            }
        }
    }

    END {
        if (max != "")
            print max
    }
    ' POTCAR
}

# ------------------------------------------------------------
# Reciprocal-density -> Gamma-centred k-point mesh.
#
# Supports:
#   - one positive POSCAR scale
#   - one negative target-volume scale
#   - three positive Cartesian scale factors
#
# ------------------------------------------------------------

get_mesh() {
    local density="$1"

    awk -v D="$density" '
    function abs(x) {
        return x < 0 ? -x : x
    }

    NR == 2 {
        ns = NF
        s1 = $1
        s2 = $2
        s3 = $3
    }

    NR >= 3 && NR <= 5 {
        i = NR - 2

        x[i] = $1
        y[i] = $2
        z[i] = $3
    }

    END {

        # ----------------------------------------------------
        # POSCAR scaling
        # ----------------------------------------------------

        if (ns == 1) {

            if (s1 == 0) {
                print "POSCAR scale factor cannot be zero" \
                    > "/dev/stderr"
                exit 1
            }

            rawV = \
                x[1] * (y[2]*z[3] - z[2]*y[3]) \
              - y[1] * (x[2]*z[3] - z[2]*x[3]) \
              + z[1] * (x[2]*y[3] - y[2]*x[3])

            if (s1 > 0)
                scale = s1
            else
                scale = ((-s1) / abs(rawV))^(1.0/3.0)

            ax = x[1] * scale
            ay = y[1] * scale
            az = z[1] * scale

            bx = x[2] * scale
            by = y[2] * scale
            bz = z[2] * scale

            cx = x[3] * scale
            cy = y[3] * scale
            cz = z[3] * scale
        }

        else if (ns == 3) {

            if (s1 <= 0 || s2 <= 0 || s3 <= 0) {
                print "Three POSCAR scale factors must be positive" \
                    > "/dev/stderr"
                exit 1
            }

            ax = x[1] * s1
            ay = y[1] * s2
            az = z[1] * s3

            bx = x[2] * s1
            by = y[2] * s2
            bz = z[2] * s3

            cx = x[3] * s1
            cy = y[3] * s2
            cz = z[3] * s3
        }

        else {
            print "Invalid POSCAR scale line" > "/dev/stderr"
            exit 1
        }

        # ----------------------------------------------------
        # Cell volume
        # ----------------------------------------------------

        V = \
            ax * (by*cz - bz*cy) \
          - ay * (bx*cz - bz*cx) \
          + az * (bx*cy - by*cx)

        if (abs(V) < 1e-12) {
            print "POSCAR lattice has zero volume" \
                > "/dev/stderr"
            exit 1
        }

        # ----------------------------------------------------
        # Reciprocal vectors without 2*pi
        # ----------------------------------------------------

        r1x = (by*cz - bz*cy) / V
        r1y = (bz*cx - bx*cz) / V
        r1z = (bx*cy - by*cx) / V

        r2x = (cy*az - cz*ay) / V
        r2y = (cz*ax - cx*az) / V
        r2z = (cx*ay - cy*ax) / V

        r3x = (ay*bz - az*by) / V
        r3y = (az*bx - ax*bz) / V
        r3z = (ax*by - ay*bx) / V

        b1 = sqrt(r1x*r1x + r1y*r1y + r1z*r1z)
        b2 = sqrt(r2x*r2x + r2y*r2y + r2z*r2z)
        b3 = sqrt(r3x*r3x + r3y*r3y + r3z*r3z)

        n1 = int(D*b1 + 0.5)
        n2 = int(D*b2 + 0.5)
        n3 = int(D*b3 + 0.5)

        if (n1 < 1) n1 = 1
        if (n2 < 1) n2 = 1
        if (n3 < 1) n3 = 1

        print n1, n2, n3
    }
    ' POSCAR
}

patch_jobname() {
    local file="$1"
    local name="$2"

    name=$(
        printf '%s' "$name" \
        | tr ' /' '__' \
        | cut -c1-100
    )

    if grep -qiE \
        '^#SBATCH[[:space:]]+--job-name=' \
        "$file"
    then

        sed -i -E \
            "s|^#SBATCH[[:space:]]+--job-name=.*|#SBATCH --job-name=${name}|I" \
            "$file"

    elif grep -qiE \
        '^#SBATCH[[:space:]]+-J([[:space:]]+|=)' \
        "$file"
    then

        sed -i -E \
            "s|^#SBATCH[[:space:]]+-J([[:space:]]+|=).*|#SBATCH -J ${name}|I" \
            "$file"
    fi
}

write_meta() {
    local dir="$1"
    local scan="$2"
    local effective="$3"

    cat > "$dir/.convgen_meta" <<EOF
AXIS=$AXIS
SCAN_VALUE=$scan
EFFECTIVE_VALUE=$effective
CALC=$CALC
PREFIX=$PREFIX
EOF
}

# ============================================================
# Argument parsing
# ============================================================

while (( $# )); do

    case "$1" in

        --encut|--kpoints|--ediff|--ediffg|--sigma|--nbands)

            [[ -z "$AXIS" ]] \
                || die "Only one convergence axis may be used"

            AXIS="${1#--}"

            shift
            ;;

        --man)

            [[ -z "$MODE" ]] \
                || die "Use only one of --man or --range"

            MODE="man"

            shift

            while (( $# )) && ! is_known_option "$1"; do
                VALUES+=("$1")
                shift
            done
            ;;

        --range)

            [[ -z "$MODE" ]] \
                || die "Use only one of --man or --range"

            (( $# >= 4 )) \
                || die "--range requires LOW HIGH STEP"

            MODE="range"

            LOW="$2"
            HIGH="$3"
            STEP="$4"

            shift 4
            ;;

        --prefix)

            (( $# >= 2 )) \
                || die "--prefix requires a value"

            PREFIX="$2"

            shift 2
            ;;

        --calc)

            (( $# >= 2 )) \
                || die "--calc requires static|ions|full|keep"

            CALC="$2"

            shift 2
            ;;

        --set)

            (( $# >= 2 )) \
                || die "--set requires TAG=VALUE"

            SETS+=("$2")

            shift 2
            ;;

        --no-submit)

            SUBMIT=false
            shift
            ;;

        -v|--verbose)

            VERBOSE=true
            shift
            ;;

        -h|--help)

            usage
            exit 0
            ;;

        *)

            die "Unknown argument: $1"
            ;;
    esac
done

[[ -n "$AXIS" ]] \
    || die "No convergence axis specified"

[[ -n "$MODE" ]] \
    || die "Use either --man or --range"

case "$CALC" in
    static|ions|full|keep)
        ;;
    *)
        die "Unknown --calc mode '$CALC'"
        ;;
esac

if [[ "$AXIS" == "ediffg" && "$CALC" == "static" ]]; then
    die "EDIFFG has no effect in a static run; use --calc ions, --calc full, or --calc keep"
fi

# ============================================================
# Required input files
# ============================================================

for f in INCAR POSCAR POTCAR submit.sh; do

    [[ -f "$f" ]] \
        || die "$f not found"

done

if [[ "$AXIS" != "kpoints" ]]; then

    [[ -f KPOINTS ]] \
        || die "KPOINTS not found"

fi

if [[ "$SUBMIT" == true ]]; then

    command -v sbatch >/dev/null 2>&1 \
        || die "sbatch not found; use --no-submit to generate only"

fi

# ============================================================
# Generate values
# ============================================================

if [[ "$MODE" == "range" ]]; then

    mapfile -t VALUES < <(
        make_range "$LOW" "$HIGH" "$STEP"
    ) || die "Invalid range"

fi

(( ${#VALUES[@]} > 0 )) \
    || die "No convergence values supplied"

ENMAX=$(get_enmax)

[[ -n "$ENMAX" ]] \
    || die "Could not extract ENMAX from POTCAR"

# ============================================================
# Default prefix
# ============================================================

if [[ -z "$PREFIX" ]]; then

    case "$AXIS" in
        encut)
            PREFIX="ENCUT_"
            ;;
        kpoints)
            PREFIX="KPOINTS_"
            ;;
        ediff)
            PREFIX="EDIFF_"
            ;;
        ediffg)
            PREFIX="EDIFFG_"
            ;;
        sigma)
            PREFIX="SIGMA_"
            ;;
        nbands)
            PREFIX="NBANDS_"
            ;;
    esac

fi

# ============================================================
# Build calculation definitions before touching disk
# ============================================================

DIRS=()
DISPLAY=()
EFFECTIVE=()

for value in "${VALUES[@]}"; do

    case "$AXIS" in

        encut)

            if [[ "$MODE" == "range" ]]; then

                encut=$(
                    awk \
                    -v e="$ENMAX" \
                    -v m="$value" '
                    BEGIN {
                        x = e*m
                        i = int(x)

                        if (x > i + 1e-10)
                            print i+1
                        else
                            print i
                    }
                    '
                )

                suffix="${value}x_${encut}"
                display="${value}x"

            else

                encut="$value"
                suffix="$encut"
                display="$value"

            fi

            effective="$encut"
            ;;

        kpoints)

            read -r n1 n2 n3 <<< "$(get_mesh "$value")"

            mesh="${n1}x${n2}x${n3}"

            suffix="D${value}_${mesh}"
            display="$value"
            effective="$mesh"
            ;;

        ediff)

            suffix="$value"
            display="$value"
            effective="$value"
            ;;

        ediffg)

            suffix="$value"
            display="$value"
            effective="$value"
            ;;

        sigma)

            suffix="$value"
            display="$value"
            effective="$value"
            ;;

        nbands)

            [[ "$value" =~ ^[0-9]+$ ]] \
                || die "NBANDS must be integer: '$value'"

            suffix="$value"
            display="$value"
            effective="$value"
            ;;
    esac

    dir="${PREFIX}${suffix}"

    [[ ! -e "$dir" ]] \
        || die "Directory already exists: $dir"

    DIRS+=("$dir")
    DISPLAY+=("$display")
    EFFECTIVE+=("$effective")

done

# ------------------------------------------------------------
# Catch duplicate names before creating anything
# ------------------------------------------------------------

declare -A seen_dirs

for dir in "${DIRS[@]}"; do

    [[ -z "${seen_dirs[$dir]:-}" ]] \
        || die "Duplicate output directory: $dir"

    seen_dirs[$dir]=1

done

# ============================================================
# Create and submit
# ============================================================

JOBIDS=()

for i in "${!VALUES[@]}"; do

    value="${VALUES[$i]}"
    dir="${DIRS[$i]}"
    effective="${EFFECTIVE[$i]}"

    $VERBOSE && echo "Creating $dir"

    mkdir "$dir"

    cp \
        INCAR \
        POSCAR \
        POTCAR \
        submit.sh \
        "$dir/"

    if [[ "$AXIS" != "kpoints" ]]; then
        cp KPOINTS "$dir/"
    fi

    # --------------------------------------------------------
    # Calculation type
    # --------------------------------------------------------

    case "$CALC" in

        static)

            set_tag "$dir/INCAR" IBRION -1
            set_tag "$dir/INCAR" NSW 0
            set_tag "$dir/INCAR" ISIF 2
            ;;

        ions)

            set_tag "$dir/INCAR" IBRION 2
            set_tag "$dir/INCAR" NSW 100
            set_tag "$dir/INCAR" ISIF 2
            ;;

        full)

            set_tag "$dir/INCAR" IBRION 2
            set_tag "$dir/INCAR" NSW 100
            set_tag "$dir/INCAR" ISIF 3
            ;;

        keep)

            ;;
    esac

    # --------------------------------------------------------
    # Generic user overrides
    # --------------------------------------------------------

    for assignment in "${SETS[@]}"; do

        [[ "$assignment" == *=* ]] \
            || die "Invalid --set '$assignment'; expected TAG=VALUE"

        tag="${assignment%%=*}"
        val="${assignment#*=}"

        set_tag \
            "$dir/INCAR" \
            "$tag" \
            "$val"

    done

    # --------------------------------------------------------
    # Scan parameter
    #
    # Applied last so the scan axis cannot accidentally be overwritten by --set.
    # --------------------------------------------------------

    case "$AXIS" in

        encut)

            set_tag \
                "$dir/INCAR" \
                ENCUT \
                "$effective"
            ;;

        kpoints)

            IFS=x read -r n1 n2 n3 <<< "$effective"

            cat > "$dir/KPOINTS" <<EOF
Automatic mesh, density=${value}
0
Gamma
${n1} ${n2} ${n3}
0 0 0
EOF
            ;;

        ediff)

            set_tag \
                "$dir/INCAR" \
                EDIFF \
                "$value"
            ;;

        ediffg)

            set_tag \
                "$dir/INCAR" \
                EDIFFG \
                "$value"
            ;;

        sigma)

            set_tag \
                "$dir/INCAR" \
                SIGMA \
                "$value"
            ;;

        nbands)

            set_tag \
                "$dir/INCAR" \
                NBANDS \
                "$value"
            ;;
    esac

    patch_jobname \
        "$dir/submit.sh" \
        "${dir%/}"

    write_meta \
        "$dir" \
        "$value" \
        "$effective"

    # --------------------------------------------------------
    # Submission
    # --------------------------------------------------------

    if [[ "$SUBMIT" == true ]]; then

        raw_jobid=$(
            cd "$dir"
            sbatch --parsable submit.sh
        )

        jobid="${raw_jobid%%;*}"

        printf '%s\n' \
            "$jobid" \
            > "$dir/.convgen_jobid"

        printf 'JOBID=%s\n' \
            "$jobid" \
            >> "$dir/.convgen_meta"

        JOBIDS+=("$jobid")

        $VERBOSE \
            && echo "  submitted as job $jobid"

    else

        JOBIDS+=("-")

    fi

done

# ============================================================
# Summary
# ============================================================

echo
echo "Convergence scan: ${AXIS^^}"

if [[ "$AXIS" == "encut" ]]; then
    echo "Detected max ENMAX: ${ENMAX} eV"
fi

echo

printf \
    '%-34s | %-14s | %-18s | %-12s\n' \
    "Directory" \
    "Scan value" \
    "Effective value" \
    "JobID"

printf \
    '%-34s-+-%-14s-+-%-18s-+-%-12s\n' \
    "----------------------------------" \
    "--------------" \
    "------------------" \
    "------------"

for i in "${!DIRS[@]}"; do

    printf \
        '%-34s | %-14s | %-18s | %-12s\n' \
        "${DIRS[$i]}" \
        "${DISPLAY[$i]}" \
        "${EFFECTIVE[$i]}" \
        "${JOBIDS[$i]}"

done

echo
echo "Created ${#DIRS[@]} calculations."

if [[ "$SUBMIT" == true ]]; then
    echo "All jobs submitted."
else
    echo "Jobs generated only."
fi

if [[ "$VERBOSE" == true ]]; then

    echo
    echo "Calculation mode: $CALC"
    echo "Prefix: $PREFIX"
    echo "max(ENMAX): $ENMAX eV"

    if (( ${#SETS[@]} )); then

        echo "INCAR overrides:"

        printf \
            '  %s\n' \
            "${SETS[@]}"

    fi
fi
