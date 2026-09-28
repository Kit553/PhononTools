#!/usr/bin/env bash
set -euo pipefail

PREFIX=""
VERBOSE=false

RELAX_MODE=""
REF_DIR=""

usage() {
    cat <<'TXT'
Usage:
  ./ConvCheck.sh --prefix PREFIX [OPTIONS]

Default/static mode:

  ./ConvCheck.sh --prefix ENCUT_

Relaxation modes:

  ./ConvCheck.sh --prefix EDIFFG_ --relax diff

      Compare each directory:
          POSCAR -> CONTCAR

  ./ConvCheck.sh --prefix EDIFFG_ --relax abs

      Report absolute final CONTCAR parameters.

  ./ConvCheck.sh --prefix EDIFFG_ --relax ref

      Compare all CONTCARs against the automatically selected
      strictest convergence setting.

  ./ConvCheck.sh --prefix EDIFFG_ --relax ref DIRECTORY

      Compare all CONTCARs against an explicitly selected
      reference directory.

Options:
  --prefix PREFIX
  --relax [diff|abs|ref]
  -v, --verbose
  -h, --help

Notes:
  Slurm accounting is obtained from sacct using the JobID stored
  by ConvGen.sh.

  If Slurm accounting is unavailable, OUTCAR wall/CPU timing is
  used as a fallback where possible.
TXT
}

die() {
    echo "Error: $*" >&2
    exit 1
}

is_known_option() {
    case "$1" in
        --prefix|--relax|-v|--verbose|-h|--help)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

meta_get() {
    local dir="$1"
    local key="$2"
    local file="$dir/.convgen_meta"

    [[ -f "$file" ]] || return 0

    awk \
        -F= \
        -v k="$key" '
        $1 == k {
            sub(/^[^=]*=/, "")
            print
            exit
        }
        ' \
        "$file"
}

format_seconds() {
    awk \
        -v s="$1" '
        BEGIN {

            if (s == "" || s < 0) {
                print "NA"
                exit
            }

            s = int(s + 0.5)

            d = int(s / 86400)
            s %= 86400

            h = int(s / 3600)
            s %= 3600

            m = int(s / 60)
            sec = s % 60

            if (d > 0)
                printf "%d-%02d:%02d:%02d", d,h,m,sec
            else
                printf "%02d:%02d:%02d", h,m,sec
        }
        '
}

time_to_seconds() {
    awk \
        -v t="$1" '
        BEGIN {

            if (
                t == "" \
                || t == "Unknown" \
                || t == "NA"
            ) {
                print ""
                exit
            }

            d = 0

            if (index(t, "-") > 0) {

                split(t, a, "-")

                d = a[1] + 0
                t = a[2]
            }

            n = split(t, b, ":")

            if (n == 3)
                s = d*86400 \
                  + b[1]*3600 \
                  + b[2]*60 \
                  + b[3]

            else if (n == 2)
                s = d*86400 \
                  + b[1]*60 \
                  + b[2]

            else
                s = d*86400 \
                  + b[1]

            printf "%.6f", s
        }
        '
}

# ============================================================
# Slurm resource information
#
# Output:
#
# wall |
# totalcpu |
# cpueff |
# maxrss_mib |
# state |
# jobid |
# alloccpus |
# reqmem
# ============================================================

get_resources() {
    local dir="$1"
    local outcar="$2"

    local jobid=""
    local acct=""
    local parent=""

    local state=""
    local elapsed=""
    local totalcpu=""
    local alloc=""
    local reqmem=""

    local maxrss=""
    local cpu_sec=""
    local eff=""
    local wall=""

    # --------------------------------------------------------
    # Recover JobID written by ConvGen.sh
    # --------------------------------------------------------

    if [[ -f "$dir/.convgen_jobid" ]]; then

        jobid=$(
            head -n1 "$dir/.convgen_jobid" \
            | cut -d';' -f1
        )

    else

        jobid=$(
            meta_get "$dir" JOBID || true
        )

    fi

    # --------------------------------------------------------
    # Slurm accounting
    # --------------------------------------------------------

    if \
        [[ -n "$jobid" ]] \
        && command -v sacct >/dev/null 2>&1
    then

        acct=$(
            sacct \
                -j "$jobid" \
                -P \
                -n \
                --units=M \
                -o \
JobIDRaw,State,ElapsedRaw,TotalCPU,AllocCPUS,MaxRSS,ReqMem \
                2>/dev/null \
            || true
        )

        parent=$(
            printf '%s\n' "$acct" \
            | awk \
                -F'|' \
                -v id="$jobid" '
                $1 == id {
                    print
                    exit
                }
                '
        )

        if [[ -n "$parent" ]]; then

            IFS='|' read -r \
                _ \
                state \
                elapsed \
                totalcpu \
                alloc \
                _ \
                reqmem \
                <<< "$parent"

            wall=$(
                format_seconds \
                    "${elapsed:-}"
            )

            cpu_sec=$(
                time_to_seconds \
                    "${totalcpu:-}"
            )

            if \
                [[ -n "$cpu_sec" ]] \
                && [[ -n "${elapsed:-}" ]] \
                && [[ -n "${alloc:-}" ]] \
                && awk \
                    -v e="$elapsed" \
                    -v a="$alloc" '
                    BEGIN {
                        exit !(e > 0 && a > 0)
                    }
                    '
            then

                eff=$(
                    awk \
                        -v c="$cpu_sec" \
                        -v e="$elapsed" \
                        -v a="$alloc" '
                        BEGIN {
                            printf "%.1f%%", \
                                100*c/(e*a)
                        }
                        '
                )

            fi

            # ------------------------------------------------
            # MaxRSS is commonly populated on job-step rows, not the parent allocation row.
            #
            # --units=M normalizes the reported memory unit.
            # ------------------------------------------------

            maxrss=$(
                printf '%s\n' "$acct" \
                | awk \
                    -F'|' '
                    {
                        x = $6

                        gsub(
                            /[[:space:]]/,
                            "",
                            x
                        )

                        sub(
                            /[Mm]$/,
                            "",
                            x
                        )

                        if (
                            x ~ /^[0-9.]+$/ \
                            && x+0 > max
                        )
                            max = x+0
                    }

                    END {
                        if (max > 0)
                            printf "%.1f", max
                    }
                    '
            )

        fi
    fi

    # --------------------------------------------------------
    # Fallback: VASP OUTCAR timing
    # --------------------------------------------------------

    if [[ -z "$wall" || "$wall" == "NA" ]]; then

        if [[ -f "$outcar" ]]; then

            sec=$(
                awk '
                /Elapsed time \(sec\)/ {
                    x=$4
                }

                END {
                    print x
                }
                ' \
                "$outcar"
            )

            if [[ -n "$sec" ]]; then
                wall=$(format_seconds "$sec")
            fi
        fi
    fi

    if [[ -z "$totalcpu" && -f "$outcar" ]]; then

        csec=$(
            awk '
            /Total CPU time used \(sec\)/ {
                x=$6
            }

            END {
                print x
            }
            ' \
            "$outcar"
        )

        if [[ -n "$csec" ]]; then
            totalcpu=$(format_seconds "$csec")
        fi
    fi

    printf \
        '%s|%s|%s|%s|%s|%s|%s|%s\n' \
        "${wall:-NA}" \
        "${totalcpu:-NA}" \
        "${eff:-NA}" \
        "${maxrss:-NA}" \
        "${state:-NA}" \
        "${jobid:-NA}" \
        "${alloc:-NA}" \
        "${reqmem:-NA}"
}

# ============================================================
# OUTCAR values
#
# Output:
#
# ENCUT |
# ETOT |
# E/atom |
# MaxForce |
# MaxStress
# ============================================================

get_outcar_values() {
    local out="$1"

    local encut
    local etot
    local nions
    local epa
    local maxforce
    local maxstress

    if [[ ! -f "$out" ]]; then

        printf \
            'NA|NA|NA|NA|NA\n'

        return
    fi

    encut=$(
        awk '
        /ENCUT[[:space:]]*=/ {

            for (i=1; i<=NF; i++) {

                if ($i == "ENCUT") {

                    for (j=i+1; j<=NF; j++) {

                        if ($(j) ~ /^[0-9.]+$/) {
                            print $(j)
                            exit
                        }
                    }
                }
            }
        }
        ' \
        "$out"
    )

    etot=$(
        awk '
        /free  energy   TOTEN/ {
            e=$5
        }

        END {
            print e
        }
        ' \
        "$out"
    )

    nions=$(
        awk '
        /NIONS/ {

            for (i=1; i<=NF; i++) {

                if ($i == "NIONS") {

                    for (j=i+1; j<=NF; j++) {

                        if ($(j) ~ /^[0-9]+$/) {
                            n=$(j)
                            break
                        }
                    }
                }
            }
        }

        END {
            print n
        }
        ' \
        "$out"
    )

    if \
        [[ -n "$etot" ]] \
        && [[ -n "$nions" ]] \
        && [[ "$nions" != "0" ]]
    then

        epa=$(
            awk \
                -v e="$etot" \
                -v n="$nions" '
                BEGIN {
                    printf "%.8f", e/n
                }
                '
        )

    else

        epa="NA"

    fi

    # --------------------------------------------------------
    # Maximum norm of the atomic force vector in the final TOTAL-FORCE block.
    # --------------------------------------------------------

    maxforce=$(
        awk '
        /TOTAL-FORCE/ {
            inside=1
            max=0
            next
        }

        inside \
        && NF >= 6 \
        && $1 ~ /^[-+0-9.]/ {

            f = sqrt(
                $4*$4 \
                + $5*$5 \
                + $6*$6
            )

            if (f > max)
                max=f

            found=1
        }

        END {

            if (found)
                printf "%.6f", max
            else
                print "NA"
        }
        ' \
        "$out"
    )

    # --------------------------------------------------------
    # Largest absolute component of the final VASP stress tensor reported on the "in kB" line.
    # --------------------------------------------------------

    maxstress=$(
        awk '
        /in kB/ {

            max=0

            for (i=3; i<=8 && i<=NF; i++) {

                x=$i+0

                if (x < 0)
                    x=-x

                if (x > max)
                    max=x
            }

            last=max
        }

        END {

            if (last != "")
                printf "%.4f", last
            else
                print "NA"
        }
        ' \
        "$out"
    )

    printf \
        '%s|%s|%s|%s|%s\n' \
        "${encut:-NA}" \
        "${etot:-NA}" \
        "$epa" \
        "$maxforce" \
        "$maxstress"
}

# ============================================================
# Absolute lattice parameters from POSCAR/CONTCAR
#
# Output:
#
# a | b | c | alpha | beta | gamma | volume
# ============================================================

structure_abs() {

    python3 \
        - "$1" \
        <<'PY'
import sys
import math

path = sys.argv[1]

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def cross(a, b):
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]

def norm(v):
    return math.sqrt(dot(v, v))

def read_cell(path):

    raw = [
        line.rstrip()
        for line in open(path)
    ]

    scale = list(
        map(
            float,
            raw[1].split()
        )
    )

    vec = [
        list(
            map(
                float,
                raw[i].split()[:3]
            )
        )
        for i in range(2, 5)
    ]

    if len(scale) == 1:

        s = scale[0]

        if s == 0:
            raise ValueError(
                "zero POSCAR scale"
            )

        if s < 0:

            rv = abs(
                dot(
                    vec[0],
                    cross(
                        vec[1],
                        vec[2]
                    )
                )
            )

            s = (
                -s / rv
            ) ** (1/3)

        cell = [
            [s*x for x in v]
            for v in vec
        ]

    elif \
        len(scale) == 3 \
        and all(x > 0 for x in scale):

        cell = [
            [
                v[0]*scale[0],
                v[1]*scale[1],
                v[2]*scale[2],
            ]
            for v in vec
        ]

    else:

        raise ValueError(
            "unsupported POSCAR scale line"
        )

    return cell

cell = read_cell(path)

a, b, c = cell

la, lb, lc = map(
    norm,
    cell
)

def angle(u, v):

    x = (
        dot(u, v)
        /
        (
            norm(u)
            * norm(v)
        )
    )

    x = max(
        -1,
        min(1, x)
    )

    return math.degrees(
        math.acos(x)
    )

alpha = angle(b, c)
beta  = angle(a, c)
gamma = angle(a, b)

V = abs(
    dot(
        a,
        cross(b, c)
    )
)

print(
    f"{la:.8f}|"
    f"{lb:.8f}|"
    f"{lc:.8f}|"
    f"{alpha:.6f}|"
    f"{beta:.6f}|"
    f"{gamma:.6f}|"
    f"{V:.6f}"
)
PY
}

# ============================================================
# Structural difference
#
# Cell:
#   direct difference in a,b,c,angles,volume
#
# Atomic positions:
#   minimum-image wrapped fractional displacement converted
#   with the average of the two cells.
#
# This separates homogeneous cell deformation from internal atomic rearrangement.
#
# Output:
#
# da | db | dc |
# dalpha | dbeta | dgamma |
# dV_percent |
# max_atom_shift |
# median_atom_shift
# ============================================================

structure_diff() {

    python3 \
        - "$1" "$2" \
        <<'PY'
import sys
import math
import statistics

p0, p1 = sys.argv[1:]

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def cross(a, b):
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]

def norm(v):
    return math.sqrt(dot(v, v))

def parse(path):

    lines = [
        x.rstrip()
        for x in open(path)
    ]

    scale = list(
        map(
            float,
            lines[1].split()
        )
    )

    rawcell = [
        list(
            map(
                float,
                lines[i].split()[:3]
            )
        )
        for i in range(2, 5)
    ]

    # --------------------------------------------------------
    # Scaling
    # --------------------------------------------------------

    if len(scale) == 1:

        s = scale[0]

        if s == 0:
            raise ValueError(
                "zero POSCAR scale"
            )

        if s < 0:

            rv = abs(
                dot(
                    rawcell[0],
                    cross(
                        rawcell[1],
                        rawcell[2]
                    )
                )
            )

            s = (
                -s / rv
            ) ** (1/3)

        cell = [
            [s*x for x in v]
            for v in rawcell
        ]

        cart_scale = [
            s, s, s
        ]

    elif \
        len(scale) == 3 \
        and all(x > 0 for x in scale):

        cell = [
            [
                v[0]*scale[0],
                v[1]*scale[1],
                v[2]*scale[2],
            ]
            for v in rawcell
        ]

        cart_scale = scale

    else:

        raise ValueError(
            "unsupported POSCAR scale line"
        )

    # --------------------------------------------------------
    # Species/count line
    # --------------------------------------------------------

    idx = 5

    tok = lines[idx].split()

    def all_int(xs):

        try:
            [int(x) for x in xs]
            return True

        except ValueError:
            return False

    if all_int(tok):

        symbols = None

        counts = list(
            map(
                int,
                tok
            )
        )

        idx += 1

    else:

        symbols = tok

        counts = list(
            map(
                int,
                lines[idx+1].split()
            )
        )

        idx += 2

    while not lines[idx].strip():
        idx += 1

    if \
        lines[idx] \
        .strip() \
        .lower() \
        .startswith("s"):

        idx += 1

    mode = (
        lines[idx]
        .strip()
        .lower()
    )

    direct = mode.startswith("d")
    cart = mode.startswith(
        ("c", "k")
    )

    if not (direct or cart):

        raise ValueError(
            "unknown coordinate mode"
        )

    idx += 1

    n = sum(counts)

    pos = []

    for j in range(n):

        pos.append(
            list(
                map(
                    float,
                    lines[idx+j]
                    .split()[:3]
                )
            )
        )

    a, b, c = cell

    V = dot(
        a,
        cross(b, c)
    )

    def cart_to_frac(r):

        return [
            dot(
                r,
                cross(b, c)
            ) / V,

            dot(
                r,
                cross(c, a)
            ) / V,

            dot(
                r,
                cross(a, b)
            ) / V,
        ]

    if direct:

        frac = pos

    else:

        scaled = [
            [
                r[0]*cart_scale[0],
                r[1]*cart_scale[1],
                r[2]*cart_scale[2],
            ]
            for r in pos
        ]

        frac = [
            cart_to_frac(r)
            for r in scaled
        ]

    return \
        cell, \
        frac, \
        counts, \
        symbols

def lattice(cell):

    a, b, c = cell

    la, lb, lc = map(
        norm,
        cell
    )

    def angle(u, v):

        x = (
            dot(u, v)
            /
            (
                norm(u)
                * norm(v)
            )
        )

        x = max(
            -1,
            min(1, x)
        )

        return math.degrees(
            math.acos(x)
        )

    return (
        la,
        lb,
        lc,
        angle(b, c),
        angle(a, c),
        angle(a, b),
        abs(
            dot(
                a,
                cross(b, c)
            )
        ),
    )

c0, f0, n0, s0 = parse(p0)
c1, f1, n1, s1 = parse(p1)

if \
    n0 != n1 \
    or len(f0) != len(f1):

    raise ValueError(
        "atom counts differ"
    )

if \
    s0 is not None \
    and s1 is not None \
    and s0 != s1:

    raise ValueError(
        "species order differs"
    )

l0 = lattice(c0)
l1 = lattice(c1)

d = [
    l1[i] - l0[i]
    for i in range(6)
]

dv = (
    100
    * (
        l1[6] - l0[6]
    )
    / l0[6]
)

# ------------------------------------------------------------
# Average cell for internal-coordinate displacement
# ------------------------------------------------------------

avg = [
    [
        (
            c0[i][j]
            + c1[i][j]
        ) / 2
        for j in range(3)
    ]
    for i in range(3)
]

shifts = []

for x0, x1 in zip(f0, f1):

    df = [
        x1[i] - x0[i]
        for i in range(3)
    ]

    # Minimum-image convention
    df = [
        x - round(x)
        for x in df
    ]

    dr = [
        df[0]*avg[0][j]
        + df[1]*avg[1][j]
        + df[2]*avg[2][j]
        for j in range(3)
    ]

    shifts.append(
        norm(dr)
    )

print(
    "|".join(
        [
            f"{d[0]:.8f}",
            f"{d[1]:.8f}",
            f"{d[2]:.8f}",
            f"{d[3]:.6f}",
            f"{d[4]:.6f}",
            f"{d[5]:.6f}",
            f"{dv:.6f}",
            f"{max(shifts):.8f}",
            f"{statistics.median(shifts):.8f}",
        ]
    )
)
PY
}

# ============================================================
# Automatically select strictest convergence calculation
# ============================================================

auto_reference() {
    local axis=""
    local best_dir=""
    local best_metric=""

    local dir
    local a
    local scan
    local effective
    local metric
    local direction
    local better

    for dir in "${DIRS[@]}"; do

        a=$(
            meta_get \
                "${dir%/}" \
                AXIS \
            || true
        )

        scan=$(
            meta_get \
                "${dir%/}" \
                SCAN_VALUE \
            || true
        )

        effective=$(
            meta_get \
                "${dir%/}" \
                EFFECTIVE_VALUE \
            || true
        )

        if [[ -z "$a" || -z "$scan" ]]; then

            die "Automatic ref needs .convgen_meta in every directory; specify a reference directory explicitly"

        fi

        if [[ -z "$axis" ]]; then

            axis="$a"

        elif [[ "$a" != "$axis" ]]; then

            die "Automatic ref requires all directories to use the same convergence axis"

        fi

        case "$axis" in

            encut)

                metric="$effective"
                direction="max"
                ;;

            kpoints)

                metric="$scan"
                direction="max"
                ;;

            ediff)

                metric="$scan"
                direction="min"
                ;;

            ediffg)

                metric=$(
                    awk \
                        -v x="$scan" '
                        BEGIN {
                            if (x < 0)
                                x=-x

                            print x
                        }
                        '
                )

                direction="min"
                ;;

            sigma)

                metric="$scan"
                direction="min"
                ;;

            nbands)

                metric="$effective"
                direction="max"
                ;;

            *)

                die "Cannot auto-select strictest setting for axis '$axis'"
                ;;
        esac

        if [[ -z "$best_dir" ]]; then

            best_dir="${dir%/}"
            best_metric="$metric"

        else

            if [[ "$direction" == "max" ]]; then

                better=$(
                    awk \
                        -v x="$metric" \
                        -v y="$best_metric" '
                        BEGIN {
                            print (x > y) ? 1 : 0
                        }
                        '
                )

            else

                better=$(
                    awk \
                        -v x="$metric" \
                        -v y="$best_metric" '
                        BEGIN {
                            print (x < y) ? 1 : 0
                        }
                        '
                )

            fi

            if [[ "$better" == "1" ]]; then

                best_dir="${dir%/}"
                best_metric="$metric"

            fi
        fi
    done

    printf '%s\n' "$best_dir"
}

# ============================================================
# Arguments
# ============================================================

while (( $# )); do

    case "$1" in

        --prefix)

            (( $# >= 2 )) \
                || die "--prefix requires a value"

            PREFIX="$2"

            shift 2
            ;;

        --relax)

            shift

            RELAX_MODE="diff"

            if \
                (( $# )) \
                && [[ "$1" =~ ^(diff|abs|ref)$ ]]
            then

                RELAX_MODE="$1"

                shift

                if \
                    [[ "$RELAX_MODE" == "ref" ]] \
                    && (( $# )) \
                    && ! is_known_option "$1"
                then

                    REF_DIR="${1%/}"

                    shift
                fi
            fi
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

[[ -n "$PREFIX" ]] \
    || die "--prefix is required"

if [[ -n "$RELAX_MODE" ]]; then

    command -v python3 >/dev/null 2>&1 \
        || die "python3 is required for --relax modes"

fi

# ============================================================
# Find matching directories
# ============================================================

shopt -s nullglob

DIRS=(
    "${PREFIX}"*/
)

(( ${#DIRS[@]} > 0 )) \
    || die "No directories matching '${PREFIX}*'"

mapfile \
    -t DIRS \
    < <(
        printf '%s\n' \
            "${DIRS[@]}" \
        | sort -V
    )

# ============================================================
# Reference mode
# ============================================================

if [[ "$RELAX_MODE" == "ref" ]]; then

    if [[ -z "$REF_DIR" ]]; then
        REF_DIR=$(auto_reference)
    fi

    [[ -f "$REF_DIR/CONTCAR" ]] \
        || die "Reference CONTCAR not found: $REF_DIR/CONTCAR"

    echo "Reference: $REF_DIR"
    echo
fi

# ============================================================
# Default/static table
# ============================================================

if [[ -z "$RELAX_MODE" ]]; then

    printf \
        '%-28s %-11s %-14s %15s %15s %13s %14s %12s %12s %9s %16s\n' \
        "Directory" \
        "Scan" \
        "Effective" \
        "ETOT(eV)" \
        "E/Atom(eV)" \
        "MaxF(eV/A)" \
        "MaxStress(kB)" \
        "WallTime" \
        "TotalCPU" \
        "CPUEff" \
        "MaxRSS/task(MiB)"

    printf \
        '%*s\n' \
        172 \
        '' \
        | tr ' ' '-'

    for d in "${DIRS[@]}"; do

        dir="${d%/}"

        scan=$(
            meta_get \
                "$dir" \
                SCAN_VALUE \
            || true
        )

        effective=$(
            meta_get \
                "$dir" \
                EFFECTIVE_VALUE \
            || true
        )

        [[ -n "$scan" ]] \
            || scan="-"

        [[ -n "$effective" ]] \
            || effective="-"

        IFS='|' read -r \
            _ \
            etot \
            epa \
            maxf \
            maxstress \
            <<< "$(
                get_outcar_values \
                    "$dir/OUTCAR"
            )"

        IFS='|' read -r \
            wall \
            totalcpu \
            cpueff \
            maxrss \
            state \
            jobid \
            alloc \
            reqmem \
            <<< "$(
                get_resources \
                    "$dir" \
                    "$dir/OUTCAR"
            )"

        printf \
            '%-28s %-11s %-14s %15s %15s %13s %14s %12s %12s %9s %16s\n' \
            "$dir" \
            "$scan" \
            "$effective" \
            "$etot" \
            "$epa" \
            "$maxf" \
            "$maxstress" \
            "$wall" \
            "$totalcpu" \
            "$cpueff" \
            "$maxrss"

        if [[ "$VERBOSE" == true ]]; then

            printf \
                '  Slurm: JobID=%s State=%s AllocCPUS=%s ReqMem=%s\n' \
                "$jobid" \
                "$state" \
                "$alloc" \
                "$reqmem"

        fi

    done

    exit 0
fi

# ============================================================
# Relaxation mode: diff
# POSCAR -> CONTCAR within each directory
# ============================================================

if [[ "$RELAX_MODE" == "diff" ]]; then

    printf \
        '%-26s %9s %9s %9s %9s %9s %9s %9s %12s %12s %12s %12s %9s %16s\n' \
        "Directory" \
        "da(A)" \
        "db(A)" \
        "dc(A)" \
        "dAlpha" \
        "dBeta" \
        "dGamma" \
        "dV(%)" \
        "MaxShift(A)" \
        "MedShift(A)" \
        "WallTime" \
        "TotalCPU" \
        "CPUEff" \
        "MaxRSS/task(MiB)"

    printf \
        '%*s\n' \
        178 \
        '' \
        | tr ' ' '-'

    for d in "${DIRS[@]}"; do

        dir="${d%/}"

        if \
            [[ ! -f "$dir/POSCAR" ]] \
            || [[ ! -f "$dir/CONTCAR" ]]
        then

            $VERBOSE \
                && echo \
                    "Skipping $dir: POSCAR/CONTCAR missing" \
                    >&2

            continue
        fi

        IFS='|' read -r \
            da \
            db \
            dc \
            dalpha \
            dbeta \
            dgamma \
            dv \
            maxshift \
            medshift \
            <<< "$(
                structure_diff \
                    "$dir/POSCAR" \
                    "$dir/CONTCAR"
            )"

        IFS='|' read -r \
            wall \
            totalcpu \
            cpueff \
            maxrss \
            state \
            jobid \
            alloc \
            reqmem \
            <<< "$(
                get_resources \
                    "$dir" \
                    "$dir/OUTCAR"
            )"

        printf \
            '%-26s %9.5f %9.5f %9.5f %9.4f %9.4f %9.4f %9.4f %12.6f %12.6f %12s %12s %9s %16s\n' \
            "$dir" \
            "$da" \
            "$db" \
            "$dc" \
            "$dalpha" \
            "$dbeta" \
            "$dgamma" \
            "$dv" \
            "$maxshift" \
            "$medshift" \
            "$wall" \
            "$totalcpu" \
            "$cpueff" \
            "$maxrss"

        if [[ "$VERBOSE" == true ]]; then

            printf \
                '  Slurm: JobID=%s State=%s AllocCPUS=%s ReqMem=%s\n' \
                "$jobid" \
                "$state" \
                "$alloc" \
                "$reqmem"

        fi
    done

    exit 0
fi

# ============================================================
# Relaxation mode: abs
# Absolute final CONTCAR parameters
# ============================================================

if [[ "$RELAX_MODE" == "abs" ]]; then

    printf \
        '%-26s %10s %10s %10s %9s %9s %9s %12s %15s %13s %12s %12s %9s %16s\n' \
        "Directory" \
        "a(A)" \
        "b(A)" \
        "c(A)" \
        "Alpha" \
        "Beta" \
        "Gamma" \
        "Volume(A3)" \
        "E/Atom(eV)" \
        "MaxF(eV/A)" \
        "WallTime" \
        "TotalCPU" \
        "CPUEff" \
        "MaxRSS/task(MiB)"

    printf \
        '%*s\n' \
        181 \
        '' \
        | tr ' ' '-'

    for d in "${DIRS[@]}"; do

        dir="${d%/}"

        if [[ ! -f "$dir/CONTCAR" ]]; then

            $VERBOSE \
                && echo \
                    "Skipping $dir: CONTCAR missing" \
                    >&2

            continue
        fi

        IFS='|' read -r \
            a \
            b \
            c \
            alpha \
            beta \
            gamma \
            vol \
            <<< "$(
                structure_abs \
                    "$dir/CONTCAR"
            )"

        IFS='|' read -r \
            _ \
            _ \
            epa \
            maxf \
            _ \
            <<< "$(
                get_outcar_values \
                    "$dir/OUTCAR"
            )"

        IFS='|' read -r \
            wall \
            totalcpu \
            cpueff \
            maxrss \
            state \
            jobid \
            alloc \
            reqmem \
            <<< "$(
                get_resources \
                    "$dir" \
                    "$dir/OUTCAR"
            )"

        printf \
            '%-26s %10.5f %10.5f %10.5f %9.4f %9.4f %9.4f %12.4f %15s %13s %12s %12s %9s %16s\n' \
            "$dir" \
            "$a" \
            "$b" \
            "$c" \
            "$alpha" \
            "$beta" \
            "$gamma" \
            "$vol" \
            "$epa" \
            "$maxf" \
            "$wall" \
            "$totalcpu" \
            "$cpueff" \
            "$maxrss"

        if [[ "$VERBOSE" == true ]]; then

            printf \
                '  Slurm: JobID=%s State=%s AllocCPUS=%s ReqMem=%s\n' \
                "$jobid" \
                "$state" \
                "$alloc" \
                "$reqmem"

        fi
    done

    exit 0
fi

# ============================================================
# Relaxation mode: ref
# CONTCAR_i -> reference CONTCAR
# ============================================================

if [[ "$RELAX_MODE" == "ref" ]]; then

    printf \
        '%-26s %9s %9s %9s %9s %9s %9s %9s %12s %12s %12s %12s %9s %16s\n' \
        "Directory" \
        "da(A)" \
        "db(A)" \
        "dc(A)" \
        "dAlpha" \
        "dBeta" \
        "dGamma" \
        "dV(%)" \
        "MaxShift(A)" \
        "MedShift(A)" \
        "WallTime" \
        "TotalCPU" \
        "CPUEff" \
        "MaxRSS/task(MiB)"

    printf \
        '%*s\n' \
        178 \
        '' \
        | tr ' ' '-'

    for d in "${DIRS[@]}"; do

        dir="${d%/}"

        if [[ ! -f "$dir/CONTCAR" ]]; then

            $VERBOSE \
                && echo \
                    "Skipping $dir: CONTCAR missing" \
                    >&2

            continue
        fi

        IFS='|' read -r \
            da \
            db \
            dc \
            dalpha \
            dbeta \
            dgamma \
            dv \
            maxshift \
            medshift \
            <<< "$(
                structure_diff \
                    "$REF_DIR/CONTCAR" \
                    "$dir/CONTCAR"
            )"

        IFS='|' read -r \
            wall \
            totalcpu \
            cpueff \
            maxrss \
            state \
            jobid \
            alloc \
            reqmem \
            <<< "$(
                get_resources \
                    "$dir" \
                    "$dir/OUTCAR"
            )"

        printf \
            '%-26s %9.5f %9.5f %9.5f %9.4f %9.4f %9.4f %9.4f %12.6f %12.6f %12s %12s %9s %16s\n' \
            "$dir" \
            "$da" \
            "$db" \
            "$dc" \
            "$dalpha" \
            "$dbeta" \
            "$dgamma" \
            "$dv" \
            "$maxshift" \
            "$medshift" \
            "$wall" \
            "$totalcpu" \
            "$cpueff" \
            "$maxrss"

        if [[ "$VERBOSE" == true ]]; then

            printf \
                '  Slurm: JobID=%s State=%s AllocCPUS=%s ReqMem=%s\n' \
                "$jobid" \
                "$state" \
                "$alloc" \
                "$reqmem"

        fi
    done

    exit 0
fi
