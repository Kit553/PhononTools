#!/usr/bin/env bash

# APGU - Automated POTCAR Generation Utility
#
# Creates a VASP POTCAR from the element symbols in a VASP5+ POSCAR.
#
# The user first selects a VASP potential library and is then prompted to choose between the available potential variants for each element.GW potentials are ignored.
#
# Existing POTCAR files are never modified during generation. The new POTCAR is assembled in a temporary file in the current directory and replaces the old POTCAR only after every required potential has been found, selected, appended, and validated successfully.
#
# A POTCAR.sources file is generated alongside the POTCAR to record the exact potential directories used for reproducibility.
#
# Configuration:
#   VASP_POTCAR_DIR can be set in the environment to override the default potential library location.
#
# Requirements:
#   - Bash 4+
#   - VASP5+ POSCAR format with element symbols on line 6
#   - Standard VASP potential directory structure
#
# Usage:
#   apgu
#
# Output:
#   POTCAR
#   POTCAR.sources

set -Eeuo pipefail

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    cat <<'EOF'
Usage: apgu.sh

Interactively select a VASP potential library and variants for the species in POSCAR.
Writes POTCAR and POTCAR.sources after validating the assembled potentials.
Existing output files are replaced after successful generation.

Environment:
  VASP_POTCAR_DIR  Library root (default: $HOME/VASP_POTCAR)
  POSCAR          Input structure (default: POSCAR)
  POTCAR          Output potential (default: POTCAR)
EOF
    exit 0
fi
if (($#)); then
    echo "Error: no positional arguments are accepted; use --help for configuration." >&2
    exit 2
fi

POSCAR=${POSCAR:-POSCAR}
POTCAR=${POTCAR:-POTCAR}

POTLOC=${VASP_POTCAR_DIR:-"$HOME/VASP_POTCAR"}

declare -A POT_PATHS=(
    [pbe64]="$POTLOC/PAW_PBE_64"
    [lda64]="$POTLOC/PAW_LDA_64"
    [pbe54]="$POTLOC/PAW_PBE_54"
    [lda54]="$POTLOC/PAW_LDA_54"
    [gga_uspp]="$POTLOC/USPP_GGA"
    [lda_uspp]="$POTLOC/USPP_LDA"
    [pw91]="$POTLOC/PAW_PW91"
)

declare -A POT_LABELS=(
    [pbe64]="PBE PAW 64"
    [lda64]="LDA PAW 64"
    [pbe54]="PBE PAW 54"
    [lda54]="LDA PAW 54"
    [gga_uspp]="GGA USPP (legacy)"
    [lda_uspp]="LDA USPP (legacy)"
    [pw91]="PW91 PAW (legacy)"
)

declare -A POT_DESCRIPTIONS=(
    [_h]="hard"
    [_s]="soft"
    [_sv]="s and p semicore"
    [_pv]="p semicore"
    [_d]="d semicore"
    [_2]="2 f electrons semicore"
    [_3]="3 f electrons semicore"
)

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}


potential_description() {
    local name=$1
    local suffix

    if [[ $name =~ (_sv|_pv|_[23hsd])$ ]]; then
        suffix=${BASH_REMATCH[1]}
        printf '%s' "${POT_DESCRIPTIONS[$suffix]}"
    else
        printf '%s' "normal valence shell"
    fi
}


choose_library() {
    local choice

    echo "APGU - Automated POTCAR Generation Utility"
    echo
    echo "Available potential libraries:"
    echo "  1) PBE PAW 64"
    echo "  2) LDA PAW 64"
    echo "  3) GGA USPP (legacy)"
    echo "  4) LDA USPP (legacy)"
    echo "  5) PW91 PAW (legacy)"
    echo "  6) PBE PAW 54"
    echo "  7) LDA PAW 54"
    echo

    while true; do
        read -rp "Select potential library: " choice ||
            die "No selection received."

        case "$choice" in
            1) potcase=pbe64 ;;
            2) potcase=lda64 ;;
            3) potcase=gga_uspp ;;
            4) potcase=lda_uspp ;;
            5) potcase=pw91 ;;
            6) potcase=pbe54 ;;
            7) potcase=lda54 ;;
            *)
                echo "Invalid input, try again."
                continue
                ;;
        esac

        break
    done

    echo
    echo "Selected: ${POT_LABELS[$potcase]}"
    echo
}

append_potential() {
    local element=$1
    local root=$2
    local target=$3
    local source_target=$4

    local choice
    local selected
    local base
    local i

    local -a options=()

    [[ -d $root ]] ||
        die "Potential directory does not exist: $root"

    mapfile -t options < <(
        find "$root" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            \( -name "$element" -o -name "${element}_*" \) \
            ! -name '*GW*' \
            -print |
            sort
    )

    if ((${#options[@]} == 0)); then
        die "No non-GW potential found for $element in $root"
    fi

    if ((${#options[@]} == 1)); then
        selected=${options[0]}

        printf '%s: %s [%s; only potential available]\n' \
            "$element" \
            "${selected##*/}" \
            "$(potential_description "${selected##*/}")"

    else
        echo "Available potentials for $element:"

        for i in "${!options[@]}"; do
            base=${options[$i]##*/}

            printf '  %d) %-12s %s\n' \
                "$((i + 1))" \
                "$base" \
                "($(potential_description "$base"))"
        done

        while true; do
            read -rp "Choose potential for $element: " choice ||
                die "No selection received."

            if [[ $choice =~ ^[0-9]+$ ]] &&
               ((choice >= 1 && choice <= ${#options[@]})); then

                selected=${options[choice-1]}
                break
            fi

            echo "Invalid input, try again."
        done

        printf '%s: %s [%s]\n' \
            "$element" \
            "${selected##*/}" \
            "$(potential_description "${selected##*/}")"
    fi

    [[ -s "$selected/POTCAR" ]] ||
        die "POTCAR is missing or empty: $selected/POTCAR"

    cat -- "$selected/POTCAR" >> "$target"
    printf '%s\n' "$selected" >> "$source_target"

    echo
}

[[ -f $POSCAR ]] ||
    die "No $POSCAR found in $PWD"

read -r -a species < <(
    awk 'NR == 6 {print; exit}' "$POSCAR"
)

((${#species[@]} > 0)) ||
    die "Could not read element symbols from line 6 of $POSCAR."

for element in "${species[@]}"; do
    [[ $element =~ ^[A-Z][a-z]?$ ]] ||
        die "POSCAR line 6 does not appear to contain element symbols: '$element'"
done

choose_library

potential_root=${POT_PATHS[$potcase]}

[[ -d $potential_root ]] ||
    die "Selected potential library does not exist: $potential_root"

if [[ -e $POTCAR ]]; then
    echo "Existing $POTCAR detected."
    echo "It will remain untouched until the replacement has been completely generated."
    echo
fi

tmp_potcar=$(mktemp "./.${POTCAR}.tmp.XXXXXX")
tmp_sources=$(mktemp "./.${POTCAR}.sources.tmp.XXXXXX")

cleanup() {
    rm -f -- "$tmp_potcar" "$tmp_sources"
}

trap cleanup EXIT INT TERM

for element in "${species[@]}"; do
    append_potential \
        "$element" \
        "$potential_root" \
        "$tmp_potcar" \
        "$tmp_sources"
done

[[ -s $tmp_potcar ]] ||
    die "Generated POTCAR is empty."

expected_potentials=${#species[@]}

actual_potentials=$(
    grep -c '^[[:space:]]*TITEL[[:space:]]*=' "$tmp_potcar" || true
)

if ((actual_potentials != expected_potentials)); then
    die "Generated POTCAR contains $actual_potentials TITEL entries; expected $expected_potentials."
fi

actual_sources=$(
    wc -l < "$tmp_sources"
)

if ((actual_sources != expected_potentials)); then
    die "Source record contains $actual_sources entries; expected $expected_potentials."
fi

mv -f -- "$tmp_potcar" "$POTCAR"
mv -f -- "$tmp_sources" "${POTCAR}.sources"

trap - EXIT INT TERM


echo "POTCAR successfully created."
echo

echo "Potential sources:"
cat "${POTCAR}.sources"

echo
echo "Potential information:"

awk '
    /^[[:space:]]*TITEL[[:space:]]*=/ ||
    /^[[:space:]]*VRHFIN[[:space:]]*=/ ||
    /^[[:space:]]*ZVAL[[:space:]]*=/ ||
    /^[[:space:]]*ENMAX[[:space:]]*=/ {
        print
    }
' "$POTCAR"

echo
echo "Finished."
