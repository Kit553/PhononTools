#!/usr/bin/env bash
set -euo pipefail

POSCAR="${1:-POSCAR}"
DENSITY="${2:-50}"

[[ -f "$POSCAR" ]] || { echo "Error: $POSCAR not found"; exit 1; }

read N1 N2 N3 < <(
awk -v D="$DENSITY" '
NR==2 { s=$1 }

NR==3 { ax=$1*s; ay=$2*s; az=$3*s }
NR==4 { bx=$1*s; by=$2*s; bz=$3*s }
NR==5 { cx=$1*s; cy=$2*s; cz=$3*s }

END {
    # Cell volume
    V = ax*(by*cz-bz*cy) \
      - ay*(bx*cz-bz*cx) \
      + az*(bx*cy-by*cx)

    # Reciprocal lattice vectors (without 2*pi)
    r1x=(by*cz-bz*cy)/V
    r1y=(bz*cx-bx*cz)/V
    r1z=(bx*cy-by*cx)/V

    r2x=(cy*az-cz*ay)/V
    r2y=(cz*ax-cx*az)/V
    r2z=(cx*ay-cy*ax)/V

    r3x=(ay*bz-az*by)/V
    r3y=(az*bx-ax*bz)/V
    r3z=(ax*by-ay*bx)/V

    b1=sqrt(r1x*r1x+r1y*r1y+r1z*r1z)
    b2=sqrt(r2x*r2x+r2y*r2y+r2z*r2z)
    b3=sqrt(r3x*r3x+r3y*r3y+r3z*r3z)

    n1=int(D*b1+0.5)
    n2=int(D*b2+0.5)
    n3=int(D*b3+0.5)

    if(n1<1) n1=1
    if(n2<1) n2=1
    if(n3<1) n3=1

    print n1, n2, n3
}' "$POSCAR"
)

cat > KPOINTS <<EOF
Automatic mesh, density=${DENSITY}
0
Gamma
$N1 $N2 $N3
0 0 0
EOF

echo "Generated KPOINTS: ${N1}x${N2}x${N3} (density=${DENSITY})"
