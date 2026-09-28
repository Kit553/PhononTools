#!/usr/bin/env bash

# cat-tools phono3py sample
#
# Generate batched SLURM jobs for q-point resolved phono3py calculations.
#
# The wrapper centralizes grid-point selection, batching, SLURM settings and
# commonly used phono3py options. It supports both legacy phono3py 2/3 and
# phono3py 4 where practical. In particular:
#
#   - legacy Wigner transport uses --wigner
#   - phono3py >= 4 uses --tt wte and requires the phono3py-wte plugin
#   - --nac is passed only to phono3py < 4; phono3py >= 4 auto-detects BORN
#   - --nonac is passed only to phono3py >= 4
#
# Supported modes:
#
#   rta        RTA thermal-transport grid-point sampling (--br)
#   wigner     Wigner thermal transport
#   lbte       Direct LBTE transport (--lbte)
#   jdos       Joint / weighted joint DOS (--jdos)
#   spf        Spectral function (--spf)
#   ise        Imaginary self-energy (--ise)
#   rse        Real self-energy (--rse)
#   gruneisen  Mode-Gruneisen parameters on a mesh (--gruneisen)
#   stp        Triplet-count estimate (--stp)
#
# Grid points are read, in order of precedence, from:
#
#   1. explicit grid points supplied after all options
#   2. --gp-file FILE
#   3. ir_grid_points.yaml in the working directory
#
# Explicit and file-based grid points are validated and deduplicated while
# preserving their original order.
#
# Usage:
#
#   cat-tools phono3py sample MODE M1 M2 M3 N_GP [OPTIONS] [GP ...]
#
# Examples:
#
#   cat-tools phono3py sample rta 21 21 21 10 --isotope --bmfp 4
#   cat-tools phono3py sample wigner 21 21 21 10 --bmfp 0.5 --isotope
#   cat-tools phono3py sample jdos 21 21 21 10 --ts "300" --num-freq-points 600
#   cat-tools phono3py sample spf 21 21 21 5 --ts "100 200 300" --bi "1 2, 3 4"
#   cat-tools phono3py sample ise 21 21 21 10 --gp-file missing_gps.txt
#   cat-tools phono3py sample rse 21 21 21 5 --sigma "0.05 0.1"
#   cat-tools phono3py sample rta 21 21 21 10 --read-phonon --write-gamma --submit
#
# Use --help for the complete option list.

set -Eeuo pipefail
shopt -s nullglob


usage() {
    cat <<'EOF'
Usage:
  cat-tools phono3py sample MODE M1 M2 M3 N_GP [OPTIONS] [GP ...]

Modes:
  rta        RTA thermal transport (--br)
  wigner     Wigner transport; legacy --wigner or phono3py>=4 --tt wte
  lbte       Direct solution of the linearized BTE (--lbte)
  jdos       Joint / weighted joint density of states (--jdos)
  spf        Spectral function (--spf)
  ise        Imaginary self-energy (--ise)
  rse        Real self-energy (--rse)
  gruneisen  Mode-Gruneisen parameters (--gruneisen)
  stp        Number of q-point triplets (--stp)

Grid-point selection:
  --gp-file FILE              Read grid points from FILE
  GP ...                      Explicit grid-point indices
                              Default: read ir_grid_points.yaml

Phono3py sampling:
  --ts "T1 T2 ..."            Explicit temperatures
  --tmin T                    Minimum temperature
  --tmax T                    Maximum temperature
  --tstep T                   Temperature step
  --num-freq-points N         Number of frequency sampling points
  --freq-pitch VALUE          Frequency sampling interval
  --bi "..."                  Band indices
  --sigma "S1 S2 ..."         Gaussian broadening value(s)
  --sigma-cutoff VALUE        Gaussian cutoff in standard deviations
  --thm                       Explicitly request tetrahedron integration
  --q-direction "X Y Z"       q->0 direction, mainly for NAC calculations

Transport and scattering:
  --bmfp VALUE                Boundary mean free path in micrometres
  --isotope                   Use natural isotope scattering
  --mv "V1 V2 ..."            Explicit mass variances
  --nu                        Separate Normal and Umklapp scattering
  --transport-type TYPE       njc23, ibdb19, smm19, or wte
  --full-pp                   Calculate full ph-ph interaction strengths
  --ave-pp                    Use averaged ph-ph interaction
  --const-ave-pp VALUE        Constant averaged ph-ph interaction
  --scattering-event-class N  ISE scattering class 1 or 2

Force constants and phonons:
  --fc2 / --no-fc2            Read / do not explicitly read fc2.hdf5
  --fc3 / --no-fc3            Read / do not explicitly read fc3.hdf5
  --pa VALUE                  Primitive axes; default AUTO
  --nac                       Request NAC
  --nonac                     Explicitly disable NAC

Phono3py I/O:
  --write-gamma               Write grid-point gamma/kappa files
  --read-gamma                Read previously written gamma files
  --write-gamma-detail        Write detailed triplet-resolved gamma
  --write-phonon              Write phonon mesh HDF5
  --read-phonon               Read phonon mesh HDF5
  --write-pp                  Write ph-ph interaction HDF5
  --read-pp                   Read ph-ph interaction HDF5
  --hdf5-compression VALUE    gzip, lzf, None, or compression level
  --output-suffix TEXT        Pass phono3py output filename suffix with -o

Expert passthrough:
  --ph3-arg ARG               Append one raw phono3py argument.
                              Repeat for multiple arguments.
                              Example:
                                --ph3-arg=--foo --ph3-arg=value

SLURM / generation:
  --job-dir DIR               Generated-job directory
  --tag TAG                   Job and output tag
  --partition NAME            SLURM partition
  --time HH:MM:SS             SLURM wall time
  --mem MEMORY                Memory per job, e.g. 4G
  --cpus N                    CPUs per task
  --nodes N                   Number of nodes
  --ntasks N                  Number of tasks
  --account NAME              SLURM account
  --qos NAME                  SLURM QoS
  --module "MODULES..."        Module stack loaded in generated jobs
  --no-module                 Do not run module load
  --submit                    Submit jobs with sbatch after generation
  --overwrite                 Replace existing generated job scripts
  --dry-run                   Print planned jobs without writing/submitting
  -h, --help                  Show this help

Environment defaults:
  JOB_DIR
  JOB_TAG
  PHONO3PY_MODULES
  SLURM_PARTITION
  SLURM_TIME
  SLURM_MEM
  SLURM_CPUS
  SLURM_ACCOUNT
  SLURM_QOS

Notes:
  * phono3py >= 4 removed --wigner. Wigner mode uses --tt wte and therefore
    requires the external phono3py-wte plugin.
  * phono3py >= 4 enables NAC automatically when BORN/nac_params are present.
    The wrapper therefore treats --nac as a compatibility request rather than
    blindly passing the removed --nac option.
  * Options not explicitly represented by this wrapper can be forwarded with
    repeated --ph3-arg arguments.
EOF
}


die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}


warn() {
    printf 'Warning: %s\n' "$*" >&2
}


is_pos_int() {
    [[ ${1:-} =~ ^[1-9][0-9]*$ ]]
}


is_nonneg_int() {
    [[ ${1:-} =~ ^[0-9]+$ ]]
}


is_number() {
    [[ ${1:-} =~ ^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]
}


require_value() {
    local opt=$1
    local count=$2
    ((count >= 2)) || die "$opt requires a value."
}


shell_join() {
    local out=""
    local arg
    local quoted

    for arg in "$@"; do
        printf -v quoted '%q' "$arg"
        out+="${out:+ }${quoted}"
    done

    printf '%s' "$out"
}


sanitize_tag() {
    printf '%s' "$1" |
        sed -E 's/[^A-Za-z0-9._-]+/_/g; s/^_+//; s/_+$//'
}


if (($# == 0)); then
    usage
    exit 1
fi

case "$1" in
    -h|--help)
        usage
        exit 0
        ;;
esac

(($# >= 5)) || {
    usage >&2
    exit 2
}

mode=$1
mesh1=$2
mesh2=$3
mesh3=$4
n_gp=$5
shift 5

case "$mode" in
    rta|wigner|lbte|jdos|spf|ise|rse|gruneisen|stp)
        ;;
    *)
        die "unknown mode '$mode'."
        ;;
esac

for mesh in "$mesh1" "$mesh2" "$mesh3"; do
    is_pos_int "$mesh" || die "mesh dimensions must be positive integers."
done

is_pos_int "$n_gp" || die "N_GP must be a positive integer."


# Mode-specific defaults reflect the historical workflow but remain overridable.

case "$mode" in
    rta|wigner)
        default_ts="1 2 3 4 5 7 10 12 15 20 25 30 40 50 70 100 125 150 175 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000"
        default_mem="1G"
        ;;
    lbte)
        default_ts="300"
        default_mem="8G"
        ;;
    jdos)
        default_ts="300"
        default_mem="4G"
        ;;
    spf)
        default_ts="1 10 20 100 200 500"
        default_mem="1G"
        ;;
    ise|rse)
        default_ts="300"
        default_mem="4G"
        ;;
    gruneisen|stp)
        default_ts=""
        default_mem="1G"
        ;;
esac

temperatures=$default_ts
temperature_style=explicit

tmin=""
tmax=""
tstep=""

num_freq_points=600
num_freq_points_set=false
freq_pitch=""

band_indices=""
sigmas=""
sigma_cutoff=""
tetrahedron=false
q_direction=""

bmfp=""
isotope=false
mass_variances=""
normal_umklapp=false
transport_type=""
full_pp=false
ave_pp=false
const_ave_pp=""
scattering_event_class=""

read_fc2=true
read_fc3=true
primitive_axes=AUTO
nac_mode=default

write_gamma=false
read_gamma=false
write_gamma_detail=false
write_phonon=false
read_phonon=false
write_pp=false
read_pp=false
hdf5_compression=""
output_suffix=""

job_dir=${JOB_DIR:-generated_jobs}
job_tag=${JOB_TAG:-}
partition=${SLURM_PARTITION:-long}
walltime=${SLURM_TIME:-48:00:00}
memory=${SLURM_MEM:-$default_mem}
cpus=${SLURM_CPUS:-1}
nodes=1
ntasks=1
account=${SLURM_ACCOUNT:-}
qos=${SLURM_QOS:-}

module_stack=${PHONO3PY_MODULES:-"palma/2023a foss/2023a phono3py/2.7.0"}
load_modules=true

submit=false
overwrite=false
dry_run=false

gp_file=""
explicit_gps=()
extra_ph3_args=()


while (($#)); do
    case "$1" in
        --gp-file)
            require_value "$1" "$#"
            [[ -z "$gp_file" ]] || die "--gp-file was specified more than once."
            gp_file=$2
            shift 2
            ;;

        --ts)
            require_value "$1" "$#"
            temperatures=$2
            temperature_style=explicit
            shift 2
            ;;

        --tmin)
            require_value "$1" "$#"
            tmin=$2
            temperature_style=range
            shift 2
            ;;

        --tmax)
            require_value "$1" "$#"
            tmax=$2
            temperature_style=range
            shift 2
            ;;

        --tstep)
            require_value "$1" "$#"
            tstep=$2
            temperature_style=range
            shift 2
            ;;

        --num-freq-points)
            require_value "$1" "$#"
            num_freq_points=$2
            num_freq_points_set=true
            shift 2
            ;;

        --freq-pitch)
            require_value "$1" "$#"
            freq_pitch=$2
            shift 2
            ;;

        --bi)
            require_value "$1" "$#"
            band_indices=$2
            shift 2
            ;;

        --sigma)
            require_value "$1" "$#"
            sigmas=$2
            shift 2
            ;;

        --sigma-cutoff)
            require_value "$1" "$#"
            sigma_cutoff=$2
            shift 2
            ;;

        --thm)
            tetrahedron=true
            shift
            ;;

        --q-direction)
            require_value "$1" "$#"
            q_direction=$2
            shift 2
            ;;

        --bmfp|--boundary-mfp)
            require_value "$1" "$#"
            bmfp=$2
            shift 2
            ;;

        --isotope)
            isotope=true
            shift
            ;;

        --mv|--mass-variances)
            require_value "$1" "$#"
            mass_variances=$2
            shift 2
            ;;

        --nu)
            normal_umklapp=true
            shift
            ;;

        --transport-type|--tt)
            require_value "$1" "$#"
            transport_type=${2,,}
            shift 2
            ;;

        --full-pp)
            full_pp=true
            shift
            ;;

        --ave-pp)
            ave_pp=true
            shift
            ;;

        --const-ave-pp)
            require_value "$1" "$#"
            const_ave_pp=$2
            shift 2
            ;;

        --scattering-event-class)
            require_value "$1" "$#"
            scattering_event_class=$2
            shift 2
            ;;

        --fc2)
            read_fc2=true
            shift
            ;;

        --no-fc2)
            read_fc2=false
            shift
            ;;

        --fc3)
            read_fc3=true
            shift
            ;;

        --no-fc3)
            read_fc3=false
            shift
            ;;

        --pa)
            require_value "$1" "$#"
            primitive_axes=$2
            shift 2
            ;;

        --nac)
            [[ "$nac_mode" != nonac ]] || die "--nac and --nonac are mutually exclusive."
            nac_mode=nac
            shift
            ;;

        --nonac)
            [[ "$nac_mode" != nac ]] || die "--nac and --nonac are mutually exclusive."
            nac_mode=nonac
            shift
            ;;

        --write-gamma)
            write_gamma=true
            shift
            ;;

        --read-gamma)
            read_gamma=true
            shift
            ;;

        --write-gamma-detail)
            write_gamma_detail=true
            shift
            ;;

        --write-phonon)
            write_phonon=true
            shift
            ;;

        --read-phonon)
            read_phonon=true
            shift
            ;;

        --write-pp)
            write_pp=true
            shift
            ;;

        --read-pp)
            read_pp=true
            shift
            ;;

        --hdf5-compression)
            require_value "$1" "$#"
            hdf5_compression=$2
            shift 2
            ;;

        --output-suffix)
            require_value "$1" "$#"
            output_suffix=$2
            shift 2
            ;;

        --ph3-arg)
            require_value "$1" "$#"
            extra_ph3_args+=("$2")
            shift 2
            ;;

        --ph3-arg=*)
            extra_ph3_args+=("${1#*=}")
            shift
            ;;

        --job-dir)
            require_value "$1" "$#"
            job_dir=$2
            shift 2
            ;;

        --tag)
            require_value "$1" "$#"
            job_tag=$2
            shift 2
            ;;

        --partition)
            require_value "$1" "$#"
            partition=$2
            shift 2
            ;;

        --time)
            require_value "$1" "$#"
            walltime=$2
            shift 2
            ;;

        --mem)
            require_value "$1" "$#"
            memory=$2
            shift 2
            ;;

        --cpus)
            require_value "$1" "$#"
            cpus=$2
            shift 2
            ;;

        --nodes)
            require_value "$1" "$#"
            nodes=$2
            shift 2
            ;;

        --ntasks)
            require_value "$1" "$#"
            ntasks=$2
            shift 2
            ;;

        --account)
            require_value "$1" "$#"
            account=$2
            shift 2
            ;;

        --qos)
            require_value "$1" "$#"
            qos=$2
            shift 2
            ;;

        --module)
            require_value "$1" "$#"
            module_stack=$2
            load_modules=true
            shift 2
            ;;

        --no-module)
            load_modules=false
            shift
            ;;

        --submit)
            submit=true
            shift
            ;;

        --overwrite)
            overwrite=true
            shift
            ;;

        --dry-run)
            dry_run=true
            shift
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        --)
            shift
            explicit_gps+=("$@")
            break
            ;;

        -*)
            die "unknown option '$1'. Use --ph3-arg to forward unsupported phono3py options."
            ;;

        *)
            explicit_gps+=("$1")
            shift
            ;;
    esac
done

if [[ -n "$gp_file" ]] && ((${#explicit_gps[@]} > 0)); then
    die "use either --gp-file or explicit grid points, not both."
fi

if [[ "$temperature_style" == range ]]; then
    [[ -n "$tmin" && -n "$tmax" && -n "$tstep" ]] ||
        die "--tmin, --tmax and --tstep must be supplied together."

    [[ -z "$temperatures" || "$temperatures" == "$default_ts" ]] || true

    for value in "$tmin" "$tmax" "$tstep"; do
        is_number "$value" || die "invalid temperature range value '$value'."
    done
fi

if [[ -n "$freq_pitch" && "$num_freq_points_set" == true ]]; then
    die "--freq-pitch and --num-freq-points are mutually exclusive."
fi

if [[ -n "$freq_pitch" ]]; then
    is_number "$freq_pitch" || die "--freq-pitch must be numeric."
fi

is_pos_int "$num_freq_points" || die "--num-freq-points must be a positive integer."

if [[ -n "$bmfp" ]]; then
    is_number "$bmfp" || die "--bmfp must be numeric."
fi

if [[ -n "$sigma_cutoff" ]]; then
    is_number "$sigma_cutoff" || die "--sigma-cutoff must be numeric."
fi

if [[ -n "$const_ave_pp" ]]; then
    is_number "$const_ave_pp" || die "--const-ave-pp must be numeric."
fi

if [[ -n "$transport_type" ]]; then
    case "$transport_type" in
        njc23|ibdb19|smm19|wte) ;;
        *) die "unsupported transport type '$transport_type'." ;;
    esac
fi

if [[ -n "$scattering_event_class" ]]; then
    [[ "$mode" == ise ]] ||
        die "--scattering-event-class is only valid in ise mode."

    [[ "$scattering_event_class" == 1 || "$scattering_event_class" == 2 ]] ||
        die "--scattering-event-class must be 1 or 2."
fi

if [[ "$mode" == rse && -z "$sigmas" ]]; then
    die "rse mode requires --sigma because real self-energy uses smearing."
fi

if [[ "$mode" != rta && "$mode" != wigner && "$mode" != lbte ]]; then
    if [[ -n "$bmfp" || "$isotope" == true || -n "$mass_variances" ||
          "$normal_umklapp" == true || -n "$transport_type" ||
          "$ave_pp" == true || -n "$const_ave_pp" ]]; then
        warn "transport/scattering options were supplied to mode '$mode'; phono3py may reject or ignore them."
    fi
fi

if [[ "$mode" == wigner && -n "$transport_type" && "$transport_type" != wte ]]; then
    die "wigner mode fixes the modern transport type to 'wte'; do not combine it with --transport-type $transport_type."
fi

if [[ "$write_gamma" == true && "$read_gamma" == true ]]; then
    warn "--write-gamma and --read-gamma were both requested."
fi

if [[ "$write_pp" == true && "$read_pp" == true ]]; then
    warn "--write-pp and --read-pp were both requested."
fi

for value_name in cpus nodes ntasks; do
    value=${!value_name}
    is_pos_int "$value" || die "--$value_name must be a positive integer."
done

gps=()

if [[ -n "$gp_file" ]]; then
    [[ -f "$gp_file" ]] || die "grid-point file not found: $gp_file"

    mapfile -t gps < <(
        awk '
        {
            gsub(/[,;]/, " ")
            for (i = 1; i <= NF; i++) {
                print $i
            }
        }
        ' "$gp_file"
    )

elif ((${#explicit_gps[@]} > 0)); then
    gps=("${explicit_gps[@]}")

else
    [[ -f ir_grid_points.yaml ]] ||
        die "ir_grid_points.yaml not found. Generate it with phono3py --mesh ... --wgp or supply grid points explicitly."

    mapfile -t gps < <(
        awk '/grid_point:/ {print $NF}' ir_grid_points.yaml
    )
fi

((${#gps[@]} > 0)) || die "no grid points were found."

declare -A seen
clean_gps=()

for gp in "${gps[@]}"; do
    is_nonneg_int "$gp" || die "invalid grid point '$gp'."

    if [[ -z "${seen[$gp]+exists}" ]]; then
        clean_gps+=("$gp")
        seen[$gp]=1
    fi
done

gps=("${clean_gps[@]}")

common_args=(
    "--pa=$primitive_axes"
    "--mesh=$mesh1 $mesh2 $mesh3"
)

$read_fc2 && common_args+=(--fc2)
$read_fc3 && common_args+=(--fc3)

case "$mode" in
    rta)
        mode_args=(--br)
        ;;
    wigner)
        mode_args=(--br)
        ;;
    lbte)
        mode_args=(--lbte)
        ;;
    jdos)
        mode_args=(--jdos)
        ;;
    spf)
        mode_args=(--spf)
        ;;
    ise)
        mode_args=(--ise)
        ;;
    rse)
        mode_args=(--rse)
        ;;
    gruneisen)
        mode_args=(--gruneisen)
        ;;
    stp)
        mode_args=(--stp)
        ;;
esac

property_args=()

if [[ "$temperature_style" == explicit && -n "$temperatures" &&
      "$mode" != gruneisen && "$mode" != stp ]]; then
    property_args+=("--ts=$temperatures")
elif [[ "$temperature_style" == range ]]; then
    property_args+=("--tmin=$tmin" "--tmax=$tmax" "--tstep=$tstep")
fi

case "$mode" in
    jdos|spf|ise|rse)
        if [[ -n "$freq_pitch" ]]; then
            property_args+=("--freq-pitch=$freq_pitch")
        else
            property_args+=("--num-freq-points=$num_freq_points")
        fi
        ;;
esac

[[ -n "$band_indices" ]] && property_args+=("--bi=$band_indices")
[[ -n "$sigmas" ]] && property_args+=("--sigma=$sigmas")
[[ -n "$sigma_cutoff" ]] && property_args+=("--sigma-cutoff=$sigma_cutoff")
$tetrahedron && property_args+=(--thm)
[[ -n "$q_direction" ]] && property_args+=("--q-direction=$q_direction")

[[ -n "$bmfp" ]] && property_args+=("--bmfp=$bmfp")
$isotope && property_args+=(--isotope)
[[ -n "$mass_variances" ]] && property_args+=("--mv=$mass_variances")
$normal_umklapp && property_args+=(--nu)
[[ -n "$transport_type" ]] && property_args+=("--tt=$transport_type")
$full_pp && property_args+=(--full-pp)
$ave_pp && property_args+=(--ave-pp)
[[ -n "$const_ave_pp" ]] && property_args+=("--const-ave-pp=$const_ave_pp")
[[ -n "$scattering_event_class" ]] &&
    property_args+=("--scattering-event-class=$scattering_event_class")

$write_gamma && property_args+=(--write-gamma)
$read_gamma && property_args+=(--read-gamma)
$write_gamma_detail && property_args+=(--write-gamma-detail)
$write_phonon && property_args+=(--write-phonon)
$read_phonon && property_args+=(--read-phonon)
$write_pp && property_args+=(--write-pp)
$read_pp && property_args+=(--read-pp)
[[ -n "$hdf5_compression" ]] &&
    property_args+=("--hdf5-compression=$hdf5_compression")
[[ -n "$output_suffix" ]] &&
    property_args+=("-o" "$output_suffix")

property_args+=("${extra_ph3_args[@]}")

if [[ -z "$job_tag" ]]; then
    job_tag=$mode

    if [[ "$mode" == wigner && -n "$bmfp" ]]; then
        job_tag+="_bmfp${bmfp}"
    fi
fi

job_tag=$(sanitize_tag "$job_tag")
[[ -n "$job_tag" ]] || die "job tag became empty after sanitization."

workdir=$PWD
total_gps=${#gps[@]}
job_nr=0
generated_scripts=()

if [[ "$dry_run" == false ]]; then
    mkdir -p -- "$job_dir"
fi


for ((start=0; start<total_gps; start+=n_gp)); do
    batch=("${gps[@]:start:n_gp}")

    first_gp=${batch[0]}
    last_index=$((${#batch[@]} - 1))
    last_gp=${batch[$last_index]}
    gp_for_job="${batch[*]}"

    script_name="${job_dir}/job_${job_tag}_${first_gp}-${last_gp}.sh"
    output_name="sampling_${job_tag}_${first_gp}-${last_gp}.out"

    batch_args=(
        "${common_args[@]}"
        "${mode_args[@]}"
        "${property_args[@]}"
        "--gp=$gp_for_job"
    )

    if [[ "$dry_run" == true ]]; then
        printf 'Job %-4d GP %-18s -> %s\n' \
            "$job_nr" "$first_gp-$last_gp" "$script_name"
        printf '  phono3py %s\n' "$(shell_join "${batch_args[@]}")"
        ((job_nr += 1))
        continue
    fi

    if [[ -e "$script_name" && "$overwrite" == false ]]; then
        die "job script already exists: $script_name (use --overwrite to replace it)."
    fi

    tmp_script=$(mktemp "${job_dir}/.qmesh_job.XXXXXX")

    {
        echo '#!/usr/bin/env bash'
        echo
        printf '#SBATCH --job-name=%s\n' "${job_tag}_${first_gp}-${last_gp}"
        printf '#SBATCH --nodes=%s\n' "$nodes"
        printf '#SBATCH --ntasks=%s\n' "$ntasks"
        printf '#SBATCH --cpus-per-task=%s\n' "$cpus"
        printf '#SBATCH --partition=%s\n' "$partition"
        printf '#SBATCH --time=%s\n' "$walltime"
        printf '#SBATCH --mem=%s\n' "$memory"
        echo '#SBATCH --mail-type=NONE'
        [[ -n "$account" ]] && printf '#SBATCH --account=%s\n' "$account"
        [[ -n "$qos" ]] && printf '#SBATCH --qos=%s\n' "$qos"
        echo
        echo 'set -Eeuo pipefail'
        echo
        if [[ "$load_modules" == true ]]; then
            printf 'module load %s\n' "$module_stack"
            echo
        fi

        printf 'cd -- %q\n' "$workdir"
        echo
        cat <<'JOBVERSION'
command -v phono3py >/dev/null 2>&1 || {
    echo "Error: phono3py is not available in PATH." >&2
    exit 127
}

ph3_version=$(phono3py --version 2>&1 | awk 'match($0, /[0-9]+([.][0-9]+)+/) {print substr($0, RSTART, RLENGTH); exit}')
ph3_major=${ph3_version%%.*}

if [[ -z "$ph3_version" || ! "$ph3_major" =~ ^[0-9]+$ ]]; then
    echo "Error: unable to determine phono3py version." >&2
    exit 1
fi

echo "phono3py version: $ph3_version"

version_args=()
JOBVERSION

        if [[ "$mode" == wigner ]]; then
            cat <<'JOBWIGNER'
if ((ph3_major >= 4)); then
    version_args+=(--tt=wte)
    echo "Wigner backend: phono3py-wte plugin (--tt wte)"
else
    version_args+=(--wigner)
    echo "Wigner backend: legacy bundled --wigner"
fi
JOBWIGNER
        fi

        case "$nac_mode" in
            nac)
                cat <<'JOBNAC'
if ((ph3_major >= 4)); then
    if [[ ! -f BORN ]]; then
        echo "Warning: --nac was requested, but no BORN file is visible. phono3py >=4 enables NAC automatically only when NAC data are available." >&2
    fi
    echo "NAC: automatic phono3py >=4 handling"
else
    version_args+=(--nac)
    echo "NAC: legacy --nac"
fi
JOBNAC
                ;;
            nonac)
                cat <<'JOBNONAC'
if ((ph3_major >= 4)); then
    version_args+=(--nonac)
    echo "NAC: explicitly disabled"
else
    echo "NAC: disabled by legacy default"
fi
JOBNONAC
                ;;
        esac

        printf '\nargs=(\n'
        for arg in "${batch_args[@]}"; do
            printf '    %q\n' "$arg"
        done
        echo ')'
        echo
        printf 'printf "Grid points: %%s\\n" %q\n' "$gp_for_job"
        echo 'printf "Command: phono3py"'
        echo 'printf " %q" "${args[@]}" "${version_args[@]}"'
        echo 'printf "\n"'
        echo
        printf 'phono3py "${args[@]}" "${version_args[@]}" > %q 2>&1\n' "$output_name"
    } > "$tmp_script"

    chmod +x "$tmp_script"
    mv -f -- "$tmp_script" "$script_name"

    generated_scripts+=("$script_name")

    printf 'Generated %-55s GP: %s\n' "$script_name" "$gp_for_job"

    ((job_nr += 1))
done


echo
echo "Mode:          $mode"
echo "Mesh:          ${mesh1} x ${mesh2} x ${mesh3}"
echo "Grid points:   $total_gps"
echo "Jobs:          $job_nr"
echo "GP per job:    $n_gp"
echo "Job tag:       $job_tag"

if [[ "$temperature_style" == explicit && -n "$temperatures" ]]; then
    echo "Temperatures:  $temperatures"
elif [[ "$temperature_style" == range ]]; then
    echo "Temperatures:  $tmin:$tstep:$tmax"
fi

[[ -n "$bmfp" ]] && echo "BMFP [um]:     $bmfp"
[[ -n "$sigmas" ]] && echo "Sigma:         $sigmas"

case "$mode" in
    jdos|spf|ise|rse)
        if [[ -n "$freq_pitch" ]]; then
            echo "Freq pitch:    $freq_pitch"
        else
            echo "Freq points:   $num_freq_points"
        fi
        ;;
esac

if [[ "$dry_run" == true ]]; then
    echo "Dry run:       no files written."
    exit 0
fi

echo "Job scripts:   $job_dir"

if [[ "$submit" == true ]]; then
    command -v sbatch >/dev/null 2>&1 ||
        die "--submit requested, but sbatch is not available."

    echo
    echo "Submitting jobs:"

    for script_name in "${generated_scripts[@]}"; do
        printf '  %s -> ' "$script_name"
        (
            cd -- "$(dirname "$script_name")"
            sbatch "$(basename "$script_name")"
        )
    done
fi
