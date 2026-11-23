# Shared helpers sourced by every step. Not meant to be run directly.
#
# All steps are run from the repository root:
#   qsub -v VAR=value,VAR2=value script/NN_step.pbs     (PBS Pro / OpenPBS / Torque)
#   bash script/NN_step.pbs -x value                       (interactive / no scheduler)
#
# Optional environment variables understood by every step:
#   MODULES   space-separated list of environment modules to load (e.g. "star/2.7.8a samtools/1.13")
#   CONDA_ENV name of a conda environment to activate

# Move to the submission directory (repository root) so that relative paths work.
cd "${PBS_O_WORKDIR:-$(pwd)}"

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*" >&2; }
die()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; exit 1; }

require_cmd() {
  local c
  for c in "$@"; do
    command -v "$c" >/dev/null 2>&1 || die "'$c' not found in PATH (set MODULES or CONDA_ENV, or install it)"
  done
}

require_file() { [ -s "$1" ] || die "file not found or empty: $1"; }
require_dir()  { [ -d "$1" ] || die "directory not found: $1"; }

# Load site-specific software, if requested.
setup_env() {
  if [ -n "${MODULES:-}" ]; then
    if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
      # shellcheck disable=SC1091
      source /etc/profile.d/modules.sh
    fi
    command -v module >/dev/null 2>&1 || die "MODULES is set but the 'module' command is not available"
    # shellcheck disable=SC2086
    module load $MODULES
  fi
  if [ -n "${CONDA_ENV:-}" ]; then
    command -v conda >/dev/null 2>&1 || die "CONDA_ENV is set but 'conda' is not available"
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
  fi
}

# Number of CPUs granted by the scheduler, falling back to a default.
ncpus() { echo "${NCPUS:-${PBS_NUM_PPN:-${1:-4}}}"; }

# Iterate over trimmed reads produced by step 02 and print one line per sample:
#   <bioproject> <sample> <R1> <R2-or-empty>
# Layout: <trimdir>/paired/<BIOPROJECT>/<SAMPLE>_1.trimmed.fastq.gz (+ _2)
#         <trimdir>/single/<BIOPROJECT>/<SAMPLE>.trimmed.fastq.gz
list_trimmed_samples() {
  local trimdir="$1" r1 r2 bp s
  if [ -d "$trimdir/paired" ]; then
    find "$trimdir/paired" -mindepth 2 -maxdepth 2 -type f -name '*_1.trimmed.fastq.gz' | sort | while read -r r1; do
      r2="${r1%_1.trimmed.fastq.gz}_2.trimmed.fastq.gz"
      [ -f "$r2" ] || { warn "mate not found for $r1, skipping"; continue; }
      bp=$(basename "$(dirname "$r1")")
      s=$(basename "$r1" _1.trimmed.fastq.gz)
      printf '%s\t%s\t%s\t%s\n' "$bp" "$s" "$r1" "$r2"
    done
  fi
  if [ -d "$trimdir/single" ]; then
    find "$trimdir/single" -mindepth 2 -maxdepth 2 -type f -name '*.trimmed.fastq.gz' | sort | while read -r r1; do
      bp=$(basename "$(dirname "$r1")")
      s=$(basename "$r1" .trimmed.fastq.gz)
      printf '%s\t%s\t%s\t%s\n' "$bp" "$s" "$r1" ""
    done
  fi
}
