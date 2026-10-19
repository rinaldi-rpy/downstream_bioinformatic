# Shared functions for Bash entry points. Caller sets SCRIPT_DIR.
set -Eeuo pipefail
trap 'echo "ERROR: ${BASH_SOURCE[0]}:$LINENO: $BASH_COMMAND" >&2' ERR
set -a
source "$SCRIPT_DIR/config/parameters.sh"
set +a
BASE_DIR=${BASE_DIR:-$SCRIPT_DIR}
CONFIG=${CONFIG:-$SCRIPT_DIR/config/parameters.sh}
[[ "$CONFIG" == "$SCRIPT_DIR/config/parameters.sh" ]] || { set -a; source "$CONFIG"; set +a; }
AUTOSOME_MAP=${AUTOSOME_MAP:-$SCRIPT_DIR/config/autosomes.tsv}
need(){ command -v "$1" >/dev/null || { echo "Missing executable: $1" >&2; exit 1; }; }
require_plink19(){
    need "$PLINK"
    local version
    version=$("$PLINK" --version)
    [[ "$version" == *"PLINK v1.9"* ]] || { echo "PLINK 1.9 required; found: $version" >&2; exit 1; }
}
require_file(){ [[ -s "$1" ]] || { echo "Missing input: $1" >&2; exit 1; }; }
export AUTOSOME_MAP GROUP_COLUMN
need "$PYTHON"
"$PYTHON" "$SCRIPT_DIR/python/validate_parameters.py"
