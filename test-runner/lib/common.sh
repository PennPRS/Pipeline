# lib/common.sh — shared helpers for the PennPRS test runner.
# Sourced by run_all.sh and the other lib/*.sh files. Not executable on its own.

# -------- colors & logging --------------------------------------------------

if [[ -t 1 ]] && [[ -z "${NO_COLOR:-}" ]]; then
    C_RED=$'\e[31m'; C_GREEN=$'\e[32m'; C_YELLOW=$'\e[33m'
    C_BLUE=$'\e[34m'; C_BOLD=$'\e[1m'; C_RESET=$'\e[0m'
else
    C_RED= C_GREEN= C_YELLOW= C_BLUE= C_BOLD= C_RESET=
fi

_ts() { date '+%Y-%m-%d %H:%M:%S'; }

log_info()  { printf '%s %sINFO%s  %s\n'  "$(_ts)" "$C_BLUE"   "$C_RESET" "$*"; }
log_ok()    { printf '%s %sPASS%s  %s\n'  "$(_ts)" "$C_GREEN"  "$C_RESET" "$*"; }
log_warn()  { printf '%s %sWARN%s  %s\n'  "$(_ts)" "$C_YELLOW" "$C_RESET" "$*" >&2; }
log_error() { printf '%s %sFAIL%s  %s\n'  "$(_ts)" "$C_RED"    "$C_RESET" "$*" >&2; }

die() { log_error "$*"; exit 1; }

# -------- yaml parsing ------------------------------------------------------
# We rely on `yq` (either mikefarah/yq v4+ or the python-based kislyuk/yq).
# Both accept jq-style queries via `-r`; we detect which flavor at runtime.

yq_flavor() {
    # Caches the detected flavor in _YQ_FLAVOR. "go" = mikefarah, "py" = python.
    if [[ -n "${_YQ_FLAVOR:-}" ]]; then printf '%s' "$_YQ_FLAVOR"; return; fi
    command -v yq >/dev/null 2>&1 || die \
        "yq not found. Install with 'conda install -c conda-forge yq' or 'pip install yq'."
    if yq --version 2>&1 | grep -qi 'mikefarah'; then
        _YQ_FLAVOR=go
    else
        _YQ_FLAVOR=py
    fi
    printf '%s' "$_YQ_FLAVOR"
}

yq_q() {
    # yq_q <query> <file>  — evaluate a jq-style query, always raw output.
    local query="$1" file="$2"
    case "$(yq_flavor)" in
        go) yq eval "$query" "$file" ;;
        py) yq -r "$query" "$file" ;;
    esac
}

# -------- variable expansion ------------------------------------------------
# Expand ${VAR} references inside a string using already-exported variables.
# Pure bash — no eval, no awk. Iterates to fixed point so ${homedir} -> expansion
# containing ${PennPRS_path} resolves in one call.

expand_vars() {
    # expand_vars <string>  — prints the expanded result.
    local s="$1" prev="" name rest
    while [[ "$s" != "$prev" ]]; do
        prev="$s"
        local out=""
        while [[ "$s" == *'${'*'}'* ]]; do
            out+="${s%%'${'*}"            # prefix before first ${
            s="${s#*'${'}"                # strip through ${
            name="${s%%'}'*}"             # variable name
            rest="${s#*'}'}"              # remainder after }
            if [[ "$name" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
                out+="${!name-}"          # expand; empty if unset
            else
                out+="\${${name}}"        # not a var reference; leave literal
            fi
            s="$rest"
        done
        s="${out}${s}"
    done
    printf '%s' "$s"
}

# -------- runner-wide paths -------------------------------------------------

: "${RUNNER_DIR:?RUNNER_DIR must be set by caller}"
LOGS_DIR="${RUNNER_DIR}/logs"
RESULTS_DIR="${RUNNER_DIR}/results"
mkdir -p "$LOGS_DIR" "$RESULTS_DIR"
