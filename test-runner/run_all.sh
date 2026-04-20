#!/usr/bin/env bash
# ============================================================================
# PennPRS test runner — submits and verifies examples 5.1 – 5.7 on SLURM.
#
# Usage:
#   ./run_all.sh                       # run every test sequentially
#   ./run_all.sh --only 5.1a,5.7       # run only these test ids
#   ./run_all.sh --only 5.1,5.2        # run all the test whose section starts with 5.1 / 5.2
#   ./run_all.sh --parallel 4          # submit up to 4 jobs concurrently
#   ./run_all.sh --dry-run             # print the sbatch command for each test, don't submit
#   ./run_all.sh --verify-only         # skip submission; re-run verification on existing output
#   ./run_all.sh --list                # list all test ids and exit
#
# Environment overrides:
#   CONFIG=/abs/path/tests.yaml        # use a different tests.yaml
#   NO_COLOR=1                         # disable ANSI colors
#
# Exit codes:
#   0   every selected test passed submission + verification
#   1   at least one test failed
#   2   invalid invocation
# ============================================================================

set -u -o pipefail

RUNNER_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RUNNER_DIR
CONFIG="${CONFIG:-${RUNNER_DIR}/tests.yaml}"

# shellcheck source=lib/common.sh
source "${RUNNER_DIR}/lib/common.sh"
# shellcheck source=lib/parse.sh
source "${RUNNER_DIR}/lib/parse.sh"
# shellcheck source=lib/submit.sh
source "${RUNNER_DIR}/lib/submit.sh"
# shellcheck source=lib/verify.sh
source "${RUNNER_DIR}/lib/verify.sh"

# -------- argument parsing --------------------------------------------------

ONLY=""
PARALLEL=1
DRY_RUN=0
VERIFY_ONLY=0
LIST_ONLY=0

usage() { sed -n '2,/^# =\+/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --only)         ONLY="$2"; shift 2 ;;
        --only=*)       ONLY="${1#*=}"; shift ;;
        --parallel)     PARALLEL="$2"; shift 2 ;;
        --parallel=*)   PARALLEL="${1#*=}"; shift ;;
        --dry-run)      DRY_RUN=1; shift ;;
        --verify-only)  VERIFY_ONLY=1; shift ;;
        --list)         LIST_ONLY=1; shift ;;
        -h|--help)      usage 0 ;;
        *)              log_error "Unknown arg: $1"; usage 2 ;;
    esac
done
export PARALLEL DRY_RUN

[[ -f "$CONFIG" ]] || die "Config not found: ${CONFIG}"

# -------- sanity checks -----------------------------------------------------

yq_flavor >/dev/null  # errors out with a helpful message if yq is absent

if [[ "$VERIFY_ONLY" -eq 0 && "$DRY_RUN" -eq 0 && "$LIST_ONLY" -eq 0 ]]; then
    command -v sbatch >/dev/null 2>&1 || die \
        "sbatch not found — this runner requires SLURM. (Try --verify-only on a login node.)"
    command -v squeue >/dev/null 2>&1 || die "squeue not found"
    command -v sacct  >/dev/null 2>&1 || log_warn \
        "sacct not found — parallel-mode state reporting will be degraded"
fi

# -------- build test list ---------------------------------------------------

mapfile -t ALL_IDS < <(list_test_ids "$CONFIG")

select_tests() {
    if [[ -z "$ONLY" ]]; then
        printf '%s\n' "${ALL_IDS[@]}"; return
    fi
    local id selector
    IFS=',' read -r -a selectors <<< "$ONLY"
    for id in "${ALL_IDS[@]}"; do
        for selector in "${selectors[@]}"; do
            # Exact match OR section prefix match ("5.1" matches "5.1a", "5.1b")
            if [[ "$id" == "$selector" ]] || [[ "$id" == "$selector"* ]]; then
                printf '%s\n' "$id"; break
            fi
        done
    done
}

mapfile -t SELECTED < <(select_tests)
if [[ ${#SELECTED[@]} -eq 0 ]]; then
    die "No tests matched --only='${ONLY}'. Available: ${ALL_IDS[*]}"
fi

if [[ "$LIST_ONLY" -eq 1 ]]; then
    printf 'Configured tests (from %s):\n' "$CONFIG"
    for id in "${ALL_IDS[@]}"; do
        name="$(yq_q ".tests[] | select(.id == \"${id}\") | .name" "$CONFIG")"
        printf '  %-6s  %s\n' "$id" "$name"
    done
    exit 0
fi

log_info "Config:   ${CONFIG}"
log_info "Selected: ${SELECTED[*]}"
log_info "Mode:     $([[ $DRY_RUN == 1 ]] && echo dry-run || ([[ $VERIFY_ONLY == 1 ]] && echo verify-only || echo submit))"
log_info "Parallel: ${PARALLEL}"

# -------- parse every selected test up front -------------------------------

for id in "${SELECTED[@]}"; do
    parse_test "$id" "$CONFIG"
done

# -------- run or verify -----------------------------------------------------

overall_rc=0

if [[ "$VERIFY_ONLY" -eq 1 ]]; then
    for id in "${SELECTED[@]}"; do
        verify_test "$id" || overall_rc=1
    done
elif [[ "$PARALLEL" == "1" ]]; then
    # Sequential: submit with sbatch --wait, verify immediately after each.
    for id in "${SELECTED[@]}"; do
        if submit_test "$id"; then
            [[ "$DRY_RUN" -eq 1 ]] || verify_test "$id" || overall_rc=1
        else
            overall_rc=1
        fi
    done
else
    # Parallel: submit everything, then wait, then verify.
    batch=()
    for id in "${SELECTED[@]}"; do
        if (( ${#batch[@]} >= PARALLEL )); then
            wait_for_jobs "${batch[@]}" || overall_rc=1
            for bid in "${batch[@]}"; do
                verify_test "$bid" || overall_rc=1
            done
            batch=()
        fi
        submit_test "$id" || overall_rc=1
        batch+=("$id")
    done
    if (( ${#batch[@]} > 0 )); then
        wait_for_jobs "${batch[@]}" || overall_rc=1
        for bid in "${batch[@]}"; do
            verify_test "$bid" || overall_rc=1
        done
    fi
fi

# -------- summary -----------------------------------------------------------

printf '\n%s========== Summary ==========%s\n' "$C_BOLD" "$C_RESET"
pad() { printf '%-6s %-14s %-s\n' "$1" "$2" "$3"; }
pad "ID" "VERIFY" "OUTPUT"
for id in "${SELECTED[@]}"; do
    v="$(cat "${RESULTS_DIR}/${id}.verify"  2>/dev/null || echo '-')"
    o="$(cat "${RESULTS_DIR}/${id}.outdir"  2>/dev/null || echo '-')"
    pad "$id" "$v" "$o"
done
printf '\nLogs: %s\n' "$LOGS_DIR"

exit "$overall_rc"
