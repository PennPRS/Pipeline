# lib/verify.sh — output-existence checks.
# Requires: lib/common.sh sourced first.

# verify_test <test_id>
#
# Verifies that the job produced the expected artifacts:
#
#   - For training tests (5.1 – 5.6):
#       * Locates the output folder by matching *${submissionID}/ under
#         the configured `base` directory.
#       * Checks that at least `min_weight_files` files matching
#         `weight_glob` exist.
#       * Confirms the first matching weight file has at least one
#         non-header row.
#
#   - For the evaluation test (5.7):
#       * Confirms the output directory exists.
#       * Confirms `required_file` (evaluation_results.txt) exists and
#         has > 1 line (header + at least one row).
#
# Returns 0 on pass, 1 on fail.

verify_test() {
    local id="$1"
    local rec="${RESULTS_DIR}/${id}.rec"
    [[ -f "$rec" ]] || { log_error "[${id}] no record file"; return 1; }
    # shellcheck disable=SC1090
    source "$rec"

    local failures=0

    if [[ -n "${EXPECT_DIR:-}" ]]; then
        # Evaluation-style: a single fixed output directory.
        _verify_evaluation_dir "$id" || failures=$((failures + 1))
    else
        # Training-style: discover output dir by submissionID.
        _verify_training_dir "$id" || failures=$((failures + 1))
    fi

    if [[ "$failures" -eq 0 ]]; then
        log_ok "[${id}] verify OK"
        printf 'PASS\n' > "${RESULTS_DIR}/${id}.verify"
        return 0
    else
        log_error "[${id}] verify FAILED"
        printf 'FAIL\n' > "${RESULTS_DIR}/${id}.verify"
        return 1
    fi
}

# --- training outputs -------------------------------------------------------

_verify_training_dir() {
    local id="$1"
    local base="${EXPECT_BASE:?missing EXPECT_BASE}"
    local sid="${TEST_SUBMISSION_ID:?missing TEST_SUBMISSION_ID}"

    # Match any folder under `base` whose name ends with _${sid}.
    # (Pipeline convention: <trait>_<races>_<method>_<submissionID>/)
    local hits=()
    while IFS= read -r -d '' d; do hits+=("$d"); done < <(
        find "$base" -maxdepth 1 -mindepth 1 -type d -name "*_${sid}" -print0 \
            2>/dev/null
    )

    if [[ ${#hits[@]} -eq 0 ]]; then
        log_error "[${id}] no output folder matching *_${sid} under ${base}"
        return 1
    fi
    if [[ ${#hits[@]} -gt 1 ]]; then
        log_warn "[${id}] multiple matches for ${sid} (${hits[*]}); using first"
    fi
    local outdir="${hits[0]}"
    log_info "[${id}] output: ${outdir}"
    printf '%s\n' "$outdir" > "${RESULTS_DIR}/${id}.outdir"

    # Count weight files.
    local glob="${EXPECT_WEIGHT_GLOB:-*.txt}"
    local min="${EXPECT_MIN_WEIGHT_FILES:-1}"
    local -a wfiles=()
    while IFS= read -r -d '' f; do wfiles+=("$f"); done < <(
        find "$outdir" -maxdepth 2 -type f -name "$glob" -print0 2>/dev/null
    )
    if [[ ${#wfiles[@]} -lt "$min" ]]; then
        log_error "[${id}] expected >= ${min} weight file(s), found ${#wfiles[@]}"
        return 1
    fi
    log_info "[${id}] found ${#wfiles[@]} weight file(s) matching ${glob}"

    # Non-empty check: the first file must have >= 1 data row.
    local sample="${wfiles[0]}"
    _verify_weight_file "$id" "$sample" || return 1
    return 0
}

_verify_weight_file() {
    local id="$1" f="$2"
    local rows
    rows="$(awk 'NR>1' "$f" | wc -l | tr -d ' ')"
    if [[ "$rows" -lt 1 ]]; then
        log_error "[${id}] ${f}: no data rows"
        return 1
    fi
    log_info "[${id}] ${f##*/}: ${rows} data rows"
    return 0
}

# --- evaluation outputs -----------------------------------------------------

_verify_evaluation_dir() {
    local id="$1"
    local outdir="${EXPECT_DIR}"
    if [[ ! -d "$outdir" ]]; then
        log_error "[${id}] output dir does not exist: ${outdir}"
        return 1
    fi
    log_info "[${id}] output: ${outdir}"
    printf '%s\n' "$outdir" > "${RESULTS_DIR}/${id}.outdir"

    local required="${outdir%/}/${EXPECT_REQUIRED_FILE:-evaluation_results.txt}"
    if [[ ! -f "$required" ]]; then
        log_error "[${id}] missing ${required}"
        return 1
    fi
    local lines
    lines="$(wc -l < "$required" | tr -d ' ')"
    if [[ "$lines" -lt 2 ]]; then
        log_error "[${id}] ${required} has only ${lines} line(s) (need header + >=1 row)"
        return 1
    fi
    log_info "[${id}] ${required##*/}: ${lines} lines"
    return 0
}
