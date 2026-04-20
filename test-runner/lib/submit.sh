# lib/submit.sh — sbatch submission helpers.
# Requires: lib/common.sh sourced first.

# submit_test <test_id>
#
# Reads all parameters for the test from the parsed-test-record file
# (${RESULTS_DIR}/<id>.rec produced by parse_test()), submits via sbatch,
# and:
#   - in sequential mode:  `sbatch --wait` (blocks, exit status = job status)
#   - in parallel mode:    `sbatch` (returns immediately; records job id)
#
# Writes:
#   ${RESULTS_DIR}/<id>.jobid      (Slurm job id)
#   ${RESULTS_DIR}/<id>.submit.log (stdout of sbatch)
#   ${LOGS_DIR}/<id>_<jobid>.out   (Slurm stdout)
#   ${LOGS_DIR}/<id>_<jobid>.err   (Slurm stderr)
#
# Returns 0 on submission success (or on sequential job success),
# nonzero on sbatch error (or sequential job failure).

submit_test() {
    local id="$1"
    local rec="${RESULTS_DIR}/${id}.rec"
    [[ -f "$rec" ]] || die "No parsed record for ${id} at ${rec}"

    # Source the per-test record — it sets TEST_* and SLURM_* variables.
    # shellcheck disable=SC1090
    source "$rec"

    local name="$TEST_NAME"
    local script="${PennPRS_path}/${TEST_SCRIPT}"
    if [[ "${DRY_RUN:-0}" != "1" && ! -f "$script" ]]; then
        die "Script not found: ${script}"
    fi

    local out="${LOGS_DIR}/${id}_%j.out"
    local err="${LOGS_DIR}/${id}_%j.err"

    local -a sbatch_args=(
        --job-name  "pennprs_${id}"
        --time      "${SLURM_TIME}"
        --nodes     "${SLURM_NODES:-1}"
        --ntasks    "${SLURM_NTASKS:-1}"
        --cpus-per-task "${SLURM_CPUS}"
        --mem-per-cpu   "${SLURM_MEM_PER_CPU}"
        --partition "${SLURM_PARTITION}"
        --qos       "${SLURM_QOS}"
        # -----------------------------------------------------------------
        # SLURM account:
        #   The line below hard-codes the account as "feixue". To use your
        #   own account, replace it with:
        #
        #       --account "${SLURM_ACCOUNT}"
        #
        #   and then:
        #     1) Add `account: <your-account-name>` under `globals.slurm`
        #        in tests.yaml, e.g.
        #
        #            globals:
        #              slurm:
        #                partition: cpu
        #                qos: normal
        #                account: myaccount
        #
        #     2) Have lib/parse.sh emit SLURM_ACCOUNT into the per-test
        #        record (alongside SLURM_PARTITION / SLURM_QOS), e.g.
        #
        #            SLURM_ACCOUNT="$(yq_q '.globals.slurm.account' "$yaml")"
        #            ...
        #            printf 'SLURM_ACCOUNT=%q\n' "$SLURM_ACCOUNT"
        #
        #   Alternatively, for a quick change, just replace "feixue" with
        #   your own account name on the line below.
        # -----------------------------------------------------------------
        -A feixue
        --output    "$out"
        --error     "$err"
    )

    # In sequential (the default) mode, block until the job finishes.
    if [[ "${PARALLEL:-1}" == "1" ]]; then
        sbatch_args+=(--wait)
    fi

    log_info "[${id}] ${name}"
    log_info "[${id}] sbatch ${sbatch_args[*]} ${script} ${TEST_ARGS[*]}"

    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        log_info "[${id}] (dry-run — not submitting)"
        printf 'DRYRUN\n' > "${RESULTS_DIR}/${id}.jobid"
        return 0
    fi

    local submit_log="${RESULTS_DIR}/${id}.submit.log"
    local status=0
    sbatch "${sbatch_args[@]}" "$script" "${TEST_ARGS[@]}" \
        > "$submit_log" 2>&1 || status=$?

    # Extract job id from "Submitted batch job 12345".
    local jobid
    jobid="$(awk '/Submitted batch job/ {print $NF}' "$submit_log" | tail -n1)"
    if [[ -z "$jobid" ]]; then
        log_error "[${id}] sbatch failed — see ${submit_log}"
        return "${status:-1}"
    fi
    printf '%s\n' "$jobid" > "${RESULTS_DIR}/${id}.jobid"

    if [[ "${PARALLEL:-1}" == "1" ]]; then
        # --wait returns the job's exit code in $status.
        if [[ "$status" -ne 0 ]]; then
            log_error "[${id}] job ${jobid} exited with status ${status}"
            return "$status"
        fi
        log_ok "[${id}] job ${jobid} finished (exit 0)"
    else
        log_info "[${id}] submitted as job ${jobid}"
    fi
    return 0
}

# wait_for_jobs <id> [<id> ...]
#
# For parallel mode: blocks until every listed job has left the Slurm queue,
# then returns nonzero if any finished with a non-COMPLETED state.

wait_for_jobs() {
    local ids=("$@")
    local -a jobids=()
    local id jobid
    for id in "${ids[@]}"; do
        jobid="$(cat "${RESULTS_DIR}/${id}.jobid" 2>/dev/null || true)"
        [[ -n "$jobid" && "$jobid" != "DRYRUN" ]] && jobids+=("$jobid")
    done
    [[ ${#jobids[@]} -eq 0 ]] && return 0

    log_info "Waiting for ${#jobids[@]} job(s): ${jobids[*]}"
    # Poll every 30s; squeue returns empty when all are gone.
    while :; do
        local alive
        alive="$(squeue --noheader --jobs="$(IFS=,; echo "${jobids[*]}")" \
                 2>/dev/null | wc -l | tr -d ' ')"
        [[ "$alive" == "0" ]] && break
        sleep 30
    done

    # Collect final states.
    local rc=0
    for id in "${ids[@]}"; do
        jobid="$(cat "${RESULTS_DIR}/${id}.jobid" 2>/dev/null || echo DRYRUN)"
        [[ "$jobid" == "DRYRUN" ]] && continue
        local state
        state="$(sacct -j "$jobid" --noheader --format=State%-20 \
                 | head -n1 | awk '{print $1}')"
        printf '%s\n' "$state" > "${RESULTS_DIR}/${id}.state"
        if [[ "$state" != "COMPLETED" ]]; then
            log_error "[${id}] job ${jobid} ended in state ${state:-UNKNOWN}"
            rc=1
        else
            log_ok "[${id}] job ${jobid} COMPLETED"
        fi
    done
    return "$rc"
}
