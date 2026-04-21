# ============================================================
# GWAS summary statistics QC pipeline with consolidated report
# ============================================================
# Key updates relative to the original QC.R:
#   1. Writes one consolidated QC report across all ancestries.
#   2. Prints the required-column / accepted-name section once.
#   3. Records early termination clearly in the report before stopping.
#   4. Fixes a few reporting / logic issues in the original script:
#      - A1 missing check now checks A1 rather than SNP.
#      - N column is renamed to N (not SE).
#      - The report file is closed once, after the loop or on early termination.
#      - Z-score reporting includes both pass and fail counts.
#      - The extreme-effect-size threshold is reported consistently.
#
# Assumptions:
#   - workdir, trait, PennPRS_path already exist in the environment.
#   - bigreadr::fread2() and readr::write_delim() are available.
# ============================================================


read_input_file <- function(file_path, ...) {
  detect_delim <- function(file) {
    lines <- readLines(file, n = 5, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) == 0) stop("Input GWAS summary data file is empty. Please re-upload file and retry.")
    
    header <- sub("^\ufeff", "", lines[1])  # remove BOM if present
    
    n_tab   <- lengths(regmatches(header, gregexpr("\t", header)))
    n_comma <- lengths(regmatches(header, gregexpr(",", header)))
    
    if (n_tab > 0 && n_tab >= n_comma) {
      return("\t")
    } else if (n_comma > 0) {
      return(",")
    } else {
      return("whitespace")
    }
  }
  
  clean_colnames <- function(dat) {
    colnames(dat) <- sub("^\ufeff", "", colnames(dat))
    colnames(dat) <- trimws(colnames(dat))
    dat
  }
  
  read_sumstats <- function(fi) {
    sep_guess <- detect_delim(fi)
    
    if (sep_guess == "\t") {
      dat <- bigreadr::fread2(fi, sep = "\t")
    } else if (sep_guess == ",") {
      dat <- bigreadr::fread2(fi, sep = ",")
    } else {
      # For space / generic whitespace-delimited files, use read.table(sep = "")
      # because it reliably splits on one-or-more spaces / tabs.
      dat <- utils::read.table(
        fi,
        header = TRUE,
        sep = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = "",
        quote = "",
        fill = TRUE
      )
    }
    
    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
    dat <- clean_colnames(dat)
    
    # Safety fallback: if the file was mistakenly read as one wide column,
    # retry explicitly as generic whitespace-delimited.
    if (ncol(dat) == 1 && grepl("[[:space:],]", colnames(dat)[1])) {
      dat <- utils::read.table(
        fi,
        header = TRUE,
        sep = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = "",
        quote = "",
        fill = TRUE
      )
      dat <- as.data.frame(dat, stringsAsFactors = FALSE)
      dat <- clean_colnames(dat)
    }
    
    dat
  }
  
  # read a plain text file with automatic delimiter detection
  read_plain_file <- function(path, ...) {
    read_sumstats(path)
  }
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install the 'data.table' package first.")
  }
  
  if (!requireNamespace("bigreadr", quietly = TRUE)) {
    stop("Please install the 'bigreadr' package first.")
  }
  
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Please install the 'readr' package first.")
  }
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Please install the 'readxl' package first for .xlsx support.")
  }
  
  fname <- tolower(basename(file_path))
  
  # -------------------------
  # Case 1: .tar.gz
  # -------------------------
  if (grepl("\\.tar\\.gz$", fname)) {
    tmpdir <- tempfile("untar_dir_")
    dir.create(tmpdir)
    
    utils::untar(file_path, exdir = tmpdir)
    extracted_files <- list.files(tmpdir, recursive = TRUE, full.names = TRUE)
    
    supported <- extracted_files[
      grepl("\\.(txt|tsv|csv|gz|gzip)$", tolower(extracted_files))
    ]
    
    if (length(supported) == 0) {
      stop("No supported file found inside tar.gz archive.")
    }
    
    return(read_input_file(supported[1], ...))
  }
  
  # -------------------------
  # Case 2: .zip
  # -------------------------
  if (grepl("\\.zip$", fname)) {
    zip_contents <- utils::unzip(file_path, list = TRUE)
    inside_files <- zip_contents$Name
    
    supported_idx <- grepl("\\.(txt|tsv|csv|gz|gzip)$", tolower(inside_files))
    
    if (!any(supported_idx)) {
      stop("No supported file found inside zip archive.")
    }
    
    target_file <- inside_files[which(supported_idx)[1]]
    tmpdir <- tempfile("zip_extract_")
    dir.create(tmpdir)
    
    extracted <- utils::unzip(file_path, files = target_file, exdir = tmpdir)
    return(read_input_file(extracted[1], ...))
  }
  
  # -------------------------
  # Case 3: .gz or .gzip
  # -------------------------
  if (grepl("\\.(gz|gzip)$", fname)) {
    # fread can usually read gz directly
    return(read_plain_file(file_path, ...))
  }
  
  # -------------------------
  # Case 4: plain-text files
  # -------------------------
  if (grepl("\\.(txt|tsv|csv)$", fname)) {
    return(read_plain_file(file_path, ...))
  }
  
  # -------------------------
  # Case 5: .xlsx
  # -------------------------
  if (grepl("\\.xlsx$", fname)) {
    return(as.data.frame(readxl::read_xlsx(file_path, ...)))
  }
  
  
  stop("Unsupported file format: ", file_path)
}

write_QC_report <- function(ancestries, trait, temp_path, workdir, PennPRS_path){
  # ----------------------------
  # Helper functions
  # ----------------------------
  write_line <- function(..., .file = zz) {
    cat(paste0(...), "\n", file = .file, append = TRUE, sep = "")
  }
  
  write_section <- function(title) {
    write_line("")
    write_line(paste(rep("=", 60), collapse = ""))
    write_line(title)
    write_line(paste(rep("=", 60), collapse = ""))
  }
  
  write_subsection <- function(title) {
    write_line("")
    write_line(title)
    write_line(paste(rep("-", nchar(title)), collapse = ""))
  }
  
  terminate_pipeline <- function(msg, ancestry, pending_ancestries = character()) {
    write_line("")
    write_line("STATUS: TERMINATED EARLY")
    write_line(paste0("Reason: ", msg))
    write_line(paste0("Termination occurred while processing ancestry: ", ancestry))
    if (length(pending_ancestries) > 0) {
      write_line(paste0(
        "The following ancestries were not processed because the job terminated early: ",
        paste(pending_ancestries, collapse = ", ")
      ))
    }
    flush(zz)
    close(zz)
    stop(msg, call. = FALSE)
  }
  
  fmt_num <- function(x) {
    if (length(x) == 0 || is.null(x) || all(is.na(x))) return("NA")
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  }
  
  
  find_trait_file <- function(directory, trait_name, full.names = TRUE) {
    files <- list.files(directory, full.names = full.names)
    matches <- files[grepl(trait_name, basename(files), fixed = TRUE)]
    return(matches)
  }
  
  
  # ----------------------------
  # Report file
  # ----------------------------
  filen <- paste0(workdir, "QC_report.txt")
  if (file.exists(filen)) file.remove(filen)
  file.create(filen)
  zz <- file(filen, open = "wt")
  on.exit({
    if (exists("zz") && isOpen(zz)) close(zz)
  }, add = TRUE)
  # on.exit({
  #   if (exists("zz")) {
  #     try(close(zz), silent = TRUE)
  #   }
  # }, add = TRUE)
  
  # ----------------------------
  # Report header
  # ----------------------------
  write_line("GWAS SUMMARY STATISTICS QC REPORT")
  write_line(paste0("Trait name: ", trait))
  write_line(paste0("Ancestries detected: ", paste(ancestries, collapse = ", ")))
  write_line(paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  write_section("1. OVERVIEW")
  write_line("This report summarizes automated quality-control checks applied to GWAS summary data files.")
  write_line("")
  write_line("Please note:")
  write_line("")
  write_line("  * We support the following file formats: .txt, .tsv, .csv, .xlsx, .zip, .gz, .gzip, and .tar.gz.")
  # write_line("")
  write_line("  * We accept file separated by tab, comma, or space.")
  # write_line("")
  write_line("  * The size of each single file cannot exceed: 800MB.")
  write_line("")
  write_line("If the pipeline terminates early, the report will explicitly record the failure point and note which QC blocks did not run.")
  write_line("")
  write_line("Please contact us at pennprs@googlegroups.com or report issues at https://groups.google.com/g/pennprs if you encounter further issues.")
  
  
  write_section("2. REQUIRED COLUMNS AND COLUMN NAMES ALLOWED (CASE-INSENSITIVE)")
  write_line("The following columns are required by PennPRS.")
  write_line("Input matching is case-insensitive, and any of the accepted names below may be used.")
  write_line("")
  write_line("  - CHR | #CHROM | Chromosome: chromosome.")
  write_line("  - SNP | ID | RSID | Snpid | Snpid_UKB | RS | Rs_id: SNP RSID (format: rsXXXX).")
  write_line("  - A1 | ALT | Effect_allele | allele1 | allele_1 | alt_allele | EA: effect allele.")
  write_line("  - A2 | REF | Allele2 | OMITTED | allele0 | allele_0 | Allele_2 | Ref_allele | Other_allele | NEA: alternative/other allele.")
  write_line("  - MAF | AF | AAF | AF1 | A1_FREQ | A1FREQ | Effect_allele_frequency | Eaf | FRQ | FRQ_U | F_U: minor allele frequency.")
  write_line("  - BETA | OR | logOR | Effect: SNP effect size.")
  write_line("  - SE | LOG(OR)_SE | Stderr | Std_Error | Stderr_Beta | SE_Beta | Stderr_B: standard error of BETA.")
  write_line("  - P | LOG10P | LP | Pvalue | P_value | Pval | P_val | GC_Pvalue: p-value.")
  write_line("  - N | OBS_CT | Neff | N_eff: GWAS sample size (default column for sample size).")
  write_line("  - Ncase | N_case: number of cases for binary trait (optional).")
  write_line("  - Ncontrol | N_control: number of controls for binary trait (optional).")
  write_line("")
  write_line("Column information and notes:")
  write_line("  - SNP: if rsID information is missing, impute rsID with position information using reference genotype data from the same genome build.")
  write_line("  - BETA: if OR is provided for a binary trait, the QC pipeline will convert OR to logOR.")
  write_line("  - MAF: minor allele frequency, effect-allele frequency, or the frequency of either A1 or A2.")
  write_line("  - SE: if only the SE of OR is available, please compute SE(OR)/OR as the approximate SE(logOR) and input it as the SE column.")
  write_line("  - For quantitative traits, N is the total sample size. For binary traits, N is the effective sample size: 4 / (1 / N_control + 1 / N_case).")
  write_line("  - If N is not provided for a binary trait, both Ncase and Ncontrol should be provided.")
  
  # Track ancestry status for optional final summary
  ancestry_status <- setNames(rep("NOT RUN", length(ancestries)), ancestries)
  
  # ----------------------------
  # Main ancestry loop
  # ----------------------------
  for (i.ans in seq_along(ancestries)) {
    ancestry <- ancestries[i.ans]
    pending_ancestries <- ancestries[(min(i.ans + 1, length(ancestries))):length(ancestries)]
    if (length(pending_ancestries) == 1 && is.na(pending_ancestries)) pending_ancestries <- character()
    
    trait_name <- paste0(ancestry, "_", trait)
    ancestry_status[ancestry] <- "STARTED"
    
    if (length(ancestries) == 1) write_section(paste0("3. QC FOR ", ancestry))
    if (length(ancestries) > 1) write_section(paste0("3.", i.ans, " QC FOR ancestry: ", ancestry))
    
    # Read input
    file.nm = find_trait_file(temp_path, trait_name, full.names = TRUE)
    if (length(file.nm) == 0){
      ancestry_status[ancestry] <- "FAILED"
      terminate_pipeline(
        msg = paste0("Input GWAS summary data with file name ", trait_name, " not detected.\nPlease check if the input GWAS summary data has been saved as ", trait_name, "{.txt, .csv, .tsv, .xlsx, .gz, .zip, .gzip, .tar.gz} in the input GWAS folder."),
        ancestry = ancestry,
        pending_ancestries = pending_ancestries
      )
    } else {
      tp = strsplit(file.nm[1], split = '/')[[1]]; file.nm1 = tp[length(tp)]; rm(tp)
      sumraw <- read_input_file(file.nm[1])
      n_input <- nrow(sumraw)
      write_line(paste0("Input file: ", file.nm1))
      write_line(paste0("Total # of input variants: ", fmt_num(n_input)))
      
      # # ----------------------------
      # # (0) Remove irrelevant columns with too many NA values:
      # # ----------------------------
      # 
      # p.missing <- 0.2
      # write_subsection(paste0("(0) Remove irrelevant columns with more than ", p.missing, " columns having NA/NaN values"))
      # na_frac_by_col <- sapply(sumraw, function(x) mean(is.na(x)))
      # cols_remove_high_na <- names(na_frac_by_col)[na_frac_by_col > p.missing]
      # n_cols_remove_high_na <- length(cols_remove_high_na)
      # if (n_cols_remove_high_na == 0) {
      #   write_line("  - No columns with >20% rows having NA/NaN values were detected.")
      # } else {
      #   sumraw <- sumraw[, !(colnames(sumraw) %in% cols_remove_high_na), drop = FALSE]
      #   write_line(paste0("  - ", n_cols_remove_high_na, " columns with >", paste0(100 * p.missing, "%"), " rows having NA/NaN values were removed."))
      #   write_line(paste0("  - Removed columns: ", paste(cols_remove_high_na, collapse = ", ")))
      #   write_line(paste0("  - Columns remaining after this step: ", ncol(sumraw)))
      # }
      
      
      # ----------------------------
      # (1) Check required GWAS column names
      # ----------------------------
      write_subsection("(1) Check required GWAS column names")
      
      col.chr <- which(tolower(colnames(sumraw)) %in% tolower(c("CHR", "#CHROM", "Chromosome")))
      if (length(col.chr) > 0) {
        original_name <- colnames(sumraw)[col.chr[1]]
        colnames(sumraw)[col.chr[1]] <- "CHR"
        write_line(paste0("  - Column '", original_name, "' detected for chromosome information."))
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required chromosome column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.snp <- which(tolower(colnames(sumraw)) %in% tolower(c("SNP", "ID", "RSID", "Snpid", "Snpid_UKB", "RS", "Rs_id")))
      if (length(col.snp) > 0) {
        original_name <- colnames(sumraw)[col.snp[1]]
        colnames(sumraw)[col.snp[1]] <- "SNP"
        write_line(paste0("  - Column '", original_name, "' detected for SNP rsID information."))
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required SNP rsID column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.a1 <- which(tolower(colnames(sumraw)) %in% tolower(c("A1", "ALT", "Effect_allele", "allele1", "allele_1", "alt_allele", "EA")))
      if (length(col.a1) > 0) {
        original_name <- colnames(sumraw)[col.a1[1]]
        colnames(sumraw)[col.a1[1]] <- "A1"
        write_line(paste0("  - Column '", original_name, "' detected for A1 (effect allele)."))
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required A1 (effect allele) column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.a2 <- which(tolower(colnames(sumraw)) %in% tolower(c("A2", "REF", "allele0", "allele_0", "OMITTED", "Allele2", "Allele_2", "Ref_allele", "Other_allele", "NEA")))
      if (length(col.a2) > 0) {
        original_name <- colnames(sumraw)[col.a2[1]]
        colnames(sumraw)[col.a2[1]] <- "A2"
        write_line(paste0("  - Column '", original_name, "' detected for A2 (other allele)."))
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required A2 (other allele) column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.maf <- which(tolower(colnames(sumraw)) %in% tolower(c("MAF", 'AF', 'AAF', 'AF1', 'A1_FREQ', 'A1FREQ', "Effect_allele_frequency", "Eaf", "FRQ", "FRQ_U", "F_U")))
      if (length(col.maf) > 0) {
        original_name <- colnames(sumraw)[col.maf[1]]
        colnames(sumraw)[col.maf[1]] <- "MAF"
        write_line(paste0("  - Column '", original_name, "' detected for MAF/allele frequency information."))
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required MAF column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.beta <- which(tolower(colnames(sumraw)) %in% tolower(c("BETA", "OR", "logOR", "Effect")))
      if (length(col.beta) > 0) {
        original_name <- colnames(sumraw)[col.beta[1]]
        if (tolower('BETA') %in% tolower(colnames(sumraw)[col.beta])) original_name <- 'BETA'
        if (tolower('logOR') %in% tolower(colnames(sumraw)[col.beta])) original_name <- 'logOR'
        if (tolower('Effect') %in% tolower(colnames(sumraw)[col.beta])) original_name <- 'Effect'
        col.beta = which(tolower(colnames(sumraw)) == tolower(original_name))
        if (tolower(original_name) == "or") {
          col.beta = which(tolower(colnames(sumraw)) == 'or')
          sumraw[, col.beta] <- log(as.numeric(sumraw[, col.beta]))
          write_line(paste0("  - Column '", original_name, "' detected for BETA (effect size). Converted to logOR."))
        } else {
          write_line(paste0("  - Column '", original_name, "' detected for BETA (effect size)."))
        }
        colnames(sumraw)[col.beta[1]] <- "BETA"
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required BETA / OR / logOR / Effect column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.se <- which(tolower(colnames(sumraw)) %in% tolower(c("SE", "LOG(OR)_SE", "Stderr", "Std_Error", "Stderr_Beta", "SE_Beta", "Stderr_B")))
      if (length(col.se) > 0) {
        original_name <- colnames(sumraw)[col.se[1]]
        colnames(sumraw)[col.se[1]] <- "SE"
        write_line(paste0("  - Column '", original_name, "' detected for SE information."))
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required SE column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.p <- which(tolower(colnames(sumraw)) %in% tolower(c("P", "LOG10P", "LP", "Pvalue", "P_value", "Pval", "P_val", "GC_Pvalue")))
      if (length(col.p) > 0) {
        original_name <- colnames(sumraw)[col.p[1]]
        colnames(sumraw)[col.p[1]] <- "P"
        write_line(paste0("  - Column '", original_name, "' detected for p-value information."))
      } else {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Required p-value column not detected.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      col.n <- which(tolower(colnames(sumraw)) %in% tolower(c("N", "OBS_CT", "Neff", "N_eff")))
      if (length(col.n) > 0) {
        original_name <- colnames(sumraw)[col.n[1]]
        colnames(sumraw)[col.n[1]] <- "N"
        write_line(paste0("  - Column '", original_name, "' detected for sample-size information."))
      } else {
        col.ncase <- which(tolower(colnames(sumraw)) %in% tolower(c("Ncase", "N_case")))
        col.ncontrol <- which(tolower(colnames(sumraw)) %in% tolower(c("Ncontrol", "N_control")))
        if ((length(col.ncase) > 0) & (length(col.ncontrol) > 0)) {
          original_name.ncase <- colnames(sumraw)[col.ncase[1]]
          original_name.ncontrol <- colnames(sumraw)[col.ncontrol[1]]
          sumraw$N = 4 / (1 / sumraw[, col.ncase[1]] + 1 / sumraw[, col.ncontrol[1]])
          write_line(paste0("  - Columns '", original_name.ncase, " and ", original_name.ncontrol, "' detected for information on Ncases and Ncontrols."))
        } 
        if ((length(col.ncase) == 0) | (length(col.ncontrol) == 0)) {
          ancestry_status[ancestry] <- "FAILED"
          msg <- if ((length(col.ncase) == 0) & (length(col.ncontrol) == 0)) {
            "Required N column or Ncase & Ncontrol information not detected.\nYour data should include either two separate columns for the number of cases and controls (Ncase & Ncontrol) or one column for the effective sample size (Neff). Please ensure your data contains these columns and input their column names."
          } else if ((length(col.ncase) > 0) & (length(col.ncontrol) == 0)) {
            "Required N column not detected. Ncase column detected but Ncontrol column is missing.\nYour data should include either two separate columns for the number of cases and controls (Ncase & Ncontrol) or one column for the effective sample size (Neff). Please ensure your data contains these columns and input their column names."
          } else {
            "Required N column not detected. Ncontrol column detected but Ncase column is missing.\nYour data should include either two separate columns for the number of cases and controls (Ncase & Ncontrol) or one column for the effective sample size (Neff). Please ensure your data contains these columns and input their column names."
          }
          terminate_pipeline(
            msg = msg,
            ancestry = ancestry,
            pending_ancestries = pending_ancestries
          )
        }
      }
      
      write_line("    Column name check and cleaning completed!")
      
      # Convert relevant columns
      sumraw$BETA <- as.numeric(sumraw$BETA)
      sumraw$SE   <- as.numeric(sumraw$SE)
      sumraw$MAF  <- as.numeric(sumraw$MAF)
      sumraw$P    <- as.numeric(sumraw$P)
      sumraw$N    <- as.numeric(sumraw$N)
      
      sumraw <- sumraw[, c('CHR', 'SNP', 'A1', 'A2', 'MAF', 'BETA', 'SE', 'P', 'N')]
      
      # ----------------------------
      # (2) Keep SNPs in HapMap3 and chromosomes 1-22
      # ----------------------------
      write_subsection("(2) Restrict to HapMap3 SNPs on autosomes (CHR 1 - 22)")
      hm3 <- bigreadr::fread2(paste0(PennPRS_path, "data/hapmap3rsid.txt"))[, 1]
      nochr <- which(!sumraw$CHR %in% c(1:22))
      nohm3 <- which(!sumraw$SNP %in% hm3)
      rm.indx0 <- unique(c(nochr, nohm3))
      write_line(paste0("  - Variants outside chromosomes 1-22 detected and removed: ", fmt_num(length(nochr))))
      write_line(paste0("  - Variants not found in HapMap 3 detected and removed: ", fmt_num(length(nohm3))))
      if (length(rm.indx0) > 0) sumraw <- sumraw[-rm.indx0, ]
      write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      
      # ----------------------------
      # (3) Missing-value check
      # ----------------------------
      write_subsection("(3) Check SNPs with missing data")
      n.na <- sum(!complete.cases(sumraw))
      if (n.na == 0) {
        write_line("  - No variants with missing required information were detected.")
      } else {
        sumraw <- sumraw[complete.cases(sumraw), ]
        write_line(paste0("  - Variants with missing required information removed: ", fmt_num(n.na)))
        write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      }
      
      # ----------------------------
      # (4) Allele-frequency range check
      # ----------------------------
      write_subsection("(4) Remove SNPs with MAF < 0.01")
      wrong.af <- which((sumraw$MAF > 0.99) | (sumraw$MAF < 0.01))
      n.wrong.af <- length(wrong.af)
      if (n.wrong.af == 0) {
        write_line("  - No variants with MAF < 0.01 were detected.")
      } else {
        sumraw <- sumraw[-wrong.af, ]
        write_line(paste0("  - Variants with MAF < 0.01 removed: ", fmt_num(n.wrong.af)))
        write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      }
      
      # ----------------------------
      # (5) Remove duplicated SNP IDs
      # ----------------------------
      write_subsection("(5) Remove duplicated SNP IDs")
      dup.id <- which(duplicated(sumraw$SNP))
      n.dup <- length(dup.id)
      if (n.dup == 0) {
        write_line("  - No duplicated SNP IDs were detected.")
      } else {
        sumraw <- sumraw[-dup.id, ]
        write_line(paste0("  - Rows with duplicated SNP IDs removed: ", fmt_num(n.dup)))
        write_line("  - For each duplicated SNP, only the first row was retained.")
        write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      }
      
      # ----------------------------
      # (6) Z-score check
      # ----------------------------
      write_subsection("(6) Check SNP z-scores")
      chi2_thr <- 30
      z_ok <- which(abs(sumraw$BETA / sumraw$SE) < sqrt(chi2_thr))
      n_ok <- length(z_ok)
      n_fail <- nrow(sumraw) - n_ok
      write_line(paste0("  - Variants with |z| < sqrt(", chi2_thr, "): ", fmt_num(n_ok)))
      write_line(paste0("  - Variants exceeding this threshold: ", fmt_num(n_fail)))
      if (n_ok < 5) {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = paste0("Fewer than 5 SNPs have |z| < sqrt(", chi2_thr, "), suggesting issues with the input GWAS data."),
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      }
      
      # ----------------------------
      # (7) Problematic BETA
      # ----------------------------
      write_subsection("(7) Check SNPs with problematic BETA values (|BETA| > 1000)")
      beta.thr <- 1e3
      rm.indx1 <- which(abs(sumraw$BETA) > beta.thr)
      if (length(rm.indx1) == 0) {
        write_line(paste0("  - No variants with |BETA| > ", beta.thr, " were detected."))
      } else {
        sumraw <- sumraw[-rm.indx1, ]
        write_line(paste0("  - Variants with |BETA| > ", beta.thr, " removed: ", fmt_num(length(rm.indx1))))
        write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      }
      
      # ----------------------------
      # (8) Extremely large effect sizes
      # ----------------------------
      write_subsection("(8) Check SNPs with extremely large effect sizes (|z-score^2| > 10000)")
      chi2_large_thr <- 1e4
      rm.indx2 <- which((sumraw$BETA / sumraw$SE)^2 > chi2_large_thr)
      if (length(rm.indx2) == 0) {
        write_line(paste0("  - No variants with z-score^2 > ", chi2_large_thr, " were detected."))
      } else {
        sumraw <- sumraw[-rm.indx2, ]
        write_line(paste0("  - Variants with z-score^2 > ", chi2_large_thr, " removed: ", fmt_num(length(rm.indx2))))
        write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      }
      
      # ----------------------------
      # (9) Problematic p-values
      # ----------------------------
      write_subsection("(9) Check SNPs with problematic p-values")
      rm.indx3 <- which((sumraw$P > 1) | (sumraw$P < 0))
      if (length(rm.indx3) == 0) {
        write_line("  - No variants with p-value < 0 or > 1 were detected.")
      } else {
        sumraw <- sumraw[-rm.indx3, ]
        write_line(paste0("  - Variants with p-value < 0 or > 1 removed: ", fmt_num(length(rm.indx3))))
        write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      }
      
      # ----------------------------
      # (10) Zero SE
      # ----------------------------
      write_subsection("(10) Check SNPs with zero or negative SE")
      rm.indx4 <- which(sumraw$SE <= 0)
      if (length(rm.indx4) == 0) {
        write_line("  - No variants with SE =/< 0 were detected.")
      } else {
        sumraw <- sumraw[-rm.indx4, ]
        write_line(paste0("  - Variants with SE =/< 0 removed: ", fmt_num(length(rm.indx4))))
        write_line(paste0("  - Variants remaining after this step: ", fmt_num(nrow(sumraw))))
      }
      
      # ----------------------------
      # Write cleaned file
      # ----------------------------
      n.rm = n.na + n.wrong.af + n.dup + length(unique(c(rm.indx0, rm.indx1, rm.indx2, rm.indx3, rm.indx4)))
      
      if (nrow(sumraw) == 0) {
        ancestry_status[ancestry] <- "FAILED"
        terminate_pipeline(
          msg = "Zero SNPs remaining after QC.",
          ancestry = ancestry,
          pending_ancestries = pending_ancestries
        )
      } else {
        readr::write_delim(sumraw, file = paste0(workdir, "sumdata/", trait_name, ".txt"), delim = "\t")
        write_line(paste0("  - Total problematic variants removed: ", fmt_num(n.rm)))
        write_line(paste0("  - Variants remaining after all QC steps: ", fmt_num(nrow(sumraw))))
      }
      
      
      ancestry_status[ancestry] <- "COMPLETED"
      write_line("")
      write_line(paste0("QC completed successfully for ", trait_name, ".txt"))
      flush(zz)
    }
  }
  
  # ----------------------------
  # Final summary
  # ----------------------------
  write_section("4. FINAL STATUS SUMMARY BY ANCESTRY")
  for (anc in names(ancestry_status)) {
    write_line(paste0("  - ", anc, ": ", ancestry_status[anc]))
  }
  
  write_section("5. OPTIONAL ADDITIONAL QC CHECKS TO CONSIDER")
  write_line("  - Review inflation or deflation in test statistics using genomic inflation factor and QQ plot diagnostics.")
  write_line("  - Verify that sample-size mapping is correct when multiple candidate columns (e.g., N, Ncase, Ncontrol) exist.")
  write_line("  - Flag variants with unusually small sample size relative to the study median.")
  
  write_line("")
  write_line("QC completed for input GWAS summary statistics.")
  flush(zz)
}




txt_to_html <- function(input_txt, output_html, title = "Text Report") {
  if (!file.exists(input_txt)) {
    stop("Input file does not exist: ", input_txt)
  }
  
  lines <- readLines(input_txt, warn = FALSE, encoding = "UTF-8")
  
  escape_html <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  
  text_html <- paste(escape_html(lines), collapse = "\n")
  
  html <- paste0(
    '<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>', title, '</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      margin: 40px;
      background: #f7f7f7;
      color: #222;
      line-height: 1.6;
    }
    .container {
      max-width: 1000px;
      margin: auto;
      background: white;
      padding: 30px;
      border-radius: 10px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.08);
    }
    h1 {
      margin-top: 0;
      border-bottom: 2px solid #e5e5e5;
      padding-bottom: 10px;
      font-size: 30px;
    }
    pre {
      white-space: pre-wrap;
      word-wrap: break-word;
      background: #fafafa;
      padding: 20px;
      border-radius: 8px;
      border: 1px solid #e5e5e5;
      overflow-x: auto;
      font-size: 14px;
    }
  </style>
</head>
<body>
  <div class="container">
    <h1>', title, '</h1>
    <pre>', text_html, '</pre>
  </div>
</body>
</html>'
  )
  
  writeLines(html, output_html, useBytes = TRUE)
  
  message("HTML file written to: ", output_html)
}




convert_to_pgs_tsv <- function(input_file,
                               output_file,
                               pgs_name = "metaGRS_CAD",
                               pgs_id = "metaGRS_CAD",
                               trait_reported = "Coronary artery disease",
                               genome_build = "GRCh37",
                               column_mapping = c(
                                 chr_name = "chr_name",
                                 chr_position = "chr_position",
                                 effect_allele = "effect_allele",
                                 other_allele = "other_allele",
                                 effect_weight = "effect_weight"
                               )) {
  
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  # Read the original tab-delimited file
  dat <- bigreadr::fread2(input_file)
  
  # Check that requested source columns exist
  missing_cols <- setdiff(unname(column_mapping), colnames(dat))
  if (length(missing_cols) > 0) {
    stop("These input columns are missing: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Keep and rename columns in the required order
  dat_out <- dat[, unname(column_mapping), drop = FALSE]
  colnames(dat_out) <- names(column_mapping)
  
  # Write metadata header
  header_lines <- c(
    paste0("#pgs_name=", pgs_name),
    paste0("#pgs_id=", pgs_id),
    paste0("#trait_reported=", trait_reported),
    paste0("#genome_build=", genome_build)
  )
  
  writeLines(header_lines, con = output_file)
  
  # Append the table in TSV format
  suppressWarnings(utils::write.table(dat_out,
                     file = output_file,
                     sep = "\t",
                     row.names = FALSE,
                     col.names = TRUE,
                     quote = FALSE,
                     append = TRUE))
  
  message("Output written to: ", output_file)
}
