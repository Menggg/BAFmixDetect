#!/usr/bin/env Rscript

## Pure base-R all-in-one BAF regression pipeline.
## Supports Illumina FinalReport and VCF inputs, multiple samples, optional BAF recomputation,
## automatic or user-specified delimiters for FinalReport/frequency files, and combined summary output.

VERSION <- "1.0.0-baseR-FinalReport-VCF"

stopf <- function(...) stop(sprintf(...), call. = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

is_flag <- function(x) grepl("^--", x)

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(as.character(x))
  if (y %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (y %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stopf("Cannot parse logical value: %s", x)
}

opt_has_value <- function(x) {
  if (is.null(x)) return(FALSE)
  if (!length(x)) return(FALSE)
  x <- as.character(x[[1L]])
  if (is.na(x)) return(FALSE)
  if (!nzchar(x)) return(FALSE)
  !x %in% c("None", "none", "NULL", "null", "NA")
}

opt_get <- function(opts, key, default = NULL) {
  if (!(key %in% names(opts))) return(default)
  val <- opts[[key]]
  if (is.null(val)) return(default)
  val
}

.qc_warning_seen <- new.env(parent = emptyenv())

emit_qc_warning <- function(options, key, message) {
  key <- as.character(key)
  if (exists(key, envir = .qc_warning_seen, inherits = FALSE)) return(invisible(FALSE))
  assign(key, TRUE, envir = .qc_warning_seen)
  line <- paste0("WARNING: ", message)
  cat(line, "\n", sep = "")
  warning_log <- opt_get(options, "warning_log", NULL)
  if (opt_has_value(warning_log)) {
    cat(paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\t", line, "\n"),
        file = warning_log, append = TRUE)
  }
  invisible(TRUE)
}

parse_cli <- function(args, defaults = list(), bool_flags = character()) {
  opts <- defaults
  pos <- character()
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      keyval <- sub("^--", "", a)
      inline <- NULL
      eq <- regexpr("=", keyval, fixed = TRUE)[[1L]]
      if (!identical(eq, -1L)) {
        inline <- substring(keyval, eq + 1L)
        keyval <- substring(keyval, 1L, eq - 1L)
      }
      key <- keyval
      if (startsWith(key, "no") && substring(key, 3L) %in% bool_flags && is.null(inline)) {
        opts[[substring(key, 3L)]] <- FALSE
        i <- i + 1L
        next
      }
      if (key %in% bool_flags) {
        if (!is.null(inline)) {
          opts[[key]] <- as_bool(inline)
          i <- i + 1L
        } else if (i < length(args) && !is_flag(args[[i + 1L]]) &&
                   tolower(args[[i + 1L]]) %in% c("true", "false", "t", "f", "1", "0", "yes", "no", "y", "n")) {
          opts[[key]] <- as_bool(args[[i + 1L]])
          i <- i + 2L
        } else {
          opts[[key]] <- TRUE
          i <- i + 1L
        }
        next
      }
      if (!is.null(inline)) {
        opts[[key]] <- inline
        i <- i + 1L
        next
      }
      if (i == length(args) || is_flag(args[[i + 1L]])) stopf("Option --%s requires a value", key)
      opts[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      pos <- c(pos, a)
      i <- i + 1L
    }
  }
  list(options = opts, args = pos)
}

normalize_sep <- function(sep) {
  if (is.null(sep)) return(NULL)
  sep <- as.character(sep)
  if (!length(sep) || is.na(sep) || sep %in% c("", "auto", "AUTO", "Auto", "NULL", "null", "None", "none", "NA")) return(NULL)
  if (sep %in% c("\\t", "tab", "TAB", "Tab", "tsv", "TSV")) return("\t")
  if (sep %in% c(",", "comma", "COMMA", "csv", "CSV")) return(",")
  if (sep %in% c("space", "SPACE", "whitespace", "white", "ws")) return(NULL)
  sep
}

split_line <- function(line, sep = NULL) {
  sep <- normalize_sep(sep)
  if (is.null(sep)) {
    z <- trimws(line)
    if (!nzchar(z)) return(character())
    return(strsplit(z, "\\s+", perl = TRUE)[[1L]])
  }
  strsplit(line, sep, fixed = TRUE)[[1L]]
}

count_fixed <- function(x, pattern) {
  m <- gregexpr(pattern, x, fixed = TRUE)[[1L]]
  if (identical(m, -1L)) 0L else length(m)
}

guess_sep_from_line <- function(line) {
  if (count_fixed(line, "\t") >= 2L) return("\t")
  if (count_fixed(line, ",") >= 2L) return(",")
  NULL
}

read_nonempty_head <- function(file, n = 80L) {
  con <- file(file, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

require_file <- function(path, label = "file") {
  if (!opt_has_value(path) || !file.exists(path) || file.info(path)$size <= 0) {
    stopf("ERROR: %s not found or empty: %s", label, as.character(path))
  }
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  if (!nzchar(x)) "sample" else x
}

make_key_env <- function(items) {
  items <- unique(as.character(items))
  items <- items[!is.na(items) & nzchar(items)]
  e <- new.env(parent = emptyenv(), hash = TRUE, size = max(29L, length(items) * 2L))
  for (v in items) assign(v, TRUE, envir = e)
  e
}

load_filter_values <- function(file, col = 1L, sep = NULL) {
  if (!opt_has_value(file)) return(character())
  lines <- readLines(file, warn = FALSE)
  col <- as.integer(col)
  vals <- character()
  if (length(lines)) {
    vals <- vapply(lines, function(line) {
      a <- split_line(line, sep)
      if (length(a) >= col) a[[col]] else NA_character_
    }, character(1L), USE.NAMES = FALSE)
  }
  unique(vals[!is.na(vals) & nzchar(vals)])
}

apply_name_filters <- function(values, include_file = NULL, include_col = 1L, include_sep = NULL,
                               exclude_file = NULL, exclude_col = 1L, exclude_sep = NULL) {
  keep <- rep(TRUE, length(values))
  if (opt_has_value(include_file)) {
    inc <- make_key_env(load_filter_values(include_file, include_col, include_sep))
    keep <- keep & vapply(as.character(values), exists, logical(1L), envir = inc, inherits = FALSE)
  }
  if (opt_has_value(exclude_file)) {
    exc <- make_key_env(load_filter_values(exclude_file, exclude_col, exclude_sep))
    keep <- keep & !vapply(as.character(values), exists, logical(1L), envir = exc, inherits = FALSE)
  }
  keep
}

parse_freq_cols <- function(freqcols) {
  x <- as.integer(strsplit(as.character(freqcols), ",", fixed = TRUE)[[1L]])
  if (length(x) < 2L || any(is.na(x[1:2]))) stopf("--freqcols must look like 1,2")
  x[1:2]
}

load_frequency <- function(freqfile, freqcols = "1,2", freqsep = NULL, min_maf = 0) {
  require_file(freqfile, "freqfile")
  cols <- parse_freq_cols(freqcols)
  sep <- normalize_sep(freqsep)
  dat <- read.table(freqfile, header = TRUE, sep = sep %||% "", quote = "", comment.char = "",
                    fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(dat) < max(cols)) stopf("freqfile has %d columns, but --freqcols requests column %d", ncol(dat), max(cols))
  marker <- as.character(dat[[cols[[1L]]]])
  maf <- suppressWarnings(as.numeric(dat[[cols[[2L]]]]))
  ok <- !is.na(marker) & nzchar(marker) & !is.na(maf) & maf > 0
  marker <- marker[ok]
  maf <- maf[ok]
  maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  min_maf <- suppressWarnings(as.numeric(min_maf))
  if (length(min_maf) != 1L || is.na(min_maf) || min_maf < 0 || min_maf > 0.5) {
    stopf("--min_maf must be a number between 0 and 0.5")
  }
  if (min_maf > 0) {
    keep_maf <- maf >= min_maf
    marker <- marker[keep_maf]
    maf <- maf[keep_maf]
  }
  out <- maf
  names(out) <- marker
  out
}


## ---------------- conservative SNP QC ----------------

num_option <- function(options, key, default = 0, min_value = -Inf, max_value = Inf) {
  x <- suppressWarnings(as.numeric(opt_get(options, key, default)))
  if (length(x) != 1L || is.na(x) || x < min_value || x > max_value)
    stopf("--%s must be between %s and %s", key, min_value, max_value)
  x
}

apply_snp_qc <- function(sample, marker, baf, abgeno, maf, options,
                         gencall = NULL, gcscore = NULL, intensity = NULL, chrom = NULL, cnv = NULL) {
  n <- length(marker)
  stopifnot(length(baf) == n, length(abgeno) == n, length(maf) == n)
  active <- rep(TRUE, n)
  counts <- list(sample = as.character(sample), N_total = n)

  remove_step <- function(flag, name) {
    flag <- active & flag
    counts[[name]] <<- sum(flag, na.rm = TRUE)
    active[flag] <<- FALSE
  }

  remove_step(is.na(marker) | !nzchar(marker), "N_missing_marker")

  if (isTRUE(opt_get(options, "remove_duplicates", TRUE))) {
    dup <- duplicated(marker)
    remove_step(dup, "N_duplicate_removed")
  } else counts$N_duplicate_removed <- 0L

  valid_baf_only <- isTRUE(opt_get(options, "valid_baf_only", TRUE))
  if (valid_baf_only)
    remove_step(!is.finite(baf) | baf < 0 | baf > 1, "N_invalid_baf_removed")
  else {
    remove_step(is.na(baf), "N_invalid_baf_removed")
  }

  remove_step(is.na(maf), "N_missing_maf_removed")

  if (isTRUE(opt_get(options, "remove_input_cnv", TRUE))) {
    if (is.null(cnv)) cnv <- rep(FALSE, n)
    remove_step(is.na(cnv) | cnv, "N_CNV_removed")
  } else counts$N_CNV_removed <- 0L

  counts$N_nocall <- sum(active & abgeno == 3L, na.rm = TRUE)
  counts$N_invalid_genotype <- sum(active & !(abgeno %in% 0:3), na.rm = TRUE)
  remove_step(!(abgeno %in% 0:3) | is.na(abgeno), "N_invalid_genotype_removed")

  min_gencall <- num_option(options, "min_gencall", 0, 0, Inf)
  if (min_gencall > 0) {
    if (is.null(gencall)) stopf("--min_gencall > 0 requires FinalReport column: %s", opt_get(options, "colgencall", "GenCall Score"))
    remove_step(is.na(gencall) | gencall < min_gencall, "N_low_gencall_removed")
  } else counts$N_low_gencall_removed <- 0L

  min_gcscore <- num_option(options, "min_gcscore", 0, 0, Inf)
  if (min_gcscore > 0) {
    if (is.null(gcscore)) stopf("--min_gcscore > 0 requires FinalReport column: %s", opt_get(options, "colgcscore", "GC Score"))
    remove_step(is.na(gcscore) | gcscore < min_gcscore, "N_low_gcscore_removed")
  } else counts$N_low_gcscore_removed <- 0L

  min_intensity <- num_option(options, "min_intensity", 0, 0, Inf)
  if (min_intensity > 0) {
    if (is.null(intensity)) stopf("--min_intensity > 0 requires usable intensity columns")
    remove_step(is.na(intensity) | !is.finite(intensity) | intensity < min_intensity, "N_low_intensity_removed")
  } else counts$N_low_intensity_removed <- 0L

  if (isTRUE(opt_get(options, "autosome_only", FALSE))) {
    if (is.null(chrom)) stopf("--autosome_only requires FinalReport chromosome column: %s", opt_get(options, "colchrom", "Chr"))
    chr <- toupper(gsub("^CHR", "", trimws(as.character(chrom))))
    is_auto <- chr %in% as.character(1:22)
    remove_step(!is_auto, "N_non_autosomal_removed")
  } else counts$N_non_autosomal_removed <- 0L

  counts$N_afterQC <- sum(active)
  counts$CallRate_after_markerQC <- if (sum(active) > 0) 1 - mean(abgeno[active] == 3L) else NA_real_
  counts$Nhom_afterQC <- sum(active & abgeno %in% c(0L, 2L), na.rm = TRUE)

  min_callrate <- num_option(options, "min_callrate", 0, 0, 1)
  counts$sample_pass_callrate <- is.na(counts$CallRate_after_markerQC) || counts$CallRate_after_markerQC >= min_callrate
  qc <- as.data.frame(counts, stringsAsFactors = FALSE, check.names = FALSE)
  list(keep = active, qc = qc)
}

## ---------------- FinalReport helpers ----------------

fr_metadata <- function(file, sep = NULL) {
  con <- file(file, open = "rt")
  on.exit(close(con), add = TRUE)
  line_no <- 0L
  data_line <- NA_integer_
  header_line <- NULL
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    line_no <- line_no + 1L
    line <- sub("\r$", "", line)
    if (trimws(line) == "[Data]") {
      data_line <- line_no
      header_line <- readLines(con, n = 1L, warn = FALSE)
      if (!length(header_line)) stopf("FinalReport [Data] found but column header is missing: %s", file)
      header_line <- sub("\r$", "", header_line)
      break
    }
  }
  if (is.na(data_line)) stopf("[Data] not found in FinalReport: %s", file)
  real_sep <- normalize_sep(sep)
  if (is.null(real_sep)) real_sep <- guess_sep_from_line(header_line)
  if (is.null(real_sep)) stopf("Could not detect FinalReport delimiter. Pass --sep comma or --sep tab")
  header <- split_line(header_line, real_sep)
  list(data_line = data_line, sep = real_sep, header = header)
}

read_finalreport_table <- function(file, needed_cols, sep = NULL, numeric_cols = character()) {
  meta <- fr_metadata(file, sep)
  header <- meta$header
  miss <- setdiff(needed_cols, header)
  if (length(miss)) stopf("Missing FinalReport column(s): %s", paste(miss, collapse = ", "))
  cc <- rep("NULL", length(header))
  names(cc) <- header
  for (nm in needed_cols) cc[which(header == nm)[1L]] <- if (nm %in% numeric_cols) "numeric" else "character"
  read.table(file, header = TRUE, sep = meta$sep, skip = meta$data_line,
             quote = "", comment.char = "", fill = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, colClasses = cc)
}

list_samples_finalreport <- function(file, colsample = "Sample ID", sep = NULL, keep_file = NULL, keep_col = 1L, keep_sep = NULL,
                                     remove_file = NULL, remove_col = 1L, remove_sep = NULL) {
  if (identical(colsample, "NONE")) return(file)
  dat <- read_finalreport_table(file, colsample, sep)
  samples <- unique(as.character(dat[[colsample]]))
  samples <- samples[!is.na(samples) & nzchar(samples)]
  samples[apply_name_filters(samples, keep_file, keep_col, keep_sep, remove_file, remove_col, remove_sep)]
}

geno_ab_to_code <- function(a1, a2) {
  g <- paste0(a1, a2)
  out <- rep(3L, length(g))
  out[g == "AA"] <- 0L
  out[g == "AB"] <- 1L
  out[g == "BA"] <- 1L
  out[g == "BB"] <- 2L
  out[g == "--"] <- 3L
  out[is.na(g) | !nzchar(g)] <- 3L
  out
}

compute_rebaf_vec <- function(dat, mode) {
  mode <- as.character(mode)
  if (mode %in% c("", "original", "Original", "ORIGINAL", "none", "NONE")) return(NULL)
  valid <- c("xy_y", "xy_x", "raw_y", "raw_x")
  if (!mode %in% valid) stopf("Invalid --rebaf_mode %s. Use original | %s", mode, paste(valid, collapse = " | "))
  if (mode %in% c("xy_y", "xy_x")) {
    x <- suppressWarnings(as.numeric(dat[["X"]]))
    y <- suppressWarnings(as.numeric(dat[["Y"]]))
    den <- x + y
    baf <- ifelse(is.na(den) | den <= 0, NaN, if (mode == "xy_y") y / den else x / den)
  } else {
    xr <- suppressWarnings(as.numeric(dat[["X Raw"]]))
    yr <- suppressWarnings(as.numeric(dat[["Y Raw"]]))
    den <- xr + yr
    baf <- ifelse(is.na(den) | den <= 0, NaN, if (mode == "raw_y") yr / den else xr / den)
  }
  baf
}

rebaf_finalreport_file <- function(input_file, output_file, mode = "xy_y", sep = NULL) {
  valid <- c("xy_y", "xy_x", "raw_y", "raw_x")
  if (!mode %in% valid) stopf("Invalid rebaf mode: %s", mode)
  meta <- fr_metadata(input_file, sep)
  header <- meta$header
  req <- c("B Allele Freq")
  if (mode %in% c("xy_y", "xy_x")) req <- c(req, "X", "Y") else req <- c(req, "X Raw", "Y Raw")
  miss <- setdiff(req, header)
  if (length(miss)) stopf("Missing FinalReport column(s) for rebaf: %s", paste(miss, collapse = ", "))
  idx_baf <- match("B Allele Freq", header)
  idx_x <- match("X", header); idx_y <- match("Y", header)
  idx_xr <- match("X Raw", header); idx_yr <- match("Y Raw", header)

  in_con <- file(input_file, open = "rt")
  out_con <- file(output_file, open = "wt")
  on.exit(close(in_con), add = TRUE)
  on.exit(close(out_con), add = TRUE)
  line_no <- 0L
  n_data <- 0L
  while (TRUE) {
    line <- readLines(in_con, n = 1L, warn = FALSE)
    if (!length(line)) break
    line_no <- line_no + 1L
    line <- sub("\r$", "", line)
    if (line_no <= meta$data_line + 1L) {
      writeLines(line, out_con, sep = "\n")
      next
    }
    a <- split_line(line, meta$sep)
    if (length(a) < length(header)) length(a) <- length(header)
    if (mode %in% c("xy_y", "xy_x")) {
      x <- suppressWarnings(as.numeric(a[[idx_x]])); y <- suppressWarnings(as.numeric(a[[idx_y]]))
      den <- x + y
      baf <- if (is.na(den) || den <= 0) NaN else if (mode == "xy_y") y / den else x / den
    } else {
      x <- suppressWarnings(as.numeric(a[[idx_xr]])); y <- suppressWarnings(as.numeric(a[[idx_yr]]))
      den <- x + y
      baf <- if (is.na(den) || den <= 0) NaN else if (mode == "raw_y") y / den else x / den
    }
    a[[idx_baf]] <- if (is.nan(baf) || is.na(baf)) "NaN" else sprintf("%.4f", baf)
    writeLines(paste(a, collapse = meta$sep), out_con, sep = "\n")
    n_data <- n_data + 1L
  }
  list(input = input_file, output = output_file, mode = mode, n_data = n_data)
}

read_finalreport_sample_data <- function(file, sample, options, maf_vec) {
  colsample <- opt_get(options, "colsample", "Sample ID")
  colmarker <- opt_get(options, "colmarker", "SNP Name")
  colbaf <- opt_get(options, "colbaf", "B Allele Freq")
  colab1 <- opt_get(options, "colab1", "Allele1 - AB")
  colab2 <- opt_get(options, "colab2", "Allele2 - AB")
  mode <- opt_get(options, "rebaf_mode", "original")
  sep <- opt_get(options, "sep", NULL)

  meta <- fr_metadata(file, sep)
  header <- meta$header
  effective_options <- options

  needed <- unique(c(colmarker, colbaf, colab1, colab2,
                     if (!identical(colsample, "NONE")) colsample else character()))
  numeric_cols <- character()

  if (mode %in% c("xy_y", "xy_x")) {
    needed <- unique(c(needed, "X", "Y"))
    numeric_cols <- c(numeric_cols, "X", "Y")
  }
  if (mode %in% c("raw_y", "raw_x")) {
    needed <- unique(c(needed, "X Raw", "Y Raw"))
    numeric_cols <- c(numeric_cols, "X Raw", "Y Raw")
  }
  if (mode %in% c("", "original", "Original", "ORIGINAL", "none", "NONE"))
    numeric_cols <- c(numeric_cols, colbaf)

  min_gencall <- num_option(options, "min_gencall", 0, 0, Inf)
  min_gcscore <- num_option(options, "min_gcscore", 0, 0, Inf)
  min_intensity <- num_option(options, "min_intensity", 0, 0, Inf)
  autosome_only <- isTRUE(opt_get(options, "autosome_only", FALSE))

  gencall_col <- opt_get(options, "colgencall", "GenCall Score")
  gcscore_col <- opt_get(options, "colgcscore", "GC Score")
  chrom_col <- opt_get(options, "colchrom", "Chr")

  gencall_available <- min_gencall <= 0 || gencall_col %in% header
  gcscore_available <- min_gcscore <= 0 || gcscore_col %in% header
  intensity_available <- min_intensity <= 0 || all(c("X", "Y") %in% header)
  chrom_available <- !autosome_only || chrom_col %in% header

  if (min_gencall > 0 && !gencall_available) {
    emit_qc_warning(options, paste0("missing_gencall:", gencall_col),
                    sprintf("FinalReport column '%s' was not found; GenCall QC was skipped.", gencall_col))
    effective_options$min_gencall <- 0
  } else if (min_gencall > 0) {
    needed <- unique(c(needed, gencall_col))
    numeric_cols <- c(numeric_cols, gencall_col)
  }

  if (min_gcscore > 0 && !gcscore_available) {
    emit_qc_warning(options, paste0("missing_gcscore:", gcscore_col),
                    sprintf("FinalReport column '%s' was not found; GC Score QC was skipped.", gcscore_col))
    effective_options$min_gcscore <- 0
  } else if (min_gcscore > 0) {
    needed <- unique(c(needed, gcscore_col))
    numeric_cols <- c(numeric_cols, gcscore_col)
  }

  if (min_intensity > 0 && !intensity_available) {
    missing_intensity <- setdiff(c("X", "Y"), header)
    emit_qc_warning(options, "missing_intensity_columns",
                    sprintf("FinalReport intensity column(s) %s were not found; intensity QC was skipped.",
                            paste(sprintf("'%s'", missing_intensity), collapse = ", ")))
    effective_options$min_intensity <- 0
  } else if (min_intensity > 0) {
    needed <- unique(c(needed, "X", "Y"))
    numeric_cols <- c(numeric_cols, "X", "Y")
  }

  if (autosome_only && !chrom_available) {
    emit_qc_warning(options, paste0("missing_chrom:", chrom_col),
                    sprintf("FinalReport chromosome column '%s' was not found; autosome filtering was skipped.", chrom_col))
    effective_options$autosome_only <- FALSE
  } else if (autosome_only) {
    needed <- unique(c(needed, chrom_col))
  }

  dat <- read_finalreport_table(file, needed, sep, numeric_cols = unique(numeric_cols))
  if (!identical(colsample, "NONE") && opt_has_value(sample))
    dat <- dat[as.character(dat[[colsample]]) == as.character(sample), , drop = FALSE]
  if (!nrow(dat)) {
    out <- data.frame(sample = character(), marker = character(), BAF = numeric(), ABGENO = integer(), MAF = numeric())
    attr(out, "qc") <- data.frame(sample = sample, N_total = 0L, N_afterQC = 0L,
                                  GenCall_QC = if (min_gencall > 0 && !gencall_available) "SKIPPED" else if (min_gencall > 0) "APPLIED" else "DISABLED",
                                  GCScore_QC = if (min_gcscore > 0 && !gcscore_available) "SKIPPED" else if (min_gcscore > 0) "APPLIED" else "DISABLED",
                                  Intensity_QC = if (min_intensity > 0 && !intensity_available) "SKIPPED" else if (min_intensity > 0) "APPLIED" else "DISABLED",
                                  Autosome_QC = if (autosome_only && !chrom_available) "SKIPPED" else if (autosome_only) "APPLIED" else "DISABLED",
                                  stringsAsFactors = FALSE)
    return(out)
  }

  marker <- as.character(dat[[colmarker]])
  baf_re <- compute_rebaf_vec(dat, mode)
  baf <- if (is.null(baf_re)) suppressWarnings(as.numeric(dat[[colbaf]])) else baf_re
  allele1 <- as.character(dat[[colab1]])
  allele2 <- as.character(dat[[colab2]])
  cnv <- grepl("CNV", toupper(allele1), fixed = TRUE) | grepl("CNV", toupper(allele2), fixed = TRUE)
  cnv[is.na(cnv)] <- FALSE
  abgeno <- geno_ab_to_code(allele1, allele2)
  maf <- maf_vec[marker]

  gencall <- if (min_gencall > 0 && gencall_available)
    suppressWarnings(as.numeric(dat[[gencall_col]])) else NULL
  gcscore <- if (min_gcscore > 0 && gcscore_available)
    suppressWarnings(as.numeric(dat[[gcscore_col]])) else NULL
  intensity <- if (min_intensity > 0 && intensity_available)
    suppressWarnings(as.numeric(dat[["X"]]) + as.numeric(dat[["Y"]])) else NULL
  chrom <- if (autosome_only && chrom_available) dat[[chrom_col]] else NULL

  q <- apply_snp_qc(sample, marker, baf, abgeno, maf, effective_options,
                    gencall, gcscore, intensity, chrom, cnv)
  q$qc$GenCall_QC <- if (min_gencall > 0 && !gencall_available) "SKIPPED" else if (min_gencall > 0) "APPLIED" else "DISABLED"
  q$qc$GCScore_QC <- if (min_gcscore > 0 && !gcscore_available) "SKIPPED" else if (min_gcscore > 0) "APPLIED" else "DISABLED"
  q$qc$Intensity_QC <- if (min_intensity > 0 && !intensity_available) "SKIPPED" else if (min_intensity > 0) "APPLIED" else "DISABLED"
  q$qc$Autosome_QC <- if (autosome_only && !chrom_available) "SKIPPED" else if (autosome_only) "APPLIED" else "DISABLED"

  list_keep <- apply_name_filters(marker,
                                  opt_get(options, "extract", NULL), opt_get(options, "extractcol", 1), opt_get(options, "extractsep", NULL),
                                  opt_get(options, "exclude", NULL), opt_get(options, "excludecol", 1), opt_get(options, "excludesep", NULL))
  q$qc$N_list_filter_removed <- sum(q$keep & !list_keep)
  keep <- q$keep & list_keep
  q$qc$N_afterQC <- sum(keep)
  q$qc$CallRate_after_markerQC <- if (sum(keep) > 0) 1 - mean(abgeno[keep] == 3L) else NA_real_
  q$qc$Nhom_afterQC <- sum(keep & abgeno %in% c(0L, 2L), na.rm = TRUE)
  min_callrate <- num_option(options, "min_callrate", 0, 0, 1)
  q$qc$sample_pass_callrate <- is.na(q$qc$CallRate_after_markerQC) || q$qc$CallRate_after_markerQC >= min_callrate

  smp <- if (identical(colsample, "NONE")) rep(sample %||% file, length(marker)) else as.character(dat[[colsample]])
  out <- data.frame(sample = smp[keep], marker = marker[keep], BAF = baf[keep],
                    ABGENO = abgeno[keep], MAF = as.numeric(maf[keep]), stringsAsFactors = FALSE)
  attr(out, "qc") <- q$qc
  out
}

## ---------------- VCF helpers ----------------

vcf_header <- function(file) {
  con <- file(file, open = "rt")
  on.exit(close(con), add = TRUE)
  line_no <- 0L
  hdr <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    line_no <- line_no + 1L
    line <- sub("\r$", "", line)
    if (startsWith(line, "#CHROM")) {
      hdr <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      break
    }
  }
  if (is.null(hdr)) stopf("VCF header line #CHROM not found: %s", file)
  list(header_line = line_no, header = hdr, samples = if (length(hdr) > 9L) hdr[10:length(hdr)] else character())
}

list_samples_vcf <- function(file, keep_file = NULL, keep_col = 1L, keep_sep = NULL,
                             remove_file = NULL, remove_col = 1L, remove_sep = NULL) {
  samples <- vcf_header(file)$samples
  samples[apply_name_filters(samples, keep_file, keep_col, keep_sep, remove_file, remove_col, remove_sep)]
}

parse_gt_to_ab <- function(gt) {
  gt <- sub(":.*$", "", gt)
  gt[is.na(gt)] <- "."
  out <- rep(3L, length(gt))
  gt0 <- gsub("\\|", "/", gt)
  out[gt0 %in% c("0/0", "0")] <- 0L
  out[gt0 %in% c("0/1", "1/0", "0/1/.", "1/0/.")] <- 1L
  out[gt0 %in% c("1/1", "1")] <- 2L
  out[grepl("\\.", gt0)] <- 3L
  out
}

extract_format_value <- function(fmt, sample_field, key) {
  keys <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
  idx <- match(key, keys)
  if (is.na(idx)) return(NA_character_)
  vals <- strsplit(sample_field, ":", fixed = TRUE)[[1L]]
  if (length(vals) < idx) return(NA_character_)
  vals[[idx]]
}

vcf_baf_one <- function(fmt, sample_field, baf_key = "BAF", ad_key = "AD") {
  keys <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
  vals <- strsplit(sample_field, ":", fixed = TRUE)[[1L]]
  idx <- match(baf_key, keys)
  if (!is.na(idx) && length(vals) >= idx) {
    z <- vals[[idx]]
    if (!is.na(z) && nzchar(z) && z != ".") {
      zz <- strsplit(z, ",", fixed = TRUE)[[1L]][[1L]]
      return(suppressWarnings(as.numeric(zz)))
    }
  }
  idx_ad <- match(ad_key, keys)
  if (!is.na(idx_ad) && length(vals) >= idx_ad) {
    ad <- suppressWarnings(as.numeric(strsplit(vals[[idx_ad]], ",", fixed = TRUE)[[1L]]))
    if (length(ad) >= 2L && !any(is.na(ad[1:2]))) {
      den <- ad[[1L]] + ad[[2L]]
      if (den > 0) return(ad[[2L]] / den)
    }
  }
  NA_real_
}

read_vcf_sample_data <- function(file, sample, options, maf_vec) {
  vh <- vcf_header(file)
  hdr <- vh$header
  samples <- vh$samples
  if (!sample %in% samples) stopf("Sample not found in VCF: %s", sample)
  sample_col <- match(sample, hdr)
  cc <- rep("NULL", length(hdr))
  names(cc) <- hdr
  keep_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "FORMAT", sample)
  for (nm in keep_cols) if (nm %in% names(cc)) cc[[nm]] <- "character"
  dat <- read.table(file, header = TRUE, sep = "\t", skip = vh$header_line - 1L,
                    quote = "", comment.char = "", fill = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, colClasses = cc)
  if (!nrow(dat)) return(data.frame(sample = character(), marker = character(), BAF = numeric(), ABGENO = integer(), MAF = numeric()))
  marker <- as.character(dat[["ID"]])
  bad_id <- is.na(marker) | !nzchar(marker) | marker == "."
  marker[bad_id] <- paste(dat[["#CHROM"]][bad_id], dat[["POS"]][bad_id], sep = ":")
  sf <- as.character(dat[[sample]])
  fmt <- as.character(dat[["FORMAT"]])
  baf_key <- opt_get(options, "vcf_baf_key", "BAF")
  ad_key <- opt_get(options, "vcf_ad_key", "AD")
  baf <- mapply(vcf_baf_one, fmt, sf, MoreArgs = list(baf_key = baf_key, ad_key = ad_key), USE.NAMES = FALSE)
  gt <- mapply(extract_format_value, fmt, sf, MoreArgs = list(key = opt_get(options, "vcf_gt_key", "GT")), USE.NAMES = FALSE)
  abgeno <- parse_gt_to_ab(gt)
  maf <- maf_vec[marker]
  vcf_qc_options <- options
  vcf_qc_options$min_gencall <- 0
  vcf_qc_options$min_gcscore <- 0
  vcf_qc_options$min_intensity <- 0
  vcf_qc_options$remove_input_cnv <- FALSE
  q <- apply_snp_qc(sample, marker, baf, abgeno, maf, vcf_qc_options, chrom = dat[["#CHROM"]])
  list_keep <- apply_name_filters(marker,
                                  opt_get(options, "extract", NULL), opt_get(options, "extractcol", 1), opt_get(options, "extractsep", NULL),
                                  opt_get(options, "exclude", NULL), opt_get(options, "excludecol", 1), opt_get(options, "excludesep", NULL))
  q$qc$N_list_filter_removed <- sum(q$keep & !list_keep)
  keep <- q$keep & list_keep
  q$qc$N_afterQC <- sum(keep)
  q$qc$CallRate_after_markerQC <- if (sum(keep) > 0) 1 - mean(abgeno[keep] == 3L) else NA_real_
  q$qc$Nhom_afterQC <- sum(keep & abgeno %in% c(0L, 2L), na.rm = TRUE)
  min_callrate <- num_option(options, "min_callrate", 0, 0, 1)
  q$qc$sample_pass_callrate <- is.na(q$qc$CallRate_after_markerQC) || q$qc$CallRate_after_markerQC >= min_callrate
  out <- data.frame(sample = rep(sample, sum(keep)), marker = marker[keep], BAF = as.numeric(baf[keep]), ABGENO = abgeno[keep], MAF = as.numeric(maf[keep]), stringsAsFactors = FALSE)
  attr(out, "qc") <- q$qc
  out
}

infer_format <- function(file, format = "auto") {
  fmt <- tolower(as.character(format))
  if (fmt %in% c("finalreport", "fr", "illumina")) return("finalreport")
  if (fmt %in% c("vcf", "bcf")) return("vcf")
  lines <- read_nonempty_head(file, 50L)
  if (any(grepl("^##fileformat=VCF", lines)) || any(grepl("^#CHROM", lines))) return("vcf")
  if (any(trimws(lines) == "[Data]")) return("finalreport")
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("vcf")) return("vcf")
  "finalreport"
}

list_input_samples <- function(file, format, options) {
  fmt <- infer_format(file, format)
  keep_file <- opt_get(options, "keep", NULL); keep_col <- opt_get(options, "keepcol", 1); keep_sep <- opt_get(options, "keepsep", NULL)
  remove_file <- opt_get(options, "remove", NULL); remove_col <- opt_get(options, "removecol", 1); remove_sep <- opt_get(options, "removesep", NULL)
  if (fmt == "vcf") {
    list_samples_vcf(file, keep_file, keep_col, keep_sep, remove_file, remove_col, remove_sep)
  } else {
    list_samples_finalreport(file, opt_get(options, "colsample", "Sample ID"), opt_get(options, "sep", NULL),
                             keep_file, keep_col, keep_sep, remove_file, remove_col, remove_sep)
  }
}

read_sample_data <- function(file, format, sample, options, maf_vec) {
  fmt <- infer_format(file, format)
  if (fmt == "vcf") read_vcf_sample_data(file, sample, options, maf_vec) else read_finalreport_sample_data(file, sample, options, maf_vec)
}

## ---------------- analysis ----------------

testsamplecontamination <- function(baf, abgeno, maf, subset = NULL, ...) {
  stopifnot(length(baf) == length(maf), length(baf) == length(abgeno))
  if (!length(baf)) stopf("No usable markers after filtering/MAF matching")
  maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  amaf <- ifelse(abgeno == 2L, -maf, maf)
  callrate <- 1 - mean(abgeno == 3L)
  subs <- abgeno == 2L | abgeno == 0L
  if (!is.null(subset)) subs <- subs & subset
  if (sum(subs) < 2L) stopf("Too few homozygous markers for regression: %d", sum(subs))
  genocat <- factor(abgeno, levels = c(0L, 2L))
  tb <- table(genocat[subs])
  fit <- if (length(tb) == 2L && all(tb > 0L)) {
    lm(baf ~ amaf + genocat, subset = subs, ...)
  } else {
    lm(baf ~ amaf, subset = subs, ...)
  }
  sm <- coefficients(summary(fit))
  if (!"amaf" %in% rownames(sm)) stopf("Regression did not return amaf coefficient")
  a <- c(sm["amaf", ], callrate, nrow(fit$model))
  names(a) <- c("estimate", "stderr", "tval", "pval", "callrate", "Nhom")
  a
}

estimate_one <- function(file, format, sample, options, maf_vec) {
  dd <- read_sample_data(file, format, sample, options, maf_vec)
  qc <- attr(dd, "qc")
  min_callrate <- num_option(options, "min_callrate", 0, 0, 1)
  if (!is.null(qc) && min_callrate > 0 && !isTRUE(qc$sample_pass_callrate[[1L]]))
    stopf("Sample %s failed --min_callrate %.4f (observed %.4f)", sample, min_callrate, qc$CallRate_after_markerQC[[1L]])
  reg <- testsamplecontamination(dd$BAF, dd$ABGENO, dd$MAF)
  ans <- data.frame(sample = sample, as.data.frame(t(reg), check.names = FALSE), stringsAsFactors = FALSE)
  attr(ans, "qc") <- qc
  ans
}

plot_one <- function(file, format, sample, options, maf_vec, outprefix) {
  dd <- read_sample_data(file, format, sample, options, maf_vec)
  png(paste(outprefix, "png", sep = "."))
  on.exit(dev.off(), add = TRUE)
  groups <- factor(dd$ABGENO, levels = 0:3, labels = c("AA", "AB", "BB", "--"))
  plot(dd$MAF, dd$BAF, xlab = "MAF", ylab = "BAF", main = sample, pch = 16, cex = 0.45)
  lev <- levels(groups)
  for (i in seq_along(lev)) {
    idx <- groups == lev[[i]]
    if (any(idx)) points(dd$MAF[idx], dd$BAF[idx], pch = i, cex = 0.45)
  }
  legend("topright", legend = lev, pch = seq_along(lev), bty = "n")
  invisible(TRUE)
}

## ---------------- pipeline ----------------

pipeline_usage <- function() {
  cat("Usage:\n")
  cat("  Rscript run_bafRegress_allinone_baseR_v2_FR_VCF.R --final FILE --freqfile FILE [options]\n\n")
  cat("Input options:\n")
  cat("  --final FILE             Input FinalReport or VCF. Alias: --input\n")
  cat("  --format auto|finalreport|vcf   Input type [default: auto]\n")
  cat("  --sep auto|comma|tab|CHAR       FinalReport delimiter [default: auto]\n")
  cat("  --freqsep auto|comma|tab|CHAR   Frequency file delimiter [default: whitespace]\n")
  cat("  --freqcols 1,2          Marker/frequency columns in freqfile [default: 1,2]\n")
  cat("  --min_maf NUM           Optional minimum MAF after folding to <=0.5 [default: 0; official behavior]\n\n")
  cat("Conservative QC options:\n")
  cat("  --remove_duplicates BOOL  Remove repeated marker IDs [default: TRUE]\n")
  cat("  --valid_baf_only BOOL     Require finite BAF in [0,1] [default: TRUE]\n")
  cat("  --min_gencall NUM         Minimum GenCall Score [default: 0.15]\n")
  cat("  --min_gcscore NUM         Minimum GC Score [default: 0.15]\n")
  cat("  --min_intensity NUM       Minimum normalized X+Y [default: 1]\n")
  cat("  --min_callrate NUM        Minimum post-marker-QC call rate [default: 0.90]\n")
  cat("  --autosome_only BOOL      Keep chromosomes 1-22 only [default: TRUE]\n")
  cat("  --remove_input_cnv BOOL   Remove records explicitly marked CNV in allele fields [default: TRUE]\n")
  cat("  --colgencall NAME         GenCall column [default: GenCall Score]\n")
  cat("  --colgcscore NAME         GC Score column [default: GC Score]\n")
  cat("  --colchrom NAME           Chromosome column [default: Chr]\n")
  cat("  Missing optional QC columns are skipped with warnings written to screen and logs/qc_warnings.log.txt.\n\n")
  cat("BAF options for FinalReport:\n")
  cat("  --rebaf_mode original|xy_y|xy_x|raw_y|raw_x [default: original]\n")
  cat("  --rebaf_out FILE        Optional output path if rebaf_mode is not original\n\n")
  cat("VCF options:\n")
  cat("  --vcf_baf_key BAF       FORMAT key for BAF [default: BAF]\n")
  cat("  --vcf_gt_key GT         FORMAT key for genotype [default: GT]\n")
  cat("  --vcf_ad_key AD         FORMAT key used to compute BAF if BAF missing [default: AD]\n\n")
  cat("Sample options:\n")
  cat("  --sample NAME           Analyze one sample only; default analyzes all samples\n")
  cat("  --keep FILE / --remove FILE  Sample include/exclude lists\n\n")
  cat("Output options:\n")
  cat("  --outdir DIR            Base output directory [default: baf_out]\n")
  cat("  --runid ID              Custom run directory name\n")
  cat("  --no_plot               Skip plots\n")
  cat("  --no_makeR              Accepted for compatibility; no input-prep is run in this all-in-one mode\n")
}

write_log <- function(path, expr_fun) {
  out <- tryCatch({
    res <- expr_fun()
    list(ok = TRUE, result = res, msg = NULL)
  }, error = function(e) list(ok = FALSE, result = NULL, msg = conditionMessage(e)))
  if (out$ok) {
    if (is.data.frame(out$result)) {
      utils::capture.output(write.table(out$result, quote = FALSE, row.names = FALSE, sep = "\t"), file = path)
    } else {
      writeLines(as.character(out$result), path)
    }
  } else {
    writeLines(out$msg, path)
  }
  out
}

run_pipeline <- function(callargs) {
  defaults <- list(
    final = "", input = "", format = "auto", sep = NULL,
    freqfile = "", freqcols = "1,2", freqsep = NULL, min_maf = 0,
    remove_duplicates = TRUE, valid_baf_only = TRUE, min_gencall = 0.15, min_gcscore = 0.15,
    min_intensity = 1, min_callrate = 0.90, autosome_only = TRUE, remove_input_cnv = TRUE,
    colgencall = "GenCall Score", colgcscore = "GC Score", colchrom = "Chr",
    outdir = "baf_out", runid = "", no_plot = FALSE, no_makeR = FALSE,
    sample = "", colsample = "Sample ID", colmarker = "SNP Name",
    colbaf = "B Allele Freq", colab1 = "Allele1 - AB", colab2 = "Allele2 - AB",
    rebaf_mode = "original", rebaf_out = "",
    vcf_baf_key = "BAF", vcf_gt_key = "GT", vcf_ad_key = "AD",
    extract = NULL, extractcol = 1, extractsep = NULL,
    exclude = NULL, excludecol = 1, excludesep = NULL,
    keep = NULL, keepcol = 1, keepsep = NULL,
    remove = NULL, removecol = 1, removesep = NULL
  )
  parsed <- parse_cli(callargs, defaults, bool_flags = c("no_plot", "no_makeR", "remove_duplicates", "valid_baf_only", "autosome_only", "remove_input_cnv", "help"))
  opts <- parsed$options
  if (isTRUE(opt_get(opts, "help", FALSE))) return(pipeline_usage())

  input <- opt_get(opts, "final", "")
  if (!nzchar(input)) input <- opt_get(opts, "input", "")
  if (!nzchar(input)) stopf("ERROR: --final or --input is required")
  require_file(input, "input")
  require_file(opt_get(opts, "freqfile", ""), "freqfile")

  input_format <- infer_format(input, opt_get(opts, "format", "auto"))
  run_plot <- !isTRUE(opt_get(opts, "no_plot", FALSE))
  outdir_base <- sub("/+$", "", opt_get(opts, "outdir", "baf_out"))
  run_id <- opt_get(opts, "runid", "")
  if (!nzchar(run_id)) run_id <- paste0("run_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  outdir <- file.path(outdir_base, run_id)
  dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_base, recursive = TRUE, showWarnings = FALSE)
  opts$warning_log <- file.path(outdir, "logs", "qc_warnings.log.txt")
  latest <- file.path(outdir_base, "latest")
  suppressWarnings(unlink(latest, recursive = TRUE, force = TRUE))
  if (!isTRUE(suppressWarnings(file.symlink(run_id, latest)))) writeLines(run_id, paste0(latest, ".txt"))

  cat(">>> Using run directory: ", outdir, "\n", sep = "")
  cat("    (base outdir: ", outdir_base, ")\n", sep = "")
  cat("    input format: ", input_format, "\n", sep = "")

  analysis_file <- input
  rebaf_mode <- opt_get(opts, "rebaf_mode", "original")
  if (input_format == "vcf" && !tolower(rebaf_mode) %in% c("", "original", "none")) {
    stopf("--rebaf_mode is only for FinalReport. For VCF, use FORMAT BAF or AD fields.")
  }
  if (input_format == "finalreport" && !tolower(rebaf_mode) %in% c("", "original", "none")) {
    rebaf_out <- opt_get(opts, "rebaf_out", "")
    if (!nzchar(rebaf_out)) {
      base <- sub("\\.[^.]*$", "", basename(input))
      rebaf_out <- file.path(outdir, "tables", paste0(base, ".", rebaf_mode, ".rebaf.csv"))
    }
    cat("[pre] Recomputing B Allele Freq using mode ", rebaf_mode, " ...\n", sep = "")
    rlog <- file.path(outdir, "logs", "rebaf.log.txt")
    rb <- tryCatch(rebaf_finalreport_file(input, rebaf_out, rebaf_mode, opt_get(opts, "sep", NULL)), error = function(e) e)
    if (inherits(rb, "error")) {
      writeLines(conditionMessage(rb), rlog)
      stopf("ERROR: BAF recomputation failed. See: %s", rlog)
    }
    writeLines(c(paste0("input\t", input), paste0("output\t", rebaf_out), paste0("mode\t", rebaf_mode), paste0("data_rows\t", rb$n_data)), rlog)
    analysis_file <- rebaf_out
    opts$rebaf_mode <- "original"
    cat("  OK. Recomputed FinalReport: ", analysis_file, "\n", sep = "")
  }

  run_info_list <- list(
    run_id = run_id,
    run_dir = outdir,
    started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    input = input,
    input_format = input_format,
    analysis_file = analysis_file,
    rebaf_mode = rebaf_mode,
    freqfile = opt_get(opts, "freqfile", ""),
    freqcols = opt_get(opts, "freqcols", "1,2"),
    sep = opt_get(opts, "sep", "auto"),
    freqsep = opt_get(opts, "freqsep", "auto"),
    min_maf = opt_get(opts, "min_maf", 0),
    remove_duplicates = opt_get(opts, "remove_duplicates", TRUE),
    valid_baf_only = opt_get(opts, "valid_baf_only", TRUE),
    min_gencall = opt_get(opts, "min_gencall", 0.15),
    min_gcscore = opt_get(opts, "min_gcscore", 0.15),
    min_intensity = opt_get(opts, "min_intensity", 1),
    min_callrate = opt_get(opts, "min_callrate", 0.90),
    autosome_only = opt_get(opts, "autosome_only", TRUE),
    remove_input_cnv = opt_get(opts, "remove_input_cnv", TRUE),
    run_plot = as.integer(run_plot)
  )
  run_info <- data.frame(
    key = names(run_info_list),
    value = unname(vapply(run_info_list, function(x) {
      if (is.null(x) || !length(x) || is.na(x[[1L]])) "" else as.character(x[[1L]])
    }, character(1))),
    stringsAsFactors = FALSE
  )
  write.table(run_info, file.path(outdir, "run_info.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  writeLines(paste(c("Rscript", shQuote(commandArgs(FALSE)[1]), shQuote(commandArgs(trailingOnly = TRUE))), collapse = " "), file.path(outdir, "cmd.sh"))
  Sys.chmod(file.path(outdir, "cmd.sh"), mode = "0755")

  cat("[1/4] Loading frequency file ...\n")
  maf_vec <- load_frequency(opt_get(opts, "freqfile", ""), opt_get(opts, "freqcols", "1,2"), opt_get(opts, "freqsep", NULL), opt_get(opts, "min_maf", 0))
  cat("  Markers with MAF: ", length(maf_vec), "\n", sep = "")

  cat("[2/4] Extracting samples ...\n")
  one_sample <- opt_get(opts, "sample", "")
  samples <- if (nzchar(one_sample)) one_sample else list_input_samples(analysis_file, input_format, opts)
  samples <- samples[!is.na(samples) & nzchar(samples)]
  samples_file <- file.path(outdir, "tables", "samples.txt")
  writeLines(samples, samples_file)
  if (!length(samples)) stopf("ERROR: No samples extracted. Check input format / sample column / VCF header.")
  cat("  Samples:\n")
  cat(paste0("    - ", samples, "\n"), sep = "")

  cat("[3/4] Running estimate for each sample ...\n")
  summary_tsv <- file.path(outdir, "tables", "estimate_summary.tsv")
  qc_tsv <- file.path(outdir, "tables", "qc_summary.tsv")
  header_written <- FALSE
  qc_header_written <- FALSE
  for (s in samples) {
    safe <- safe_name(s)
    elog <- file.path(outdir, "logs", paste0(safe, ".estimate.log.txt"))
    cat("  - estimate ", s, "\n", sep = "")
    ans <- write_log(elog, function() estimate_one(analysis_file, input_format, s, opts, maf_vec))
    if (!ans$ok) {
      cat("    ERROR: estimate failed for ", s, ". See ", elog, "\n", sep = "")
    } else {
      write.table(ans$result, summary_tsv, quote = FALSE, row.names = FALSE, sep = "\t",
                  append = header_written, col.names = !header_written)
      header_written <- TRUE
      qc <- attr(ans$result, "qc")
      if (!is.null(qc)) {
        write.table(qc, qc_tsv, quote = FALSE, row.names = FALSE, sep = "\t",
                    append = qc_header_written, col.names = !qc_header_written)
        qc_header_written <- TRUE
      }
    }
    if (run_plot && ans$ok) {
      plog <- file.path(outdir, "logs", paste0(safe, ".plot.log.txt"))
      outprefix <- file.path(outdir, "plots", safe)
      pok <- tryCatch({ plot_one(analysis_file, input_format, s, opts, maf_vec, outprefix); TRUE }, error = function(e) { writeLines(conditionMessage(e), plog); FALSE })
      if (!pok) cat("    WARNING: plot failed for ", s, ". See ", plog, "\n", sep = "")
    }
  }
  if (!file.exists(summary_tsv)) writeLines("sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom", summary_tsv)
  if (!file.exists(qc_tsv)) writeLines("sample\tN_total\tN_afterQC", qc_tsv)

  cat("[4/4] Done.\n")
  cat("  - Run dir:  ", outdir, "\n", sep = "")
  cat("  - Latest:   ", latest, "\n", sep = "")
  cat("  - Samples:  ", samples_file, "\n", sep = "")
  cat("  - Summary:  ", summary_tsv, "\n", sep = "")
  cat("  - QC:       ", qc_tsv, "\n", sep = "")
  cat("  - Plots:    ", file.path(outdir, "plots"), "/\n", sep = "")
  cat("  - Logs:     ", file.path(outdir, "logs"), "/\n", sep = "")
  invisible(outdir)
}

## ---------------- compatibility actions ----------------

run_estimate_action <- function(callargs) {
  defaults <- list(format = "auto", sep = NULL, freqfile = "", freqcols = "1,2", freqsep = NULL, min_maf = 0,
                   remove_duplicates = TRUE, valid_baf_only = TRUE, min_gencall = 0.15, min_gcscore = 0.15,
                   min_intensity = 1, min_callrate = 0.90, autosome_only = TRUE, remove_input_cnv = TRUE,
                   colgencall = "GenCall Score", colgcscore = "GC Score", colchrom = "Chr",
                   sample = "", colsample = "Sample ID", colmarker = "SNP Name", colbaf = "B Allele Freq",
                   colab1 = "Allele1 - AB", colab2 = "Allele2 - AB", rebaf_mode = "original",
                   vcf_baf_key = "BAF", vcf_gt_key = "GT", vcf_ad_key = "AD",
                   extract = NULL, extractcol = 1, extractsep = NULL, exclude = NULL, excludecol = 1, excludesep = NULL,
                   keep = NULL, keepcol = 1, keepsep = NULL, remove = NULL, removecol = 1, removesep = NULL)
  parsed <- parse_cli(callargs, defaults, bool_flags = c("remove_duplicates", "valid_baf_only", "autosome_only", "remove_input_cnv"))
  opts <- parsed$options
  if (length(parsed$args) < 1L) stopf("estimate requires input file")
  file <- parsed$args[[1L]]
  fmt <- infer_format(file, opt_get(opts, "format", "auto"))
  maf <- load_frequency(opt_get(opts, "freqfile", ""), opt_get(opts, "freqcols", "1,2"), opt_get(opts, "freqsep", NULL), opt_get(opts, "min_maf", 0))
  samples <- if (nzchar(opt_get(opts, "sample", ""))) opt_get(opts, "sample", "") else list_input_samples(file, fmt, opts)
  res <- do.call(rbind, lapply(samples, function(s) estimate_one(file, fmt, s, opts, maf)))
  write.table(res, stdout(), quote = FALSE, row.names = FALSE, sep = "\t")
}

run_listsamples_action <- function(callargs) {
  defaults <- list(format = "auto", sep = NULL, colsample = "Sample ID",
                   keep = NULL, keepcol = 1, keepsep = NULL, remove = NULL, removecol = 1, removesep = NULL)
  parsed <- parse_cli(callargs, defaults)
  opts <- parsed$options
  if (length(parsed$args) < 1L) stopf("listsamplesraw requires input file")
  file <- parsed$args[[1L]]
  fmt <- infer_format(file, opt_get(opts, "format", "auto"))
  samples <- list_input_samples(file, fmt, opts)
  cat("sample\n")
  cat(paste0(samples, "\n"), sep = "")
}

run_rebaf_action <- function(callargs) {
  defaults <- list(sep = NULL)
  parsed <- parse_cli(callargs, defaults)
  if (length(parsed$args) != 3L) stopf("Usage: rebaf input_finalreport output_finalreport MODE")
  ans <- rebaf_finalreport_file(parsed$args[[1L]], parsed$args[[2L]], parsed$args[[3L]], opt_get(parsed$options, "sep", NULL))
  cat("Done. Wrote: ", ans$output, "\n", sep = "")
}

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  actions <- c("estimate", "listsamplesraw", "rebaf")
  if (length(argv) && argv[[1L]] %in% actions) {
    action <- argv[[1L]]
    rest <- argv[-1L]
    switch(action,
           estimate = run_estimate_action(rest),
           listsamplesraw = run_listsamples_action(rest),
           rebaf = run_rebaf_action(rest))
  } else {
    run_pipeline(argv)
  }
}

main()
