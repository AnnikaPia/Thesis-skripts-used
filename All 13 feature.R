suppressPackageStartupMessages({
  library(DNAshapeR)
  library(Biostrings)
  library(data.table)
})

# Convert one DNAshapeR feature into a data.table
# rows = sequences, cols = positions
# If a feature is shorter than expected, pad with NA
# If longer than expected, trim to expected length
shape_to_dt <- function(x, shape_name, L, n_seq, is_step = NULL) {
  if (is.null(x)) return(NULL)
  
  # normalize to matrix
  if (is.list(x)) {
    lens <- unique(lengths(x))
    if (length(lens) != 1) {
      stop("Inconsistent lengths in feature: ", shape_name)
    }
    mat <- t(vapply(x, function(v) as.numeric(v), numeric(lens[1])))
  } else if (is.data.frame(x)) {
    mat <- as.matrix(x)
  } else if (is.matrix(x)) {
    mat <- x
  } else {
    stop("Unsupported type for feature: ", shape_name)
  }
  
  # orient so rows = sequences
  if (nrow(mat) != n_seq && ncol(mat) == n_seq) {
    mat <- t(mat)
  }
  
  if (nrow(mat) != n_seq) {
    stop(sprintf(
      "Unexpected number of rows for %s: got %d, expected %d sequences",
      shape_name, nrow(mat), n_seq
    ))
  }
  
  # expected length
  target_cols <- if (isTRUE(is_step)) (L - 1) else L
  observed_cols <- ncol(mat)
  
  # pad shorter features with NA
  if (observed_cols < target_cols) {
    pad_n <- target_cols - observed_cols
    pad <- matrix(NA_real_, nrow = nrow(mat), ncol = pad_n)
    mat <- cbind(mat, pad)
    message("Padded ", shape_name, " from ", observed_cols, " to ", target_cols, " columns with NA")
  }
  
  # trim longer features if needed
  if (ncol(mat) > target_cols) {
    mat <- mat[, seq_len(target_cols), drop = FALSE]
    message("Trimmed ", shape_name, " from ", observed_cols, " to ", target_cols, " columns")
  }
  
  colnames(mat) <- sprintf("%s_%02d", shape_name, seq_len(ncol(mat)))
  as.data.table(mat)
}

make_shape_features <- function(fasta_path,
                                out_csv = "shape_features.csv",
                                add_simple_seq_feats = TRUE) {
  
  # ---- check file ----
  if (!file.exists(fasta_path)) {
    stop("File not found: ", fasta_path)
  }
  
  # ---- read FASTA ----
  faset <- readDNAStringSet(fasta_path)
  seq_ids <- names(faset)
  if (is.null(seq_ids) || all(!nzchar(seq_ids))) {
    seq_ids <- paste0("seq", seq_along(faset))
  }
  
  seqs <- toupper(gsub("U", "T", as.character(faset)))
  
  # ---- checks ----
  if (length(unique(nchar(seqs))) != 1) {
    stop("All sequences must have the same length.")
  }
  
  if (any(grepl("[^ACGT]", seqs))) {
    bad <- which(grepl("[^ACGT]", seqs))
    print(data.frame(id = seq_ids[bad], seq = seqs[bad]))
    stop("Found non-ACGT characters. Fix the FASTA.")
  }
  
  L <- nchar(seqs[1])
  n_seq <- length(seqs)
  message("Sequences: ", n_seq, " | Length: ", L)
  
  # ---- explicitly request the 13 features ----
  wanted_features <- c(
    "MGW", "ProT", "HelT", "Roll",
    "Stretch", "Buckle", "Tilt", "Shear", "Opening",
    "Rise", "Shift", "Stagger", "Slide"
  )
  
  pred <- getShape(fasta_path, shapeType = wanted_features)
  available_features <- names(pred)
  message("Returned features: ", paste(available_features, collapse = ", "))
  
  # ---- which features are step parameters ----
  feature_is_step <- c(
    MGW = FALSE,
    ProT = FALSE,
    HelT = TRUE,
    Roll = TRUE,
    Stretch = FALSE,
    Buckle = FALSE,
    Tilt = FALSE,
    Shear = FALSE,
    Opening = FALSE,
    Rise = TRUE,
    Shift = TRUE,
    Stagger = FALSE,
    Slide = TRUE
  )
  
  # ---- base output table ----
  DT <- data.table(sequence_id = seq_ids, sequence = seqs)
  
  missing_features <- character(0)
  added_features <- character(0)
  
  # ---- add features one by one ----
  for (nm in wanted_features) {
    if (!nm %in% available_features || is.null(pred[[nm]])) {
      warning("Feature not returned by getShape(): ", nm)
      missing_features <- c(missing_features, nm)
      next
    }
    
    message("Processing feature: ", nm)
    
    feat_dt <- shape_to_dt(
      x = pred[[nm]],
      shape_name = nm,
      L = L,
      n_seq = n_seq,
      is_step = feature_is_step[[nm]]
    )
    
    DT <- cbind(DT, feat_dt)
    DT <- as.data.table(DT)
    added_features <- c(added_features, nm)
  }
  
  # ---- optional simple sequence feature ----
  if (add_simple_seq_feats) {
    GC_frac <- vapply(seqs, function(s) {
      nchar(gsub("[^GC]", "", s)) / nchar(s)
    }, numeric(1))
    DT$GC_frac <- GC_frac
  }
  
  # ---- write output ----
  fwrite(DT, out_csv)
  
  message("Added features: ", paste(added_features, collapse = ", "))
  if (length(missing_features) > 0) {
    message("Missing features: ", paste(missing_features, collapse = ", "))
  }
  
  message("Wrote: ", normalizePath(out_csv, mustWork = FALSE))
  invisible(DT)
}

# ---- RUN ----
input_fasta <- "panel20_with_flanks.fasta"
output_csv  <- "Up_library_features_unseen.csv"

result <- make_shape_features(
  fasta_path = input_fasta,
  out_csv = output_csv,
  add_simple_seq_feats = TRUE
)

