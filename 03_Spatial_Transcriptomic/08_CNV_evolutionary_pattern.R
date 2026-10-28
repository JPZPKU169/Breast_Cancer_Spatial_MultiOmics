# Assign three DCIS-to-IDC evolutionary patterns from SpatialInferCNV clone
# results:
#   1.  single-lineage invasion
#   2a. lineage-divergent invasion
#   2b. Multiclonal cross-compartment sharing
#   Indeterminate. insufficient reliable support for pattern assignment
#
# Compared with assign_CNV_clone.R, this version recognizes an IDC-private
# clone as further diversification of a shared clone when its CNV events are
# directionally nested and the required tree/distance structure is satisfied.
# It does not force a pattern when either DCIS or IDC lacks a reliably
# supported clone.
#
# Required input files in each sample directory:
#   0.3.2.metacell_clone_info_final_dynamic_tree_cut<analysis_suffix>.csv
#   where a v3 suffix is, for example:
#   _split0.1_cnvContinuous_treeNJ
#   Legacy clone files without an analysis suffix are also supported.
#   0.3.infercnv.observations.txt
#   0.3.infercnv.reference.txt
#   0.3.GeneOrderFile.txt
#
# Output files are written to output_dir:
#   CNV_clone_evolution_mode_summary.csv
#   CNV_clone_evolution_clone_metrics.csv
#   CNV_clone_evolution_parameters.csv
#   CNV_clone_evolution_trees.rds
#   CNV_clone_evolution_nestedness_metrics.csv

suppressPackageStartupMessages({
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required. Install it before running this script.")
  }
})

# ==============================================================================
# User settings
# ==============================================================================

infercnv_root <- "E:/Projects/BRCA/results/ST/infercnv_metacell"
output_dir <- "E:/Projects/BRCA/results/ST/infercnv_metacell/"

# NULL means: analyze every sample directory containing a final clone file.
sample_ids <- c(
  "P1","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13","P14","P15","P16","P26","P27","P30","P32","P33","P34","P35", 
  "P36","P37","P38","P39","P40_1","P42","P43_2","P44","P45","P47"
)

# Select the st_SpatialInferCNV_dynamic_tree_cut_v3.R result to read. These
# settings must match the arguments used to generate the clone assignments.
clone_result_dynamic_deep_split <- 0.1
clone_result_discretize_cnv_pattern <- FALSE
clone_result_tree_method <- "NJ"

if (length(clone_result_dynamic_deep_split) != 1 ||
    !is.numeric(clone_result_dynamic_deep_split) ||
    !is.finite(clone_result_dynamic_deep_split)) {
  stop("clone_result_dynamic_deep_split must be one finite numeric value.")
}
if (length(clone_result_discretize_cnv_pattern) != 1 ||
    !is.logical(clone_result_discretize_cnv_pattern) ||
    is.na(clone_result_discretize_cnv_pattern)) {
  stop("clone_result_discretize_cnv_pattern must be TRUE or FALSE.")
}
clone_result_tree_method <- toupper(as.character(clone_result_tree_method))
if (length(clone_result_tree_method) != 1 ||
    !(clone_result_tree_method %in% c("NJ", "MP"))) {
  stop("clone_result_tree_method must be 'NJ' or 'MP'.")
}

clone_result_split_label <- format(
  clone_result_dynamic_deep_split,
  trim = TRUE,
  scientific = FALSE
)
clone_result_cnv_label <- if (clone_result_discretize_cnv_pattern) {
  "cnv3State"
} else {
  "cnvContinuous"
}
clone_result_analysis_suffix <- paste0(
  "_split", clone_result_split_label,
  "_", clone_result_cnv_label,
  "_tree", clone_result_tree_method
)

clone_file_candidates <- c(
  paste0(
    "0.3.2.metacell_clone_info_final_dynamic_tree_cut",
    clone_result_analysis_suffix,
    ".csv"
  ),
  paste0(
    "Clone_info_full_dynamic_tree_cut",
    clone_result_analysis_suffix,
    ".csv"
  ),
  "0.3.2.metacell_clone_info_final_dynamic_tree_cut.csv",
  "Clone_info_full_dynamic_tree_cut.csv"
)

##### Set parameters #####
params <- list(
  min_pathology_metacells = 50,
  min_clone_pathology_count = 50,
  min_shared_fraction = 0.05,
  infercnv_zero_threshold = log(1.02, 2),
  cnv_retention_threshold = 0.95,
  cnv_conflict_threshold = 0.05,
  cnv_directionality_threshold = 0.90,
  min_cnv_event_run_genes = 10L,
  min_idc_novel_event_genes = 10L,
  nearest_distance_tolerance = 1e-8
)

# Please write a script named cnv_clone_metab.R to compare metabolic differences between tumor cells in DCIS and IDC regions for clones shared by both regions across all samples. Tumor cells from all samples have already been scored by GSVA using KEGG metabolic pathways. Read the KEGG pathway GSVA scores with the following code:
# gsva.res = fread('E:/Projects/BRCA/results/ST/GSVA/epi_metacell_gsva_filtered_kegg_res.csv',
#                  nThread = 8,
#                  showProgress = TRUE,
#                  verbose = TRUE,
#                  sep = ",",  # Specify the separator explicitly
#                  header = TRUE)


pattern_labels <- c(
  "1" = "single-lineage invasion",
  "2a" = "lineage-divergent invasion",
  "2b" = "Multiclonal cross-compartment sharing",
  "Indeterminate" = "Indeterminate"
)

# Numeric modes are retained only for compatibility with existing downstream
# scripts that currently expect Mode 1/2/3.
legacy_mode_by_pattern <- c(
  "1" = 1L,
  "2a" = 2L,
  "2b" = 3L,
  "Indeterminate" = NA_integer_
)

# General helpers

first_existing_file <- function(directory, candidates) {
  paths <- file.path(directory, candidates)
  found <- paths[file.exists(paths)]
  if (length(found) == 0) return(NA_character_)
  found[1]
}

read_numeric_matrix <- function(path) {
  x <- read.delim(
    path,
    sep = " ",
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

effective_clone_number <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0 || sum(x) == 0) return(NA_real_)
  p <- x / sum(x)
  exp(-sum(p * log(p)))
}

safe_ratio <- function(numerator, denominator) {
  if (!is.finite(numerator) || !is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

collapse_names <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
  if (length(x) == 0) return("")
  paste(x, collapse = ";")
}

descendant_tip_labels <- function(tree, node) {
  if (node <= length(tree$tip.label)) return(tree$tip.label[node])
  ape::extract.clade(tree, node = node)$tip.label
}

state_is_monophyletic <- function(tree, target_tips, informative_tips) {
  target_tips <- intersect(target_tips, tree$tip.label)
  informative_tips <- intersect(informative_tips, tree$tip.label)

  if (length(target_tips) == 0) return(NA)
  if (length(target_tips) == 1) return(TRUE)

  mrca <- ape::getMRCA(tree, target_tips)
  if (is.null(mrca) || is.na(mrca)) return(FALSE)

  descendants <- descendant_tip_labels(tree, mrca)
  informative_descendants <- intersect(descendants, informative_tips)
  setequal(informative_descendants, target_tips)
}

# Return all tips in the branch directly sister to a target tip or target
# monophyletic clade. The tree must already be rooted (here it is rooted by
# the reference profile). An empty vector means that no direct sister branch
# could be resolved.
direct_sister_tip_labels <- function(tree, target_tips) {
  if (is.null(tree) || length(tree$tip.label) < 2L) return(character(0))

  target_tips <- intersect(as.character(target_tips), tree$tip.label)
  if (length(target_tips) == 0L) return(character(0))

  target_node <- if (length(target_tips) == 1L) {
    match(target_tips, tree$tip.label)
  } else {
    ape::getMRCA(tree, target_tips)
  }
  if (is.null(target_node) || is.na(target_node)) return(character(0))

  parent_rows <- which(tree$edge[, 2] == target_node)
  if (length(parent_rows) != 1L) return(character(0))
  parent_node <- tree$edge[parent_rows, 1]
  sibling_nodes <- tree$edge[
    tree$edge[, 1] == parent_node & tree$edge[, 2] != target_node,
    2
  ]
  if (length(sibling_nodes) == 0L) return(character(0))

  unique(unlist(
    lapply(sibling_nodes, function(node) descendant_tip_labels(tree, node)),
    use.names = FALSE
  ))
}

# Remove short non-neutral runs so isolated noisy genes are not counted as CNV
# bands. Runs are evaluated independently within each chromosome.
filter_short_cnv_runs <- function(states, chromosomes, min_run_genes) {
  states <- as.integer(states)
  chromosomes <- as.character(chromosomes)
  if (length(states) != length(chromosomes)) {
    stop("states and chromosomes must have the same length.")
  }
  if (min_run_genes <= 1L || length(states) == 0L) return(states)

  filtered <- states
  chromosome_levels <- unique(chromosomes)
  for (chromosome in chromosome_levels) {
    indices <- which(chromosomes == chromosome)
    if (length(indices) == 0L) next

    runs <- rle(states[indices])
    run_ends <- cumsum(runs$lengths)
    run_starts <- c(1L, head(run_ends, -1L) + 1L)
    short_non_neutral <- which(
      runs$values != 0L & runs$lengths < min_run_genes
    )
    for (run_index in short_non_neutral) {
      local_indices <- seq.int(run_starts[run_index], run_ends[run_index])
      filtered[indices[local_indices]] <- 0L
    }
  }
  filtered
}

# Fitch parsimony score for one binary character. Mixed tips are treated as
# ambiguous and may take either DCIS or IDC without an extra change.
fitch_parsimony_score <- function(tree, tip_states) {
  if (is.null(tree) || length(tree$tip.label) < 2) return(NA_integer_)

  tip_states <- tip_states[tree$tip.label]
  if (all(is.na(tip_states)) || length(unique(na.omit(tip_states))) < 2) {
    return(0L)
  }

  edge <- tree$edge
  children <- split(edge[, 2], edge[, 1])
  root <- setdiff(edge[, 1], edge[, 2])[1]
  score <- 0L

  visit <- function(node) {
    if (node <= length(tree$tip.label)) {
      state <- tip_states[tree$tip.label[node]]
      if (is.na(state) || state == "Mixed" || state == "Unknown") {
        return(c("DCIS", "IDC"))
      }
      return(state)
    }

    child_sets <- lapply(children[[as.character(node)]], visit)
    current_set <- child_sets[[1]]

    if (length(child_sets) > 1) {
      for (i in 2:length(child_sets)) {
        overlap <- intersect(current_set, child_sets[[i]])
        if (length(overlap) > 0) {
          current_set <- overlap
        } else {
          current_set <- union(current_set, child_sets[[i]])
          score <<- score + 1L
        }
      }
    }
    current_set
  }

  visit(root)
  score
}

# Input preparation and tree reconstruction

read_clone_info <- function(sample_dir) {
  clone_path <- first_existing_file(sample_dir, clone_file_candidates)
  if (is.na(clone_path)) {
    stop(
      "No clone assignment file found for v3 suffix ",
      clone_result_analysis_suffix,
      " or for the legacy unsuffixed names."
    )
  }
  message("  Clone assignment: ", basename(clone_path))

  clone_info <- read.csv(
    clone_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_columns <- c("Barcode", "Node", "region")
  missing_columns <- setdiff(required_columns, colnames(clone_info))
  if (length(missing_columns) > 0) {
    stop("Clone file is missing columns: ", paste(missing_columns, collapse = ", "))
  }

  clone_info <- clone_info[, required_columns, drop = FALSE]
  clone_info$Barcode <- as.character(clone_info$Barcode)
  clone_info$Node <- as.character(clone_info$Node)
  clone_info$region <- as.character(clone_info$region)
  clone_info <- clone_info[
    !is.na(clone_info$Barcode) & nzchar(clone_info$Barcode) &
      !is.na(clone_info$Node) & grepl("^Clone", clone_info$Node),
    ,
    drop = FALSE
  ]
  clone_info <- clone_info[!duplicated(clone_info$Barcode), , drop = FALSE]

  if (nrow(clone_info) == 0) stop("Clone assignment file contains no valid rows.")
  clone_info
}

build_clone_cnv_tree <- function(sample_dir, clone_info) {
  observation_path <- file.path(sample_dir, "0.3.infercnv.observations.txt")
  reference_path <- file.path(sample_dir, "0.3.infercnv.reference.txt")
  gene_order_path <- file.path(sample_dir, "0.3.GeneOrderFile.txt")
  if (!file.exists(observation_path)) stop("Missing inferCNV observation file.")
  if (!file.exists(reference_path)) stop("Missing inferCNV reference file.")
  if (!file.exists(gene_order_path)) stop("Missing inferCNV gene-order file.")

  observation_gene_by_cell <- read_numeric_matrix(observation_path)
  reference_gene_by_cell <- read_numeric_matrix(reference_path)

  colnames(observation_gene_by_cell) <- gsub(
    "\\.", "-", colnames(observation_gene_by_cell)
  )
  colnames(reference_gene_by_cell) <- gsub(
    "\\.", "-", colnames(reference_gene_by_cell)
  )

  cnv <- t(observation_gene_by_cell)
  cnv <- cnv[!duplicated(rownames(cnv)), , drop = FALSE]

  common_barcodes <- intersect(clone_info$Barcode, rownames(cnv))
  if (length(common_barcodes) == 0) {
    stop("No clone barcodes match inferCNV observation barcodes.")
  }

  matched_info <- clone_info[match(common_barcodes, clone_info$Barcode), , drop = FALSE]
  cnv <- cnv[common_barcodes, , drop = FALSE]

  clone_names <- sort(unique(matched_info$Node))
  clone_means <- matrix(
    NA_real_,
    nrow = length(clone_names),
    ncol = ncol(cnv),
    dimnames = list(clone_names, colnames(cnv))
  )
  for (clone_name in clone_names) {
    clone_means[clone_name, ] <- colMeans(
      cnv[matched_info$Node == clone_name, , drop = FALSE],
      na.rm = TRUE
    )
  }

  reference <- t(reference_gene_by_cell)
  common_genes <- intersect(colnames(clone_means), colnames(reference))
  if (length(common_genes) == 0) {
    stop("No genes are shared by observation and reference matrices.")
  }

  clone_means <- clone_means[, common_genes, drop = FALSE]
  reference_mean <- colMeans(reference[, common_genes, drop = FALSE], na.rm = TRUE)

  gene_order <- read.delim(
    gene_order_path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(gene_order) < 2L) {
    stop("Gene-order file must contain at least gene and chromosome columns.")
  }
  gene_order <- gene_order[!duplicated(gene_order[[1]]), , drop = FALSE]
  gene_chromosomes <- as.character(
    gene_order[[2]][match(common_genes, gene_order[[1]])]
  )
  missing_chromosome <- is.na(gene_chromosomes) | !nzchar(gene_chromosomes)
  gene_chromosomes[missing_chromosome] <- paste0(
    "unmapped_", seq_len(sum(missing_chromosome))
  )

  # Match the symmetric display clipping and zero threshold used in v2 before
  # constructing the clone-level NJ tree.
  upper_limit <- max(clone_means, na.rm = TRUE)
  lower_limit <- abs(min(clone_means, na.rm = TRUE))
  symmetric_limit <- min(upper_limit, lower_limit)
  if (is.finite(symmetric_limit) && symmetric_limit > 0) {
    clone_means[clone_means > symmetric_limit] <- symmetric_limit
    clone_means[clone_means < -symmetric_limit] <- -symmetric_limit
  }

  cnv_for_tree <- rbind(clone_means, Ref = reference_mean)
  cnv_for_tree[abs(cnv_for_tree) < params$infercnv_zero_threshold] <- 0

  relative_clone_profiles <- sweep(
    clone_means,
    2,
    reference_mean,
    FUN = "-"
  )
  relative_clone_profiles[
    abs(relative_clone_profiles) < params$infercnv_zero_threshold
  ] <- 0
  clone_event_states <- sign(relative_clone_profiles)
  clone_event_states <- t(vapply(
    seq_len(nrow(clone_event_states)),
    function(row_index) {
      filter_short_cnv_runs(
        clone_event_states[row_index, ],
        gene_chromosomes,
        params$min_cnv_event_run_genes
      )
    },
    integer(ncol(clone_event_states))
  ))
  dimnames(clone_event_states) <- dimnames(relative_clone_profiles)

  tree <- NULL
  distance_matrix <- as.matrix(stats::dist(cnv_for_tree))
  if (nrow(cnv_for_tree) >= 3) {
    tree <- ape::nj(stats::as.dist(distance_matrix))
    tree <- ape::root(tree, outgroup = "Ref", resolve.root = TRUE)
  }

  list(
    tree = tree,
    clone_means = cnv_for_tree[setdiff(rownames(cnv_for_tree), "Ref"), , drop = FALSE],
    reference_mean = cnv_for_tree["Ref", ],
    relative_clone_profiles = relative_clone_profiles,
    clone_event_states = clone_event_states,
    gene_chromosomes = gene_chromosomes,
    distance_matrix = distance_matrix,
    matched_barcodes = common_barcodes
  )
}

# Metrics and classification

calculate_clone_metrics <- function(sample_id, clone_info) {
  clone_names <- sort(unique(clone_info$Node))
  pathology_info <- clone_info[clone_info$region %in% c("DCIS", "IDC"), , drop = FALSE]

  counts <- table(
    factor(pathology_info$Node, levels = clone_names),
    factor(pathology_info$region, levels = c("DCIS", "IDC"))
  )
  counts <- as.matrix(counts)

  total_dcis <- sum(counts[, "DCIS"])
  total_idc <- sum(counts[, "IDC"])
  fraction_dcis <- if (total_dcis > 0) counts[, "DCIS"] / total_dcis else rep(NA_real_, nrow(counts))
  fraction_idc <- if (total_idc > 0) counts[, "IDC"] / total_idc else rep(NA_real_, nrow(counts))

  preference_denominator <- fraction_dcis + fraction_idc
  idc_preference <- fraction_idc / preference_denominator
  idc_preference[!is.finite(idc_preference)] <- NA_real_

  dcis_supported <-
    counts[, "DCIS"] >= params$min_clone_pathology_count &
    fraction_dcis >= params$min_shared_fraction
  idc_supported <-
    counts[, "IDC"] >= params$min_clone_pathology_count &
    fraction_idc >= params$min_shared_fraction
  dcis_supported[is.na(dcis_supported)] <- FALSE
  idc_supported[is.na(idc_supported)] <- FALSE

  is_shared <- dcis_supported & idc_supported
  pathology_state <- rep("Unknown", length(clone_names))
  pathology_state[dcis_supported & !idc_supported] <- "DCIS"
  pathology_state[idc_supported & !dcis_supported] <- "IDC"
  pathology_state[is_shared] <- "Mixed"

  compartment_class <- rep("Unsupported", length(clone_names))
  compartment_class[dcis_supported & !idc_supported] <- "DCIS-private"
  compartment_class[idc_supported & !dcis_supported] <- "IDC-private"
  compartment_class[is_shared] <- "Shared"

  data.frame(
    Sample = sample_id,
    Clone = clone_names,
    DCIS_count = as.integer(counts[, "DCIS"]),
    IDC_count = as.integer(counts[, "IDC"]),
    DCIS_fraction = as.numeric(fraction_dcis),
    IDC_fraction = as.numeric(fraction_idc),
    IDC_preference = as.numeric(idc_preference),
    Pathology_state = pathology_state,
    DCIS_supported = dcis_supported,
    IDC_supported = idc_supported,
    Compartment_class = compartment_class,
    Shared_clone = is_shared,
    stringsAsFactors = FALSE
  )
}

calculate_cnv_metrics <- function(clone_metrics, cnv_result) {
  clone_means <- cnv_result$clone_means
  available_clones <- intersect(clone_metrics$Clone, rownames(clone_means))
  if (length(available_clones) == 0) {
    return(list(
      centroid_distance = NA_real_,
      centroid_distance_ratio = NA_real_,
      dcis_cnv_burden = NA_real_,
      idc_cnv_burden = NA_real_,
      idc_to_dcis_burden_ratio = NA_real_
    ))
  }

  metrics <- clone_metrics[match(available_clones, clone_metrics$Clone), , drop = FALSE]
  means <- clone_means[available_clones, , drop = FALSE]

  dcis_weights <- metrics$DCIS_fraction
  idc_weights <- metrics$IDC_fraction
  dcis_weights[!is.finite(dcis_weights)] <- 0
  idc_weights[!is.finite(idc_weights)] <- 0

  if (sum(dcis_weights) == 0 || sum(idc_weights) == 0) {
    return(list(
      centroid_distance = NA_real_,
      centroid_distance_ratio = NA_real_,
      dcis_cnv_burden = NA_real_,
      idc_cnv_burden = NA_real_,
      idc_to_dcis_burden_ratio = NA_real_
    ))
  }

  dcis_weights <- dcis_weights / sum(dcis_weights)
  idc_weights <- idc_weights / sum(idc_weights)
  dcis_centroid <- colSums(means * dcis_weights)
  idc_centroid <- colSums(means * idc_weights)

  centroid_distance <- sqrt(sum((dcis_centroid - idc_centroid)^2))
  pairwise_clone_distance <- as.matrix(stats::dist(means))
  typical_clone_distance <- if (nrow(means) > 1) {
    stats::median(pairwise_clone_distance[lower.tri(pairwise_clone_distance)], na.rm = TRUE)
  } else {
    0
  }
  centroid_distance_ratio <- if (typical_clone_distance > 0) {
    centroid_distance / typical_clone_distance
  } else if (centroid_distance == 0) {
    0
  } else {
    NA_real_
  }

  reference_mean <- cnv_result$reference_mean
  dcis_burden <- mean(abs(dcis_centroid - reference_mean), na.rm = TRUE)
  idc_burden <- mean(abs(idc_centroid - reference_mean), na.rm = TRUE)

  list(
    centroid_distance = centroid_distance,
    centroid_distance_ratio = centroid_distance_ratio,
    dcis_cnv_burden = dcis_burden,
    idc_cnv_burden = idc_burden,
    idc_to_dcis_burden_ratio = safe_ratio(idc_burden, dcis_burden)
  )
}

# Quantify whether every IDC-private clone is a directionally nested extension
# of the single shared clone and satisfies the requested structural conditions.
calculate_nestedness_metrics <- function(
  sample_id,
  clone_metrics,
  cnv_result
) {
  shared_clones <- clone_metrics$Clone[clone_metrics$Shared_clone]
  idc_private_clones <- clone_metrics$Clone[
    clone_metrics$Compartment_class == "IDC-private"
  ]

  if (length(shared_clones) != 1L || length(idc_private_clones) == 0L) {
    return(data.frame())
  }

  event_states <- cnv_result$clone_event_states
  available_shared <- intersect(shared_clones, rownames(event_states))
  available_idc <- intersect(idc_private_clones, rownames(event_states))
  if (length(available_shared) == 0L || length(available_idc) == 0L) {
    return(data.frame())
  }

  informative_clones <- clone_metrics$Clone[
    clone_metrics$DCIS_supported | clone_metrics$IDC_supported
  ]
  idc_private_monophyletic <- if (length(idc_private_clones) == 1L) {
    TRUE
  } else if (is.null(cnv_result$tree)) {
    NA
  } else {
    tree_without_ref <- ape::drop.tip(cnv_result$tree, "Ref")
    state_is_monophyletic(
      tree_without_ref,
      idc_private_clones,
      informative_clones
    )
  }

  # A distance-nearest clone is not necessarily the direct sister in a rooted
  # topology. Pattern 1 requires the single shared clone to be the only
  # pathology-supported clone in the branch immediately sister to the entire
  # monophyletic IDC-private clade. Thus a nearer intervening DCIS-private
  # branch forces Pattern 2a even when the CNV event sets are nested.
  direct_sister_tips <- if (is.null(cnv_result$tree)) {
    character(0)
  } else {
    direct_sister_tip_labels(cnv_result$tree, idc_private_clones)
  }
  direct_sister_supported_clones <- intersect(
    direct_sister_tips,
    informative_clones
  )
  shared_is_direct_sister <-
    length(direct_sister_supported_clones) == 1L &&
    identical(direct_sister_supported_clones, shared_clones[1])

  dcis_supported_clones <- intersect(
    clone_metrics$Clone[clone_metrics$DCIS_supported],
    rownames(cnv_result$distance_matrix)
  )

  rows <- list()
  row_index <- 0L
  for (shared_clone in available_shared) {
    shared_state <- event_states[shared_clone, ]
    for (idc_clone in available_idc) {
      idc_state <- event_states[idc_clone, ]

      preserved <- sum(shared_state != 0L & idc_state == shared_state)
      shared_only <- sum(shared_state != 0L & idc_state == 0L)
      idc_new <- sum(shared_state == 0L & idc_state != 0L)
      opposite <- sum(
        shared_state != 0L & idc_state != 0L & idc_state != shared_state
      )
      shared_event_total <- preserved + shared_only + opposite
      union_event_total <- preserved + shared_only + idc_new + opposite

      retention_rate <- safe_ratio(preserved, shared_event_total)
      conflict_rate <- safe_ratio(opposite, union_event_total)
      directionality <- safe_ratio(idc_new, idc_new + shared_only)

      nearest_clones <- character(0)
      shared_is_unique_nearest <- FALSE
      if (idc_clone %in% rownames(cnv_result$distance_matrix) &&
          length(dcis_supported_clones) > 0L) {
        dcis_distances <- cnv_result$distance_matrix[
          idc_clone,
          dcis_supported_clones,
          drop = TRUE
        ]
        minimum_distance <- min(dcis_distances, na.rm = TRUE)
        tolerance <- params$nearest_distance_tolerance * max(
          1,
          abs(minimum_distance)
        )
        nearest_clones <- names(dcis_distances)[
          abs(dcis_distances - minimum_distance) <= tolerance
        ]
        shared_is_unique_nearest <-
          length(nearest_clones) == 1L && identical(nearest_clones, shared_clone)
      }

      metric_pass <-
        is.finite(retention_rate) &&
        retention_rate >= params$cnv_retention_threshold &&
        is.finite(conflict_rate) &&
        conflict_rate <= params$cnv_conflict_threshold &&
        is.finite(directionality) &&
        directionality >= params$cnv_directionality_threshold &&
        idc_new >= params$min_idc_novel_event_genes

      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        Sample = sample_id,
        Shared_clone = shared_clone,
        IDC_private_clone = idc_clone,
        Preserved_event_genes = preserved,
        Shared_only_event_genes = shared_only,
        IDC_new_event_genes = idc_new,
        Opposite_event_genes = opposite,
        Union_event_genes = union_event_total,
        CNV_retention_rate = retention_rate,
        CNV_conflict_rate = conflict_rate,
        CNV_directionality = directionality,
        Nearest_DCIS_supported_clone = collapse_names(nearest_clones),
        Shared_is_unique_nearest_DCIS_clone = shared_is_unique_nearest,
        Direct_sister_tip_clones = collapse_names(direct_sister_tips),
        Direct_sister_supported_clones = collapse_names(
          direct_sister_supported_clones
        ),
        Shared_is_direct_sister_of_IDC_clade = shared_is_direct_sister,
        IDC_private_clones_monophyletic = idc_private_monophyletic,
        CNV_metric_pass = metric_pass,
        Pair_structure_pass =
          shared_is_unique_nearest && shared_is_direct_sister,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  result <- do.call(rbind, rows)
  single_shared_clone <- length(shared_clones) == 1L
  all_idc_pairs_present <-
    length(unique(result$IDC_private_clone)) == length(idc_private_clones)
  same_shared_clone <-
    single_shared_clone && all(result$Shared_clone == shared_clones[1])
  sample_pass <-
    single_shared_clone &&
    all_idc_pairs_present &&
    same_shared_clone &&
    isTRUE(idc_private_monophyletic) &&
    isTRUE(shared_is_direct_sister) &&
    all(result$CNV_metric_pass) &&
    all(result$Pair_structure_pass)

  result$All_IDC_nested_single_lineage <- sample_pass
  result
}

classify_pattern <- function(metrics) {
  sparse_pathology <-
    metrics$total_dcis < params$min_pathology_metacells ||
    metrics$total_idc < params$min_pathology_metacells

  # Do not force a three-pattern assignment when either pathological
  # compartment is insufficiently sampled or no clone passes both the
  # absolute-count and within-compartment fraction thresholds.
  if (sparse_pathology) {
    sparse_compartments <- c(
      if (metrics$total_dcis < params$min_pathology_metacells) "DCIS",
      if (metrics$total_idc < params$min_pathology_metacells) "IDC"
    )
    return(list(
      pattern = "Indeterminate",
      confidence = "Not assigned",
      rule = paste0(
        "Insufficient total metacell support in ",
        paste(sparse_compartments, collapse = " and "),
        " (required >= ", params$min_pathology_metacells,
        " metacells per compartment)."
      )
    ))
  }

  if (metrics$n_dcis_supported_clones == 0L ||
      metrics$n_idc_supported_clones == 0L) {
    unsupported_compartments <- c(
      if (metrics$n_dcis_supported_clones == 0L) "DCIS",
      if (metrics$n_idc_supported_clones == 0L) "IDC"
    )
    return(list(
      pattern = "Indeterminate",
      confidence = "Not assigned",
      rule = paste0(
        "No clone in ", paste(unsupported_compartments, collapse = " or "),
        " passed both support thresholds (count >= ",
        params$min_clone_pathology_count,
        "; fraction >= ", params$min_shared_fraction, ")."
      )
    ))
  }

  # This direct observation takes precedence over tree-based interpretation:
  # at least two clones independently occupy both pathological compartments.
  if (metrics$n_shared_clones >= 2) {
    return(list(
      pattern = "2b",
      confidence = if (sparse_pathology) "Low" else "High",
      rule = "At least two supported CNV clones are shared across DCIS and IDC."
    ))
  }

  if (metrics$n_shared_clones == 1) {
    if (metrics$n_idc_private_clones >= 1) {
      if (isTRUE(metrics$all_idc_nested_single_lineage)) {
        return(list(
          pattern = "1",
          confidence = if (sparse_pathology) "Low" else "Medium",
          rule = paste0(
            "The monophyletic IDC-private clade is a directionally nested CNV ",
            "extension of its direct sister shared clone (retention >= ",
            params$cnv_retention_threshold,
            ", conflict <= ", params$cnv_conflict_threshold,
            ", directionality >= ", params$cnv_directionality_threshold,
            "). The direct sister branch contains no other pathology-supported ",
            "clone, and the shared clone is also every IDC-private clone's ",
            "unique nearest DCIS-supported clone in CNV-distance space."
          )
        ))
      }
      if (identical(metrics$shared_is_direct_sister_of_idc_clade, FALSE)) {
        return(list(
          pattern = "2a",
          confidence = if (sparse_pathology) "Low" else "High",
          rule = paste0(
            "One shared clone coexists with at least one IDC-private clone, ",
            "but the shared clone is not the direct sister lineage of the ",
            "monophyletic IDC-private clade; an intervening supported lineage ",
            "is therefore consistent with lineage-divergent evolution."
          )
        ))
      }
      return(list(
        pattern = "2a",
        confidence = if (sparse_pathology) "Low" else "Medium",
        rule = paste0(
          "One cross-compartment shared clone coexists with at least one ",
          "IDC-private clone, but the complete CNV-nestedness and structural ",
          "criteria for a single invasion-associated lineage are not satisfied."
        )
      ))
    }
    return(list(
      pattern = "1",
      confidence = if (sparse_pathology) "Low" else "High",
      rule = paste0(
        "Exactly one clone spans DCIS and IDC and no additional IDC-private ",
        "clone is supported; additional DCIS-private clones do not add an invasion lineage."
      )
    ))
  }

  # With no shared clone, one IDC-private clone or one monophyletic IDC clade is
  # interpreted as a single invasion lineage. Multiple separated IDC branches
  # imply repeated independent invasion.
  if (metrics$n_idc_private_clones <= 1) {
    return(list(
      pattern = "1",
      confidence = if (sparse_pathology) "Low" else "Medium",
      rule = "No shared clone is detected and IDC is represented by no more than one supported private clone."
    ))
  }

  if (isTRUE(metrics$idc_monophyletic) &&
      (!is.finite(metrics$parsimony_changes) || metrics$parsimony_changes <= 1)) {
    return(list(
      pattern = "1",
      confidence = if (sparse_pathology) "Low" else "Medium",
      rule = "Multiple IDC-private clones form one monophyletic clade, consistent with one invasion lineage followed by diversification."
    ))
  }

  if (identical(metrics$idc_monophyletic, FALSE) ||
      (is.finite(metrics$parsimony_changes) && metrics$parsimony_changes >= 2)) {
    return(list(
      pattern = "2a",
      confidence = if (sparse_pathology) "Low" else "High",
      rule = "IDC-private clones occupy separated tree branches, requiring at least two parsimonious DCIS-IDC transitions."
    ))
  }

  list(
    pattern = "2a",
    confidence = "Low",
    rule = "Multiple IDC-private clones are present, but their phylogenetic separation cannot be resolved."
  )
}

empty_summary_row <- function(sample_id, error_message = "") {
  data.frame(
    Sample = sample_id,
    Pattern = NA_character_,
    Pattern_label = "Input/analysis failure",
    Mode = NA_integer_,
    Mode_label = "Input/analysis failure",
    Classification_status = "Input/analysis failure",
    Confidence = "Low",
    Rule = if (nzchar(error_message)) paste("Input/analysis error:", error_message) else "Indeterminate",
    Total_metacells = NA_integer_,
    DCIS_metacells = NA_integer_,
    IDC_metacells = NA_integer_,
    Clone_count = NA_integer_,
    Informative_clone_count = NA_integer_,
    DCIS_supported_clone_count = NA_integer_,
    DCIS_supported_clones = "",
    IDC_supported_clone_count = NA_integer_,
    IDC_supported_clones = "",
    DCIS_supported_metacells = NA_integer_,
    DCIS_supported_metacell_fraction = NA_real_,
    IDC_supported_metacells = NA_integer_,
    IDC_supported_metacell_fraction = NA_real_,
    Unsupported_clone_count = NA_integer_,
    Unsupported_clones = "",
    Shared_clone_count = NA_integer_,
    Shared_clones = "",
    DCIS_private_clone_count = NA_integer_,
    IDC_private_clone_count = NA_integer_,
    DCIS_state_clone_count = NA_integer_,
    IDC_state_clone_count = NA_integer_,
    Effective_DCIS_clones = NA_real_,
    Effective_IDC_clones = NA_real_,
    Bottleneck_ratio = NA_real_,
    Dominant_DCIS_clone = "",
    Dominant_IDC_clone = "",
    Same_dominant_clone = NA,
    Max_IDC_clone_fraction = NA_real_,
    IDC_dominant_clone_DCIS_fraction = NA_real_,
    IDC_monophyletic = NA,
    DCIS_monophyletic = NA,
    Parsimony_changes = NA_integer_,
    Nested_single_lineage_pass = NA,
    Nested_shared_clone = "",
    Nested_IDC_private_clones = "",
    IDC_private_clones_monophyletic_with_shared = NA,
    All_IDC_nearest_to_shared = NA,
    IDC_clade_direct_sister_tip_clones = "",
    IDC_clade_direct_sister_supported_clones = "",
    Shared_is_direct_sister_of_IDC_clade = NA,
    Minimum_CNV_retention_rate = NA_real_,
    Maximum_CNV_conflict_rate = NA_real_,
    Minimum_CNV_directionality = NA_real_,
    DCIS_IDC_centroid_distance = NA_real_,
    Centroid_distance_ratio = NA_real_,
    DCIS_CNV_burden = NA_real_,
    IDC_CNV_burden = NA_real_,
    IDC_to_DCIS_burden_ratio = NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

analyze_sample <- function(sample_id) {
  sample_dir <- file.path(infercnv_root, sample_id)
  message("Analyzing ", sample_id, " ...")

  clone_info <- read_clone_info(sample_dir)
  clone_metrics <- calculate_clone_metrics(sample_id, clone_info)
  cnv_result <- build_clone_cnv_tree(sample_dir, clone_info)
  cnv_metrics <- calculate_cnv_metrics(clone_metrics, cnv_result)
  nestedness_metrics <- calculate_nestedness_metrics(
    sample_id,
    clone_metrics,
    cnv_result
  )

  total_dcis <- sum(clone_metrics$DCIS_count)
  total_idc <- sum(clone_metrics$IDC_count)
  informative <- clone_metrics$DCIS_supported | clone_metrics$IDC_supported
  dcis_supported <- clone_metrics$DCIS_supported
  idc_supported <- clone_metrics$IDC_supported
  unsupported <- !informative

  n_dcis_supported_clones <- sum(dcis_supported)
  n_idc_supported_clones <- sum(idc_supported)
  dcis_supported_metacells <- sum(
    clone_metrics$DCIS_count[dcis_supported]
  )
  idc_supported_metacells <- sum(
    clone_metrics$IDC_count[idc_supported]
  )
  dcis_supported_metacell_fraction <- safe_ratio(
    dcis_supported_metacells,
    total_dcis
  )
  idc_supported_metacell_fraction <- safe_ratio(
    idc_supported_metacells,
    total_idc
  )

  effective_dcis <- effective_clone_number(clone_metrics$DCIS_count)
  effective_idc <- effective_clone_number(clone_metrics$IDC_count)
  bottleneck_ratio <- safe_ratio(effective_idc, effective_dcis)

  dominant_dcis_index <- if (total_dcis > 0) which.max(clone_metrics$DCIS_fraction) else NA_integer_
  dominant_idc_index <- if (total_idc > 0) which.max(clone_metrics$IDC_fraction) else NA_integer_
  dominant_dcis_clone <- if (!is.na(dominant_dcis_index)) clone_metrics$Clone[dominant_dcis_index] else ""
  dominant_idc_clone <- if (!is.na(dominant_idc_index)) clone_metrics$Clone[dominant_idc_index] else ""
  same_dominant <- nzchar(dominant_dcis_clone) && dominant_dcis_clone == dominant_idc_clone

  dcis_tips <- clone_metrics$Clone[clone_metrics$Pathology_state == "DCIS"]
  idc_tips <- clone_metrics$Clone[clone_metrics$Pathology_state == "IDC"]
  informative_tips <- c(dcis_tips, idc_tips)

  dcis_monophyletic <- NA
  idc_monophyletic <- NA
  parsimony_changes <- NA_integer_
  classification_tree <- cnv_result$tree

  if (!is.null(classification_tree)) {
    tree_without_ref <- ape::drop.tip(classification_tree, "Ref")
    dcis_monophyletic <- state_is_monophyletic(
      tree_without_ref, dcis_tips, informative_tips
    )
    idc_monophyletic <- state_is_monophyletic(
      tree_without_ref, idc_tips, informative_tips
    )
    tip_states <- setNames(clone_metrics$Pathology_state, clone_metrics$Clone)
    parsimony_changes <- fitch_parsimony_score(tree_without_ref, tip_states)
  }

  idc_dominant_clone_dcis_count <- if (!is.na(dominant_idc_index)) {
    clone_metrics$DCIS_count[dominant_idc_index]
  } else {
    NA_integer_
  }
  idc_dominant_clone_dcis_fraction <- if (!is.na(dominant_idc_index)) {
    clone_metrics$DCIS_fraction[dominant_idc_index]
  } else {
    NA_real_
  }

  nested_single_lineage_pass <-
    nrow(nestedness_metrics) > 0L &&
    all(nestedness_metrics$All_IDC_nested_single_lineage)
  nested_shared_clone <- if (nrow(nestedness_metrics) > 0L) {
    collapse_names(nestedness_metrics$Shared_clone)
  } else {
    ""
  }
  nested_idc_private_clones <- if (nrow(nestedness_metrics) > 0L) {
    collapse_names(nestedness_metrics$IDC_private_clone)
  } else {
    ""
  }
  nested_idc_monophyletic <- if (nrow(nestedness_metrics) > 0L) {
    isTRUE(nestedness_metrics$IDC_private_clones_monophyletic[1])
  } else {
    NA
  }
  all_idc_nearest_to_shared <- if (nrow(nestedness_metrics) > 0L) {
    all(nestedness_metrics$Shared_is_unique_nearest_DCIS_clone)
  } else {
    NA
  }
  direct_sister_tip_clones <- if (nrow(nestedness_metrics) > 0L) {
    collapse_names(nestedness_metrics$Direct_sister_tip_clones)
  } else {
    ""
  }
  direct_sister_supported_clones <- if (nrow(nestedness_metrics) > 0L) {
    collapse_names(nestedness_metrics$Direct_sister_supported_clones)
  } else {
    ""
  }
  shared_is_direct_sister <- if (nrow(nestedness_metrics) > 0L) {
    isTRUE(nestedness_metrics$Shared_is_direct_sister_of_IDC_clade[1])
  } else {
    NA
  }
  minimum_retention <- if (nrow(nestedness_metrics) > 0L) {
    safe_min(nestedness_metrics$CNV_retention_rate)
  } else {
    NA_real_
  }
  maximum_conflict <- if (nrow(nestedness_metrics) > 0L) {
    safe_max(nestedness_metrics$CNV_conflict_rate)
  } else {
    NA_real_
  }
  minimum_directionality <- if (nrow(nestedness_metrics) > 0L) {
    safe_min(nestedness_metrics$CNV_directionality)
  } else {
    NA_real_
  }

  metrics <- list(
    total_dcis = total_dcis,
    total_idc = total_idc,
    n_informative_clones = sum(informative),
    n_dcis_supported_clones = n_dcis_supported_clones,
    n_idc_supported_clones = n_idc_supported_clones,
    n_shared_clones = sum(clone_metrics$Shared_clone),
    n_dcis_private_clones = sum(clone_metrics$Compartment_class == "DCIS-private"),
    n_idc_private_clones = sum(clone_metrics$Compartment_class == "IDC-private"),
    n_dcis_state_clones = length(dcis_tips),
    n_idc_state_clones = length(idc_tips),
    effective_dcis_clones = effective_dcis,
    effective_idc_clones = effective_idc,
    bottleneck_ratio = bottleneck_ratio,
    max_idc_clone_fraction = if (total_idc > 0) max(clone_metrics$IDC_fraction, na.rm = TRUE) else NA_real_,
    idc_dominant_clone_dcis_count = idc_dominant_clone_dcis_count,
    idc_dominant_clone_dcis_fraction = idc_dominant_clone_dcis_fraction,
    same_dominant_clone = same_dominant,
    dcis_monophyletic = dcis_monophyletic,
    idc_monophyletic = idc_monophyletic,
    parsimony_changes = parsimony_changes,
    all_idc_nested_single_lineage = nested_single_lineage_pass,
    shared_is_direct_sister_of_idc_clade = shared_is_direct_sister,
    centroid_distance = cnv_metrics$centroid_distance,
    centroid_distance_ratio = cnv_metrics$centroid_distance_ratio
  )

  assignment <- classify_pattern(metrics)
  legacy_mode <- unname(legacy_mode_by_pattern[assignment$pattern])
  pattern_label <- unname(pattern_labels[assignment$pattern])
  classification_status <- if (assignment$pattern == "Indeterminate") {
    "Indeterminate"
  } else {
    "Assigned"
  }

  summary_row <- data.frame(
    Sample = sample_id,
    Pattern = assignment$pattern,
    Pattern_label = pattern_label,
    Mode = legacy_mode,
    Mode_label = pattern_label,
    Classification_status = classification_status,
    Confidence = assignment$confidence,
    Rule = assignment$rule,
    Total_metacells = nrow(clone_info),
    DCIS_metacells = total_dcis,
    IDC_metacells = total_idc,
    Clone_count = nrow(clone_metrics),
    Informative_clone_count = sum(informative),
    DCIS_supported_clone_count = n_dcis_supported_clones,
    DCIS_supported_clones = collapse_names(
      clone_metrics$Clone[dcis_supported]
    ),
    IDC_supported_clone_count = n_idc_supported_clones,
    IDC_supported_clones = collapse_names(
      clone_metrics$Clone[idc_supported]
    ),
    DCIS_supported_metacells = dcis_supported_metacells,
    DCIS_supported_metacell_fraction = dcis_supported_metacell_fraction,
    IDC_supported_metacells = idc_supported_metacells,
    IDC_supported_metacell_fraction = idc_supported_metacell_fraction,
    Unsupported_clone_count = sum(unsupported),
    Unsupported_clones = collapse_names(
      clone_metrics$Clone[unsupported]
    ),
    Shared_clone_count = sum(clone_metrics$Shared_clone),
    Shared_clones = collapse_names(clone_metrics$Clone[clone_metrics$Shared_clone]),
    DCIS_private_clone_count = sum(clone_metrics$Compartment_class == "DCIS-private"),
    IDC_private_clone_count = sum(clone_metrics$Compartment_class == "IDC-private"),
    DCIS_state_clone_count = length(dcis_tips),
    IDC_state_clone_count = length(idc_tips),
    Effective_DCIS_clones = effective_dcis,
    Effective_IDC_clones = effective_idc,
    Bottleneck_ratio = bottleneck_ratio,
    Dominant_DCIS_clone = dominant_dcis_clone,
    Dominant_IDC_clone = dominant_idc_clone,
    Same_dominant_clone = same_dominant,
    Max_IDC_clone_fraction = metrics$max_idc_clone_fraction,
    IDC_dominant_clone_DCIS_fraction = idc_dominant_clone_dcis_fraction,
    IDC_monophyletic = idc_monophyletic,
    DCIS_monophyletic = dcis_monophyletic,
    Parsimony_changes = parsimony_changes,
    Nested_single_lineage_pass = nested_single_lineage_pass,
    Nested_shared_clone = nested_shared_clone,
    Nested_IDC_private_clones = nested_idc_private_clones,
    IDC_private_clones_monophyletic_with_shared = nested_idc_monophyletic,
    All_IDC_nearest_to_shared = all_idc_nearest_to_shared,
    IDC_clade_direct_sister_tip_clones = direct_sister_tip_clones,
    IDC_clade_direct_sister_supported_clones = direct_sister_supported_clones,
    Shared_is_direct_sister_of_IDC_clade = shared_is_direct_sister,
    Minimum_CNV_retention_rate = minimum_retention,
    Maximum_CNV_conflict_rate = maximum_conflict,
    Minimum_CNV_directionality = minimum_directionality,
    DCIS_IDC_centroid_distance = cnv_metrics$centroid_distance,
    Centroid_distance_ratio = cnv_metrics$centroid_distance_ratio,
    DCIS_CNV_burden = cnv_metrics$dcis_cnv_burden,
    IDC_CNV_burden = cnv_metrics$idc_cnv_burden,
    IDC_to_DCIS_burden_ratio = cnv_metrics$idc_to_dcis_burden_ratio,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(
    summary = summary_row,
    clone_metrics = clone_metrics,
    nestedness_metrics = nestedness_metrics,
    tree = cnv_result$tree
  )
}

##### Batch analysis #####

if (!dir.exists(infercnv_root)) {
  stop("infercnv_root does not exist: ", infercnv_root)
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (is.null(sample_ids)) {
  candidate_dirs <- list.dirs(infercnv_root, recursive = FALSE, full.names = TRUE)
  has_clone_file <- vapply(
    candidate_dirs,
    function(path) !is.na(first_existing_file(path, clone_file_candidates)),
    logical(1)
  )
  sample_ids <- basename(candidate_dirs[has_clone_file])
}
sample_ids <- sort(unique(sample_ids))

if (length(sample_ids) == 0) stop("No samples with clone result files were found.")

results <- vector("list", length(sample_ids))
names(results) <- sample_ids

for (sample_id in sample_ids) {
  results[[sample_id]] <- tryCatch(
    analyze_sample(sample_id),
    error = function(e) {
      message("  Failed: ", conditionMessage(e))
      list(
        summary = empty_summary_row(sample_id, conditionMessage(e)),
        clone_metrics = NULL,
        nestedness_metrics = NULL,
        tree = NULL
      )
    }
  )
}

summary_table <- do.call(rbind, lapply(results, `[[`, "summary"))
rownames(summary_table) <- NULL

clone_metric_list <- lapply(results, `[[`, "clone_metrics")
clone_metric_list <- clone_metric_list[!vapply(clone_metric_list, is.null, logical(1))]
clone_metric_table <- if (length(clone_metric_list) > 0) {
  do.call(rbind, clone_metric_list)
} else {
  data.frame()
}
rownames(clone_metric_table) <- NULL

nestedness_metric_list <- lapply(results, `[[`, "nestedness_metrics")
nestedness_metric_list <- nestedness_metric_list[
  !vapply(
    nestedness_metric_list,
    function(x) is.null(x) || nrow(x) == 0L,
    logical(1)
  )
]
nestedness_metric_table <- if (length(nestedness_metric_list) > 0L) {
  do.call(rbind, nestedness_metric_list)
} else {
  data.frame()
}
rownames(nestedness_metric_table) <- NULL

tree_list <- lapply(results, `[[`, "tree")
tree_list <- tree_list[!vapply(tree_list, is.null, logical(1))]

parameter_table <- data.frame(
  Parameter = names(params),
  Value = vapply(params, as.character, character(1)),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "CNV_clone_evolution_mode_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  clone_metric_table,
  file.path(output_dir, "CNV_clone_evolution_clone_metrics.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  parameter_table,
  file.path(output_dir, "CNV_clone_evolution_parameters.csv"),
  row.names = FALSE
)
write.csv(
  nestedness_metric_table,
  file.path(output_dir, "CNV_clone_evolution_nestedness_metrics.csv"),
  row.names = FALSE,
  na = ""
)
saveRDS(tree_list, file.path(output_dir, "CNV_clone_evolution_trees.rds"))

message("Finished. Pattern counts:")
print(table(summary_table$Pattern, summary_table$Pattern_label, useNA = "ifany"))
message("Summary: ", file.path(output_dir, "CNV_clone_evolution_mode_summary.csv"))
