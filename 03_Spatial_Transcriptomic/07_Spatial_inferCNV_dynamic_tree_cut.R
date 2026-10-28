library(ape)
library(phylogram)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialInferCNV)
library(tidyverse)
library(ggtree)
library(cluster)
library(dynamicTreeCut)

sample_list_cnv = c(
  "P1","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13","P14","P15","P16","P26","P27","P30","P32","P33","P34","P35", 
  "P36","P37","P38","P39","P40_1","P42","P43_2","P44","P45","P47"
)


nodes_to_delete_all = 
  list(
    "P1" = c(3282),
    "P3" = c(3),
    "P6"=c(),
    "P8"=c(3),
    "P10"=c(),
    "P11"=c(),
    "P12"=c(),
    "P15"=c(),
    "P4"=c(),
    "P5"=c(),
    "P7"=c(7509),
    "P8"=c(),
    "P13"=c(),
    "P14"=c(),
    "P16"=c(),
    "P26"=c(),
    "P27"=c(),
    "P30"=c(),
    "P32"=c(3),
    "P33"=c(197),
    "P34"=c(3),
    "P35"=c(),
    "P36"=c(2),
    "P37"=c(7901),
    "P38"=c(2),
    "P39"=c(),
    "P40_1"=c(),
    "P42"=c(2208),
    "P43_2"=c(2),
    "P44"=c(),
    ""=c(4319),
    "P47"=c()
  )

for(sample_id in sample_list_cnv){
  print(sample_id)
  X_outpath=paste0("E:\\Projects\\BRCA\\results\\ST\\infercnv_metacell\\",sample_id,"\\0.3.X_plot.txt")
  X_plot=read.table(X_outpath,row.names = 1,header = T)
  
  setwd(paste0('E:\\Projects\\BRCA\\results\\ST\\infercnv_metacell\\',sample_id))
  dendrogram <- read.dendrogram(file = "./0.3.obs_dendro_Zmod.txt")
  mytree <- as.phylo(dendrogram)
  my.subtrees <- subtrees(mytree)
  
  save(X_plot, dendrogram, mytree, my.subtrees,
       file = paste0("E:\\Projects\\BRCA\\results\\ST\\infercnv_metacell\\", sample_id, "_subtrees.Rdata"))
}




for(spID in sample_list_cnv[1:length(sample_list_cnv)]){
  delete_nodes_input = nodes_to_delete_all[[spID]]
  clone_definition_and_plot(sample_id = spID, delete_clone_from_full_clones = c(), delete_tree_nodes = delete_nodes_input,
                            dynamic_deep_split = 0.1, min_clone_metacells = 50)

}



get_deleted_barcodes_from_tree_nodes <- function(my.subtrees, mytree, delete_tree_nodes) {
  if (length(delete_tree_nodes) == 0) {
    return(data.frame(DeleteNode = integer(), RowIndex = integer(), Barcode = character()))
  }
  
  delete_tree_nodes <- unique(delete_tree_nodes)
  tree_tip_labels <- as.character(mytree$tip.label)
  
  deleted_barcodes_list <- lapply(delete_tree_nodes, function(root_id) {
    subtree_df <- tryCatch(
      SelectingSubTreeData(my.subtrees, root_id),
      error = function(e) {
        stop(
          "Failed to select tree node ", root_id,
          ". Please check whether this node id exists in my.subtrees. Original error: ",
          conditionMessage(e)
        )
      }
    )
    
    if (!"Barcode" %in% colnames(subtree_df)) {
      stop("SelectingSubTreeData() result for tree node ", root_id, " does not contain a Barcode column.")
    }
    
    subtree_barcodes <- as.character(subtree_df$Barcode)
    row_index <- match(subtree_barcodes, tree_tip_labels)
    
    if (any(is.na(row_index))) {
      barcode_as_index <- suppressWarnings(as.integer(subtree_barcodes))
      valid_index <- !is.na(barcode_as_index) &
        barcode_as_index >= 1 &
        barcode_as_index <= length(tree_tip_labels)
      row_index[is.na(row_index) & valid_index] <- barcode_as_index[is.na(row_index) & valid_index]
    }
    
    if (any(is.na(row_index))) {
      stop(
        "Could not map all barcodes from tree node ", root_id,
        " to tree tip order. Please check mytree$tip.label and SelectingSubTreeData() output."
      )
    }
    
    data.frame(
      DeleteNode = root_id,
      RowIndex = row_index,
      Barcode = tree_tip_labels[row_index],
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(deleted_barcodes_list) %>%
    distinct(DeleteNode, RowIndex, Barcode)
}


# Extract the numeric part of labels such as "Clone C1". Non-clone labels are
# returned as Inf so that they can be placed after numbered clones.
get_clone_display_number <- function(labels) {
  labels <- as.character(labels)
  is_clone <- grepl("^Clone C[0-9]+$", labels)
  clone_number <- rep(Inf, length(labels))
  clone_number[is_clone] <- as.numeric(
    sub("^Clone C", "", labels[is_clone])
  )
  clone_number
}


# Rotate internal nodes only, keeping topology and branch lengths unchanged,
# so lower-numbered clones appear as low as the topology permits. ape::
# rotateConstr() interprets the constraint from bottom to top.
rotate_clone_tree_for_display <- function(tree, reference_label = "Ref") {
  if (is.null(tree) || length(tree$tip.label) <= 1) return(tree)

  tip_labels <- as.character(tree$tip.label)
  clone_labels <- tip_labels[is.finite(get_clone_display_number(tip_labels))]
  clone_labels <- clone_labels[order(get_clone_display_number(clone_labels))]
  other_labels <- setdiff(tip_labels, c(clone_labels, reference_label))
  desired_bottom_to_top <- c(
    clone_labels,
    sort(other_labels),
    intersect(reference_label, tip_labels)
  )

  if (!setequal(desired_bottom_to_top, tip_labels)) {
    stop("Failed to construct a complete clone-tree display order.")
  }

  rotated_tree <- ape::rotateConstr(tree, desired_bottom_to_top)
  attr(rotated_tree, "desired_tip_order_bottom_to_top") <- desired_bottom_to_top
  rotated_tree
}


clone_definition_and_plot = function(sample_id, delete_clone_from_full_clones = c(), delete_tree_nodes = c(),
                                     dynamic_deep_split = 0.9, min_clone_metacells = 50,
                                     discretize_cnv_pattern = FALSE, tree_method = "NJ"){
  message(paste0("Processing sample ", sample_id))
  
  if (length(dynamic_deep_split) != 1 || !is.finite(dynamic_deep_split)) {
    stop("dynamic_deep_split must be one finite numeric value.")
  }
  if (length(min_clone_metacells) != 1 || !is.finite(min_clone_metacells) || min_clone_metacells < 1) {
    stop("min_clone_metacells must be one positive number.")
  }
  if (length(discretize_cnv_pattern) != 1 || !is.logical(discretize_cnv_pattern) || is.na(discretize_cnv_pattern)) {
    stop("discretize_cnv_pattern must be either TRUE or FALSE.")
  }
  if (length(tree_method) != 1 || is.na(tree_method)) {
    stop("tree_method must be either 'NJ' or 'MP'.")
  }
  tree_method <- toupper(as.character(tree_method))
  if (!tree_method %in% c("NJ", "MP")) {
    stop("tree_method must be either 'NJ' or 'MP'.")
  }
  if (tree_method == "MP" && !discretize_cnv_pattern) {
    stop(
      "Maximum parsimony requires categorical CNV states. Set ",
      "discretize_cnv_pattern = TRUE when tree_method = 'MP'."
    )
  }
  min_clone_metacells <- as.integer(min_clone_metacells)
  dynamic_split_label <- format(dynamic_deep_split, trim = TRUE, scientific = FALSE)
  dynamic_split_suffix <- paste0("_split", dynamic_split_label)
  cnv_pattern_suffix <- if (discretize_cnv_pattern) "_cnv3State" else "_cnvContinuous"
  tree_method_suffix <- paste0("_tree", tree_method)
  analysis_suffix <- paste0(dynamic_split_suffix, cnv_pattern_suffix, tree_method_suffix)

  region = read.csv(paste0("E:/Projects/BRCA/results/ST/SuperCell/metacell_region_all.csv"), 
                    row.names = 1, header = T)
  rownames(region) = region$meta_cell
  
  if(length(delete_clone_from_full_clones) == 0){
    message("First run: calculating the optimal number of clones from scratch")
    
    setwd(paste0('E:/Projects/BRCA/results/ST/infercnv_metacell/',sample_id))
    
    load(paste0("E:/Projects/BRCA/results/ST/infercnv_metacell/", sample_id, "_subtrees.Rdata"))
    
    hc <- as.hclust(dendrogram)
    if (nrow(X_plot) == length(mytree$tip.label)) {
      barcode_for_rows <- as.character(mytree$tip.label)
    } else if (nrow(X_plot) == length(hc$labels)) {
      barcode_for_rows <- as.character(hc$labels)
    } else {
      stop(
        "Cannot map X_plot rows to true barcodes for sample ", sample_id,
        ": nrow(X_plot) = ", nrow(X_plot),
        ", length(mytree$tip.label) = ", length(mytree$tip.label),
        ", length(hc$labels) = ", length(hc$labels), "."
      )
    }
    rownames(X_plot) <- barcode_for_rows
    missing_region_barcodes <- setdiff(barcode_for_rows, rownames(region))
    if (length(missing_region_barcodes) > 0) {
      stop(
        "Some mapped barcodes are not found in region for sample ", sample_id,
        ". Please check whether X_plot row order matches the dendrogram tip order. Examples: ",
        paste(head(missing_region_barcodes, 5), collapse = ", ")
      )
    }
    
    if (length(delete_tree_nodes) > 0) {
      deleted_node_barcodes <- get_deleted_barcodes_from_tree_nodes(my.subtrees, mytree, delete_tree_nodes)
      deleted_row_indices <- sort(unique(deleted_node_barcodes$RowIndex))
      deleted_barcodes <- unique(deleted_node_barcodes$Barcode)
      
      if (length(deleted_row_indices) == 0) {
        stop("No metacell rows from delete_tree_nodes were found in X_plot for sample ", sample_id, ".")
      }
      
      keep_rows <- setdiff(seq_len(nrow(X_plot)), deleted_row_indices)
      if (length(keep_rows) < 2) {
        stop(
          "delete_tree_nodes removed too many metacells for sample ", sample_id,
          ". Remaining metacell number is ", length(keep_rows), "."
        )
      }
      
      write.csv(
        deleted_node_barcodes,
        paste0("./Deleted_tree_node_barcodes_dynamic_tree_cut_v3", analysis_suffix, ".csv"),
        row.names = FALSE
      )
      
      X_plot <- X_plot[keep_rows, , drop = FALSE]
      barcode_for_rows <- barcode_for_rows[keep_rows]
      rownames(X_plot) <- barcode_for_rows
      deleted_tree_tips <- intersect(deleted_barcodes, mytree$tip.label)
      if (length(deleted_tree_tips) > 0) {
        mytree <- drop.tip(mytree, deleted_tree_tips)
        my.subtrees <- subtrees(mytree)
      }
      message(
        "Deleted ", length(deleted_barcodes), " metacells from tree node(s): ",
        paste(delete_tree_nodes, collapse = ", "),
        ". Remaining metacells: ", nrow(X_plot)
      )
      
      dist_mat <- dist(X_plot)
      hclust_method <- hc$method
      valid_hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
      if (is.null(hclust_method) || is.na(hclust_method) || !(hclust_method %in% valid_hclust_methods)) {
        hclust_method <- "ward.D"
      }
      hc <- hclust(dist_mat, method = hclust_method)
      hc$labels <- barcode_for_rows
      message("Rebuilt hclust after tree-node deletion using method: ", hclust_method)
    } else {
      dist_mat <- dist(X_plot)
    }
    
    dynamic_min_cluster_size <- max(5, ceiling(length(hc$labels) * 0.02))
    
    clone_assignments <- cutreeDynamic(
      dendro = hc,
      distM = as.matrix(dist_mat),
      method = "hybrid",
      deepSplit = dynamic_deep_split,
      minClusterSize = dynamic_min_cluster_size,
      pamRespectsDendro = FALSE
    )
    if (length(clone_assignments) != length(hc$labels)) {
      stop(
        "Dynamic Tree Cut returned ", length(clone_assignments),
        " assignments, but hclust has ", length(hc$labels),
        " labels for sample ", sample_id, "."
      )
    }
    names(clone_assignments) <- as.character(hc$labels)
    missing_assignment_barcodes <- setdiff(barcode_for_rows, names(clone_assignments))
    if (length(missing_assignment_barcodes) > 0) {
      stop(
        "Some true barcodes are missing from Dynamic Tree Cut assignments for sample ", sample_id,
        ". Examples: ", paste(head(missing_assignment_barcodes, 5), collapse = ", ")
      )
    }
    clone_assignments <- clone_assignments[barcode_for_rows]
    
    if (any(clone_assignments == 0)) {
      assigned_clusters <- sort(unique(clone_assignments[clone_assignments > 0]))
      if (length(assigned_clusters) == 0) {
        clone_assignments[] <- 1
      } else {
        dist_mat_full <- as.matrix(dist_mat)
        unassigned_barcodes <- names(clone_assignments)[clone_assignments == 0]
        for (barcode in unassigned_barcodes) {
          mean_dist_to_cluster <- sapply(
            assigned_clusters,
            function(cluster_id) {
              cluster_barcodes <- names(clone_assignments)[clone_assignments == cluster_id]
              mean(dist_mat_full[barcode, cluster_barcodes], na.rm = TRUE)
            }
          )
          clone_assignments[barcode] <- assigned_clusters[which.min(mean_dist_to_cluster)]
        }
      }
    }
    
    clone_assignments <- as.numeric(factor(clone_assignments, levels = sort(unique(clone_assignments))))
    names(clone_assignments) <- barcode_for_rows
    final_clone_assignments <- clone_assignments
    message(
      "Based on Dynamic Tree Cut, selected clone number: ", length(unique(clone_assignments)),
      " (deepSplit = ", dynamic_deep_split,
      ", minClusterSize = ", dynamic_min_cluster_size, ")"
    )
    
    # Construct the MergedNodes data frame exactly as in the previous version
    MergedNodes <- data.frame(
      Barcode = barcode_for_rows,
      Node    = paste0("Clone C", clone_assignments),
      stringsAsFactors = FALSE
    )

    # Remove Dynamic Tree Cut clones that contain too few metacells. Filtering
    # here keeps all downstream heatmaps, trees and spatial plots synchronized.
    clone_sizes_before_filter <- table(MergedNodes$Node)
    small_clone_names <- names(clone_sizes_before_filter)[
      clone_sizes_before_filter < min_clone_metacells
    ]

    if (length(small_clone_names) > 0) {
      deleted_small_clone_metacells <- MergedNodes %>%
        filter(Node %in% small_clone_names) %>%
        mutate(CloneSize = as.integer(clone_sizes_before_filter[Node]))

      write.csv(
        deleted_small_clone_metacells,
        paste0("./Deleted_small_clones_dynamic_tree_cut", analysis_suffix, ".csv"),
        row.names = FALSE
      )

      MergedNodes <- MergedNodes %>%
        filter(!Node %in% small_clone_names)

      message(
        "Removed ", length(small_clone_names), " clone(s) with fewer than ",
        min_clone_metacells, " metacells: ",
        paste(
          paste0(small_clone_names, " (n=", clone_sizes_before_filter[small_clone_names], ")"),
          collapse = ", "
        )
      )
    }

    if (nrow(MergedNodes) == 0) {
      stop(
        "All clones were removed because they contain fewer than ",
        min_clone_metacells, " metacells for sample ", sample_id, "."
      )
    }
    final_clone_assignments <- clone_assignments[MergedNodes$Barcode]
    
    # ==============================================================================
    # Step 3: retain the optional logic for manually merging selected clones
    # ==============================================================================
    combine_clone <- c() # Leave empty if no merging is needed; otherwise specify clones such as c("C4", "C5")
    
    if (length(combine_clone) > 1) {
      target_clone <- paste0("Clone ", combine_clone[which.min(as.numeric(gsub("C", "", combine_clone)))])
      
      MergedNodes <- MergedNodes %>%
        mutate(Node = ifelse(Node %in% paste0("Clone ", combine_clone), 
                             target_clone, 
                             Node))
    }
    
    MergedNodes$region = region[MergedNodes$Barcode,"region"]
    
    region_level = c("Normal", "Fibroadenoma", "ADH", "IP", "EPC", "DCIS", "IDC")
    clone_dominant_region <- MergedNodes %>%
      count(Node, region) %>%
      group_by(Node) %>%
      slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
      ungroup()
    clone_mapping <- clone_dominant_region %>%
      mutate(
        region_factor = factor(region, levels = region_level)
      ) %>%
      arrange(region_factor, desc(n)) %>%
      mutate(
        New_Node = paste0("Clone C", row_number())
      )
    message("--- Clone sorting and renaming mapping table ---")
    print(clone_mapping[, c("Node", "New_Node", "region", "n")])
    
    MergedNodes <- MergedNodes %>%
      left_join(clone_mapping %>% dplyr::select(Node, New_Node), by = "Node") %>%
      mutate(Node = New_Node) %>%
      dplyr::select(-New_Node) %>%
      mutate(Node = factor(Node, levels = paste0("Clone C", 1:nrow(clone_mapping))))
    
    MergedNodes_full = MergedNodes
    
    write.csv(
      MergedNodes,
      paste0("./Clone_info_full_dynamic_tree_cut", analysis_suffix, ".csv")
    )
  }else{
    setwd(paste0('E:/Projects/BRCA/results/ST/infercnv_metacell/',sample_id))
    
    load(paste0("E:/Projects/BRCA/results/ST/infercnv_metacell/", sample_id, "_subtrees.Rdata"))
    
    MergedNodes_full = read.csv(
      paste0("./Clone_info_full_dynamic_tree_cut", analysis_suffix, ".csv"),
      row.names = 1
    )
    
    MergedNodes <- MergedNodes_full %>%
      filter(!as.character(Node) %in% delete_clone_from_full_clones)
    
    remaining_clones <- unique(as.character(MergedNodes$Node))
    remaining_clones <- remaining_clones[order(as.numeric(gsub("Clone C", "", remaining_clones)))]
    
    new_clone_names <- paste0("Clone C", seq_along(remaining_clones))
    
    rename_dict <- setNames(new_clone_names, remaining_clones)
    
    # Print the dictionary to verify the renaming logic
    message("--- Clone renaming mapping dictionary after exclusion ---")
    print(rename_dict)
    
    # ==============================================================================
    # Step 3: apply the mapping dictionary and fix the new factor levels
    # ==============================================================================
    MergedNodes <- MergedNodes %>%
      mutate(
        Node = rename_dict[as.character(Node)],
        Node = factor(Node, levels = new_clone_names)
      )
  }
  
  alternative_clone_color = c("#fec02f","#8dcaef","#b889b9",
                              "#e377c2","#d25756","#b3dad3",
                              "#d4ae53","#cbd0dc","#bbcf5c","#f9dcc4","#ffef9f",
                              "#f48a3d")
  
  site.color <- c(
    "#1d8dba","#f4d47f","#7ac49d","#e76b41","#f6cade","#d3d3d3","#a0003e","#e7f5a0"
  )
  names(site.color) <- c("Normal","ADH","DCIS","IDC","IP","Undefined","Fibroadenoma","EPC")
  

  ##### Calculate the region composition for each subclone #####
  row_anno_clean = region
  row_anno_clean$region <- factor(row_anno_clean$region, levels = c("Normal","Fibroadenoma","ADH","IP","EPC","DCIS","IDC","Undefined"))
  row_anno_clean=row_anno_clean[MergedNodes$Barcode,]
  row_anno_clean$Clone <- MergedNodes$Node[match(rownames(row_anno_clean), MergedNodes$Barcode)]
  
  clone_composition <- table(row_anno_clean$region, row_anno_clean$Clone)
  clone_composition = t(clone_composition)
  
  clone_composition_df <- data.frame(matrix(as.vector(clone_composition), ncol = ncol(clone_composition)))
  rownames(clone_composition_df) <- rownames(clone_composition)
  colnames(clone_composition_df) <- colnames(clone_composition)
  
  clone_composition_df$node = rownames(clone_composition_df)
  pie_df <- nodepie(clone_composition_df, cols=1:8, alpha=0.9, color = site.color )
  
  # Load CNV data
  data1 = read.delim("./0.3.infercnv.observations.txt", sep = " ", row.names = 1)
  colnames(data1) <- gsub("\\.", "-", colnames(data1))
  cnv <- t(data1)
  data <- cnv[MergedNodes[grep("Clone", (MergedNodes$Node)),]$Barcode,]
  
  # Load chromatin annotations
  col_anno <- read.delim("./0.3.GeneOrderFile.txt", head = F)
  col_anno <- col_anno[col_anno$V1 %in% colnames(cnv),]
  col_anno <- col_anno[col_anno$V2 %in% paste0("chr", 1:22),]
  rownames(col_anno) <- col_anno[,1]
  
  col_anno_clean <- data.frame(col_anno[,1:2], row.names = rownames(col_anno))
  intersect_gene <- intersect(colnames(cnv), rownames(col_anno_clean))
  cnv_hdata_clean <- cnv[, intersect_gene]
  col_anno_clean <- col_anno_clean[rownames(col_anno_clean) %in% intersect_gene,]
  col_anno_clean <- data.frame(Chromosome = col_anno_clean[,2], row.names = rownames(col_anno_clean))
  col_anno_clean$Chromosome <- factor(col_anno_clean$Chromosome, levels = paste0("chr", c(1:22)))
  
  # Calculate chromatin boundaries
  tb = table(col_anno_clean);
  tb = tb[ paste0("chr", 1:22)]
  seg.chr = cumsum(tb)
  
  # Plot the heatmap corresponding to MergeNodes_full
  if(length(delete_clone_from_full_clones) == 0){
    row_anno_full = region
    row_anno_full$region <- factor(row_anno_full$region, levels = c("Normal","Fibroadenoma","ADH","IP","EPC","DCIS","IDC","Undefined"))
    row_anno_full=row_anno_full[MergedNodes_full$Barcode,]
    row_anno_full$Clone <- MergedNodes_full$Node[match(rownames(row_anno_full), MergedNodes_full$Barcode)]
    row_anno_full$meta_cell = NULL
    
    allnode = mytree$tip.label
    cluster = MergedNodes_full$Node
    
    names(cluster) = MergedNodes_full$Barcode
    
    groupInfo <- lapply(
      levels(as.factor(cluster)),
      function(cc, cluster){
        names(cluster)[cluster == cc]
      },
      cluster
    )
    
    names(groupInfo) <- names(table(cluster))
    
    clone.color.full = alternative_clone_color[1:length(unique(MergedNodes_full$Node))]
    names(clone.color.full) = paste0("Clone C", 1:length(unique(MergedNodes_full$Node)))
    
    # Rotate the original metacell tree for display only. Metacells belonging
    # to smaller clone numbers are requested first in the bottom-to-top order;
    # the tree topology and all branch lengths remain unchanged.
    full_tip_clone <- as.character(
      MergedNodes_full$Node[
        match(mytree$tip.label, MergedNodes_full$Barcode)
      ]
    )
    full_tip_clone_number <- get_clone_display_number(full_tip_clone)
    desired_full_tip_order <- mytree$tip.label[
      order(full_tip_clone_number, seq_along(mytree$tip.label))
    ]
    mytree_display <- ape::rotateConstr(mytree, desired_full_tip_order)

    MPT = groupOTU(mytree_display, groupInfo)
    gp1 <- ggtree(MPT, aes(color = group)) + theme_tree() + geom_tippoint(alpha = 0.7)
    gp1 <- gp1 + scale_color_manual(values=clone.color.full[grep("Clone", names(clone.color.full))])
    gp1 + geom_text(aes(label = node))
    dd = subset(gp1@data, isTip)
    spot_order1 <- dd$label[order(dd$y, decreasing=TRUE)]
    
    message("Plotting heatmaps for all clones")
    data_full = cnv[!(duplicated(rownames(cnv))),]
    data_full = data_full[MergedNodes_full$Barcode,]
    
    
    library(pheatmap)
    library(pals)
    pals::cols25()[1:22] -> chr.color
    names(chr.color) = paste0("chr", c(1:22))
    
    aCol = list(Chromosome = chr.color, region = site.color, Clone = na.omit(clone.color.full))
    
    data_clean <- data_full[,rownames(col_anno_clean)]
    
    #Generate CNV heatmap
    mid_chr <- cumsum(tb)			
    mid_tb <- round(tb/2)			
    mid_chr_mid <- round(mid_chr/2)
    mid_chr_mid[1] <- round(mid_chr[1]/2)
    mid_chr_mid[2:length(mid_chr_mid)] <- round(mid_chr[-1] - mid_tb[-1])
    
    data_clean_map <- data_clean
    
    for (k in 1:length(colnames(data_clean_map))){
      colnames(data_clean_map)[k] <- " "
    }
    colnames(data_clean_map)[mid_chr_mid] <- gsub("chr", "", names(seg.chr))
    
    mat <- data_clean_map[intersect(spot_order1, rownames(data_clean_map)), , drop = FALSE]
    
    max_abs <- max(abs(range(mat, na.rm = TRUE)))
    breaks  <- seq(-max_abs, max_abs, length.out = 101)              # 101 breakpoints define 100 color intervals
    colors  <- colorRampPalette(c("navy","white","firebrick"))(100)   # Use white at the midpoint
    
    pheatmap(
      mat,
      scale = "none",
      color = colors,
      breaks = breaks,                # Required: pass the breaks explicitly
      angle_col = "90",
      cluster_cols = FALSE,
      gaps_col = seg.chr,
      annotation_colors = aCol,
      annotation_row = row_anno_full[rownames(mat), , drop = FALSE],
      annotation_col = NA,
      fontsize = 8,
      cluster_rows = FALSE,
      show_colnames = TRUE,
      show_rownames = FALSE,
      border_color = NA,
      width = 10, height = 5,
      cutree_col = 2,
      filename = paste0("0.3.2.heatmap_full_dynamic_tree_cut", analysis_suffix, ".png"),
      legend_breaks = c(-max_abs, 0, max_abs),               # Optional: legend ticks explicitly show 0 at the midpoint
      legend_labels = sprintf(c("%.2f","0","%.2f"), -max_abs, max_abs)
    )
  }
  
  data = data[!(duplicated(rownames(data))),]
  MergedNodes = MergedNodes[!(duplicated(MergedNodes$Barcode)),]
  MergedNodes_clean <- MergedNodes[ MergedNodes$Barcode %in% rownames(data), ]
  rownames(MergedNodes_clean) <- MergedNodes_clean$Barcode
  all(rownames(data) == rownames(MergedNodes_clean))
  
  write.csv(
    MergedNodes,
    file = paste0('./0.3.2.metacell_clone_info_final_dynamic_tree_cut', analysis_suffix, '.csv')
  )
  
  data_merge <- data.frame(matrix(0, nrow = length(table(MergedNodes_clean$Node)), ncol = ncol(data)))
  rownames(data_merge) <- names(table(MergedNodes_clean$Node))
  colnames(data_merge) <- colnames(data)
  
  for (i in 1:ncol(data_merge)){
    data_merge[,i] <- tapply(data[,i], MergedNodes_clean$Node, mean)
  }
  data_merge <- data_merge[,rownames(col_anno_clean)]
  
  fig_thres <- c(max((data_merge)), min((data_merge)))
  fig_thres <- min(abs(fig_thres))
  data_merge[data_merge > 0 + fig_thres] <- 0 + fig_thres
  data_merge[data_merge < 0 - fig_thres] <- 0 - fig_thres
  
  #Load the infercnv output of reference
  data2 = read.delim("./0.3.infercnv.reference.txt", sep = " ", row.names = 1)
  colnames(data2) <- gsub("\\.", "-", colnames(data2))
  cnv_ref <- t(data2)
  #cnv_ref = cnv_ref[grepl(sample_id, rownames(cnv_ref)),]
  
  data_merge_ref <- rbind(data_merge, "Ref" = apply(cnv_ref, 2, mean))
  data_merge_ref_clean <- data_merge_ref
  data_merge_ref_clean[abs(data_merge_ref_clean) < log(1.02, 2)] <- 0

  # Optionally convert the noise-filtered CNV profile into three directional
  # states. sign() maps losses to -1, neutral bins to 0 and gains to 1.
  data_merge_ref_tree <- data_merge_ref_clean
  if (discretize_cnv_pattern) {
    data_merge_ref_tree <- sign(data_merge_ref_tree)
    message("Clone-level tree uses discretized CNV states: -1, 0 and 1.")
  } else {
    message("Clone-level tree uses continuous noise-filtered CNV values.")
  }
  
  data_dist <- dist(data_merge_ref_tree)
  
  if(length(unique(MergedNodes$Node))> 1){
    # Clone-level tree construction by neighbor joining (NJ) or maximum
    # parsimony (MP). MP treats every CNV bin as an unordered 3-state character.
    library(ggtree)
    library(phangorn)
    library(ggtree)

    if (tree_method == "NJ") {
      clone_tree <- phangorn::NJ(data_dist)
    } else {
      mp_character_matrix <- matrix(
        as.character(data_merge_ref_tree),
        nrow = nrow(data_merge_ref_tree),
        ncol = ncol(data_merge_ref_tree),
        dimnames = dimnames(data_merge_ref_tree)
      )
      mp_character_matrix[mp_character_matrix == "-1"] <- "loss"
      mp_character_matrix[mp_character_matrix == "0"] <- "neutral"
      mp_character_matrix[mp_character_matrix == "1"] <- "gain"
      mp_data <- phangorn::phyDat(
        mp_character_matrix,
        type = "USER",
        levels = c("loss", "neutral", "gain")
      )

      if (nrow(data_merge_ref_tree) == 3) {
        # With three tips there is only one unrooted binary topology. NJ is used
        # only to instantiate that topology; ACCTRAN below assigns MP lengths.
        clone_tree <- phangorn::NJ(data_dist)
      } else {
        mp_start_tree <- phangorn::NJ(data_dist)
        clone_tree <- phangorn::pratchet(
          mp_data,
          start = mp_start_tree,
          method = "fitch",
          trace = 0
        )
      }
      clone_tree <- phangorn::acctran(clone_tree, mp_data)
      attr(clone_tree, "parsimony_score") <- phangorn::parsimony(
        clone_tree,
        mp_data,
        method = "fitch"
      )
      message(
        "Maximum-parsimony score: ",
        attr(clone_tree, "parsimony_score")
      )
    }

    my_tree <- as.phylo(clone_tree)
    allnode = my_tree$tip.label
    cluster = allnode
    names(cluster) = allnode
    
    groupInfo <- lapply(levels(as.factor(cluster)), function(cc, cluster, allnode){
      return(names(cluster)[cluster == cc])
    },	cluster, allnode)					 
    names(groupInfo) <- names(table(cluster))
    
    #Specify the root clone as Ref.
    stopifnot("Ref" %in% clone_tree$tip.label)
    clone_tree <- ape::root(clone_tree, outgroup = "Ref", resolve.root = TRUE)
    clone_tree <- rotate_clone_tree_for_display(
      clone_tree,
      reference_label = "Ref"
    )
    message(
      "Requested clone-tree tip order from bottom to top: ",
      paste(
        attr(clone_tree, "desired_tip_order_bottom_to_top"),
        collapse = " -> "
      )
    )
    MPT5 = groupOTU(clone_tree, groupInfo)
    
    clone.color = alternative_clone_color[1:length(unique(MergedNodes$Node))]
    names(clone.color) = paste0("Clone C", 1:length(unique(MergedNodes$Node)))
    
    #Calculate the spot number
    dot_size = rep(1.95, length(clone.color))
    names(dot_size) = names(clone.color)
    dot_size = c(dot_size, "Ref" = 1.95)
    
    clone.color["Ref"] <- "darkgrey"
    gp5 <- ggtree(MPT5, aes(color = group), size = 1 / .pt) + theme_tree() 
    gp5 <- gp5 + scale_color_manual(values=clone.color[grep("Clone|NT", names(clone.color))])
    gp5 <- gp5 + geom_tippoint(size = dot_size)
    gp5$data$label <- gsub("Clone", "", gp5$data$label)
    
    gp5 <- gp5 + geom_tiplab(hjust = 0.55, vjust = 0.45, color = "black", size = 3.8 / .pt) +  
      theme(legend.position="none")
    
    clone_composition_df_ref <- rbind(clone_composition_df, "Ref" = 0)
    clone_composition_df_ref$node <- c(1:nrow(clone_composition_df_ref))
    
    pie_df <- nodepie(clone_composition_df_ref, cols=1:8, alpha=0.9, color = site.color )
    
    dis_node <- gp5$data[gp5$data$isTip, ]
    xmin <- min(dis_node$x)
    xmax <- max(dis_node$x)
    
    left_pad  <- 0.75                
    right_pad <- 0.60                
    x_range_data <- (xmax + right_pad) - (xmin - left_pad)
    
    num_leaves <- length(clone_tree$tip.label) 
    y_range_data <- num_leaves 
    plot_area_height_mm <- (num_leaves - 1) * 4.01 + 5.2
    margin_y_mm <- (10 + 10) / .pt 
    pdf_height_mm <- plot_area_height_mm + margin_y_mm
    pdf_height_inch <- pdf_height_mm / 25.4
    
    pdf_width_mm <- 35
    pdf_width_inch <- pdf_width_mm / 25.4
    
    margin_x_mm_single <- 1 
    margin_x_pt <- margin_x_mm_single * .pt # Approximately 2.84 pt
    
    plot_area_width_mm <- pdf_width_mm - (margin_x_mm_single * 2)
    
    pie_width_data  <- x_range_data * (3.5 / plot_area_width_mm)
    pie_height_data <- y_range_data * (3.5 / plot_area_height_mm)
    
    # ---------------------------------------------------------
    
    p <- gp5
    p <- inset(p, pie_df, 
               height = pie_height_data, 
               width  = pie_width_data,
               hjust  = c(-0.25 - xmax + dis_node$x), 
               vjust  = 0.05)   
    
    p <- p +
      xlim(xmin - left_pad, xmax + right_pad) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = unit(c(10, margin_x_pt, 10, margin_x_pt), "pt"))  
    
    pdf(
      paste0("0.3.4.clone_tree_main_dynamic_tree_cut", analysis_suffix, ".pdf"),
      width = pdf_width_inch,
      height = pdf_height_inch
    )
    print(p)
    dev.off()
    
    left_pad  <- 0.75                
    right_pad <- 0.60                
    x_range_data <- (xmax + right_pad) - (xmin - left_pad)
    
    num_leaves <- length(clone_tree$tip.label) 
    y_range_data <- num_leaves 
    plot_area_height_mm <- (num_leaves - 1) * 4.01 + 5.2
    margin_y_mm <- (10 + 10) / .pt 
    pdf_height_mm <- plot_area_height_mm + margin_y_mm
    pdf_height_inch <- pdf_height_mm / 25.4
    
    pdf_width_mm <- 24
    pdf_width_inch <- pdf_width_mm / 25.4
    
    margin_x_mm_single <- 1 
    margin_x_pt <- margin_x_mm_single * .pt # Approximately 2.84 pt
    
    plot_area_width_mm <- pdf_width_mm - (margin_x_mm_single * 2)
    
    pie_width_data  <- x_range_data * (3.5 / plot_area_width_mm)
    pie_height_data <- y_range_data * (3.5 / plot_area_height_mm)
    
    # ---------------------------------------------------------
    
    p <- gp5
    p <- inset(p, pie_df, 
               height = pie_height_data, 
               width  = pie_width_data,
               hjust  = c(-0.25 - xmax + dis_node$x), 
               vjust  = 0.05)   
    
    p <- p +
      xlim(xmin - left_pad, xmax + right_pad) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = unit(c(10, margin_x_pt, 10, margin_x_pt), "pt"))  
    
    pdf(
      paste0("0.3.4.clone_tree_supp_dynamic_tree_cut", analysis_suffix, ".pdf"),
      width = pdf_width_inch,
      height = pdf_height_inch
    )
    print(p)
    dev.off()
    
  }
  
  
  # Plot the heatmap
  if(length(unique(MergedNodes$Node)) > 1){
    d5 = fortify(MPT5)
    dd = subset(d5, isTip)
    # ComplexHeatmap draws its first row at the top. Reading the rotated tree
    # from high y to low y therefore keeps smaller clone numbers near the
    # bottom of both the tree and heatmap without changing tree topology.
    spot_order <- as.character(dd$label[order(dd$y, decreasing = TRUE)])
  }else{
    spot_order = c("Ref", "Clone C1")
  }
  
  row_anno_clean_merge <- clone_composition_df[,1:2]
  colnames(row_anno_clean_merge)[1:2] <- "Clone"
  row_anno_clean_merge[,1:2] <- rownames(row_anno_clean_merge)
  row_anno_clean_merge <- as.data.frame(row_anno_clean_merge[,1])
  colnames(row_anno_clean_merge) <- "Clone"
  rownames(row_anno_clean_merge) <- row_anno_clean_merge[,1]
  row_anno_clean_merge[nrow(row_anno_clean_merge)+1,] <- "Ref"
  rownames(row_anno_clean_merge)[nrow(row_anno_clean_merge)] <- "Ref"

  # ComplexHeatmap aligns a standalone row annotation by position. Reorder the
  # annotation explicitly so that each color bar row matches the CNV matrix row.
  missing_heatmap_annotations <- setdiff(spot_order, rownames(row_anno_clean_merge))
  if (length(missing_heatmap_annotations) > 0) {
    stop(
      "Heatmap clone annotation is missing: ",
      paste(missing_heatmap_annotations, collapse = ", ")
    )
  }
  row_anno_clean_merge <- row_anno_clean_merge[spot_order, , drop = FALSE]

  data_merge <- data_merge[,rownames(col_anno_clean)]
  
  # Keep CNV amplitudes in the heatmap. Its row order still follows the clone tree,
  # which may have been inferred from either continuous or three-state profiles.
  data_merge_map <- data_merge_ref_clean
  colnames(data_merge_map)[1:length(data_merge_map)] <- " "
  colnames(data_merge_map)[mid_chr_mid] <- gsub("chr", "", names(seg.chr))
  #data_merge_map[abs(log(data_merge_map, 2)) < log(1.02, 2)] <- 1
  
  clone.color = alternative_clone_color[1:length(unique(MergedNodes$Node))]
  names(clone.color) = paste0("Clone C", 1:length(unique(MergedNodes$Node)))
  aCol = list(Chromosome = chr.color, region = site.color, Clone = na.omit(clone.color))
  
  aCol2 <- aCol
  aCol2$Clone[names(aCol2$Clone) == "Ref"] <- "darkgrey"
  aCol2$Clone <- aCol2$Clone[grep("Clone|Ref", names(aCol2$Clone))]
  aCol2$Clone <- c(aCol2$Clone,'darkgrey')
  names(aCol2$Clone)[length(aCol2$Clone)]='Ref'
  missing_heatmap_colors <- setdiff(spot_order, names(aCol2$Clone))
  if (length(missing_heatmap_colors) > 0) {
    stop(
      "Heatmap clone color is missing: ",
      paste(missing_heatmap_colors, collapse = ", ")
    )
  }
  heatmap_clone_colors <- aCol2$Clone[spot_order]
  
  library(ComplexHeatmap)
  library(circlize)
  
  # --- 1. Convert pheatmap gap indices (seg.chr) to a grouping vector ---
  num_cols <- ncol(data_merge_map)
  block_lengths <- diff(c(0, seg.chr, num_cols))
  col_split_factor <- rep(1:length(block_lengths), times = block_lengths)
  
  val_min <- min(data_merge_map[spot_order,], na.rm = TRUE)
  val_max <- max(data_merge_map[spot_order,], na.rm = TRUE)
  col_fun <- colorRamp2(
    c(val_min, 0, val_max),
    c("navy", "white", "firebrick")
  )

  row_ha <- rowAnnotation(
    Clone = row_anno_clean_merge$Clone,
    col = list(Clone = heatmap_clone_colors),
    show_annotation_name = FALSE, # Hide the "Clone" title above the annotation bar
    simple_anno_size = unit(1.32, "mm")
  )
  
  ht_opt(ROW_ANNO_PADDING = unit(0, "mm"))
  
  ht <- Heatmap(
    as.matrix(data_merge_map[spot_order, ]), 
    name = "CNV",
    col = col_fun,
    
    width = unit(53.8, "mm"),
    height = nrow(data_merge_map) * unit(4.01, "mm"),
    column_split = col_split_factor,
    column_gap = unit(0.66, "mm"),
    column_title = NULL,
    
    # --- Clustering and display options ---
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    # --- Border settings ---
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.25),
    
    rect_gp = gpar(col = NA),
    
    # --- Annotation and alignment ---
    left_annotation = row_ha
  )
  
  pdf_width_inch  <- (52.5 + 50) / 25.4
  pdf_height_inch <- (nrow(data_merge_map) * 4.01 + 15) / 25.4 
  
  pdf(
    paste0("./0.3.5.heatmap_subclone_final_Complex_main_dynamic_tree_cut", analysis_suffix, ".pdf"),
    width = pdf_width_inch,
    height = pdf_height_inch
  )
  # Render the heatmap with a physical outer margin of 2 mm on all sides
  ComplexHeatmap::draw(ht, padding = unit(c(2, 2, 2, 2), "mm")) 
  dev.off()
  
  
  # Plot the supplementary heatmap
  row_ha <- rowAnnotation(
    Clone = row_anno_clean_merge$Clone,
    col = list(Clone = heatmap_clone_colors),
    show_annotation_name = FALSE, # Hide the "Clone" title above the annotation bar
    simple_anno_size = unit(0.98, "mm")
  )
  
  ht <- Heatmap(
    as.matrix(data_merge_map[spot_order, ]), 
    name = "CNV",
    col = col_fun,
    
    width = unit(39.02, "mm"),
    height = nrow(data_merge_map) * unit(4.01, "mm"),
    column_split = col_split_factor,
    column_gap = unit(0.49, "mm"),
    column_title = NULL,
    
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.25),
    
    rect_gp = gpar(col = NA),
    
    left_annotation = row_ha
  )
  
  pdf_width_inch  <- (50 + 40) / 25.4
  pdf_height_inch <- (nrow(data_merge_map) * 4.01 + 15) / 25.4 
  
  pdf(
    paste0("./0.3.5.heatmap_subclone_final_Complex_supp_dynamic_tree_cut", analysis_suffix, ".pdf"),
    width = pdf_width_inch,
    height = pdf_height_inch
  )

  ComplexHeatmap::draw(ht, padding = unit(c(2, 2, 2, 2), "mm")) 
  dev.off()
  

  stobj = readRDS(paste0("E:/Projects/BRCA/Intermediate_Data/ST/object_unified/", sample_id, "_seg.rds"))

  all_members = read.csv("E:/Projects/BRCA/results/ST/SuperCell/all_metacell_members.csv", row.names = 1)
  sp_member = filter(all_members, sample == sample_id)
  sp_member = filter(sp_member, meta_cell %in% MergedNodes_clean$Barcode)
  sp_member$Clone = as.character(MergedNodes_clean$Node[match(sp_member$meta_cell, MergedNodes_clean$Barcode)])
  
  stobj$cnv_clone = "non_epi"
  stobj$cnv_clone[sp_member$cellid] = sp_member$Clone
  
  stobj$cnv_clone = factor(stobj$cnv_clone, levels = c(paste0("Clone C", 1:(length(unique(stobj$cnv_clone))-1)),
                                                       "non_epi"))
  color_spatial = aCol2$Clone
  names(color_spatial)[length(color_spatial)] = "non_epi"
  color_spatial["non_epi"] = "#d3d3d3"
  plt = SpatialDimPlot(stobj ,group.by = "cnv_clone", image.alpha = 0, pt.size.factor = seurat_spot_size[sample_id]+2) + 
    scale_fill_manual(values = color_spatial) + 
    theme(plot.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = "NA", size = 2))
  ggsave(
    paste0("./clone_spatial_final_dynamic_tree_cut", analysis_suffix, ".pdf"),
    plot = plt, dpi = 700, width = 20, height = 17, units = "cm"
  )
  ggsave(
    paste0("./clone_spatial_final_dynamic_tree_cut", analysis_suffix, ".png"),
    plot = plt, dpi = 700, width = 20, height = 17, units = "cm"
  )
}
