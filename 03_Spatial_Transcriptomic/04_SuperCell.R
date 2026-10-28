library(SuperCell)

###### Supercell ######
gamma = 12
k.knn = 5

for(spID in sample_list_st){

  object = readRDS(file = paste0("E:/Projects/BRCA/Intermediate_Data/ST/object_unified/", spID, "_seg.rds"))
  swed_info = read.csv(paste0("E:/Projects/BRCA/results/ST/WNN/", spID, "/swed_", spID, ".csv"))
  colnames(swed_info)[1] = "cell"
  object = subset(object, cells = swed_info$cell)
  
  # Plot the SWED spatial distribution using Seurat
  object@meta.data$swed = swed_info$swed[match(colnames(object), swed_info$cell)]
  alternative_color_celltype = c("#f4a261","#2a9d8f","#a8dadc","#f5cac3","#ef476f",
                                 "#20854ECC","#9d4edd","#adc178","#E64B35CC","#5e60ce",
                                 "#B09C85CC","#DC0000CC","#8491B4CC","#F39B7FCC","#4DBBD5CC",
                                 "#A20056B2","#FFDC91CC","#5a189a","#3c6e71","#eaf4f4",
                                 "#ffd60a","#b7e4c7","#ff7b00","#001233","#9f86c0",
                                 "#ccdbfd","#41ead4","#495867","#004e98","#603808",
                                 "#ff5d8f","#70e000")
  SpatialDimPlot(subset(object, swed != "others"), group.by = "swed", image.alpha = 0, pt.size.factor = seurat_spot_size[spID]) +
    scale_fill_manual(values = alternative_color_celltype)
  ggsave(paste0("E:/Projects/BRCA/results/ST/WNN/", spID, "/swed_spatial_dimplot_", spID, ".png"),
         width = 20, height = 17, units = "cm")
  
  # Supercell
  gamma = 12
  k.knn = 5
  
  object = FindVariableFeatures(object)
  swed_unique = unique(object$swed)
  swed_unique = swed_unique[str_detect(swed_unique, "DCIS|IDC|ADH|IP|EPC|Normal")]
  genes_used = union(intersect(model_genes, rownames(object)), VariableFeatures(object))
  
  for(sp_swed in swed_unique){
    cell_index = colnames(object)[object$swed == sp_swed]
    sub_obj = subset(object, cells = cell_index)
    sub_obj = ScaleData(sub_obj, features = genes_used)
    sub_obj = RunPCA(sub_obj, features = genes_used)
    matrix = sub_obj@reductions$pca@cell.embeddings
    
    simp = SCimplify_from_embedding(
      X =  matrix,
      k.knn = k.knn,
      gamma = gamma
    )
    
    ge = supercell_GE(as.matrix(GetAssayData(sub_obj, assay = "Spatial", layer = "counts")),
                      simp$membership)
    
    member_info = data.frame(cellid = colnames(sub_obj), meta_cell = paste0(spID, "_", sp_swed, "_", simp$membership))
    colnames(ge) = paste0(spID, "_", sp_swed, "_", 1:ncol(ge))
    
    output.dir = paste0("E:/Projects/BRCA/results/ST/SuperCell/", spID, "/")
    if(!(dir.exists(output.dir))){
      dir.create(output.dir)
    }
    
    write.csv(ge, file = paste0(output.dir, spID, "_", sp_swed, "_expr.csv"), sep = ",")
    write.csv(member_info, file = paste0(output.dir, spID, "_", sp_swed, "_member.csv"), sep = ",")
  }
}

metacell_matrix = matrix()
for(spID in sample_list_st){
  print(spID)
  swed_info = read.csv(paste0("E:/Projects/BRCA/results/ST/WNN/", spID, "/swed_", spID, ".csv"))
  colnames(swed_info)[1] = "cell"
  
  swed_unique = unique(swed_info$swed)
  swed_unique = swed_unique[str_detect(swed_unique, "DCIS|IDC|ADH|IP|EPC|Normal")]
  
  output.dir = paste0("E:/Projects/BRCA/results/ST/SuperCell/", spID, "/")
  for(sp_swed in swed_unique){
    print(sp_swed)
    sp_matrix = read.csv(paste0(output.dir, spID, "_", sp_swed, "_expr.csv"), row.names = 1)
    metacell_matrix = cbind(metacell_matrix, sp_matrix)
  }
}

metacell_matrix = metacell_matrix[,-1]
write.csv(metacell_matrix, "E:/Projects/BRCA/results/ST/SuperCell/all_metacell_matrix.csv")


##### Assign metacell regions by voting #####
metacell_region_all = data.frame(meta_cell = c(), region = c())
for(spID in sample_list_st){
  cell_region = read.csv(paste0("E:/Projects/BRCA/results/ST/cluster_ljy/all_sample_region/", spID, "_region.csv"))
  member_files_all = list.files(path = paste0("E:/Projects/BRCA/results/ST/SuperCell/",spID,"/"), 
                                pattern = "_member\\.csv$", full.names = TRUE)
  sp_member_merged = data.frame()
  for(fname in member_files_all){
    swed_member = read.csv(fname, row.names = 1)
    if(nrow(sp_member_merged) == 0){
      sp_member_merged = swed_member
    }else{
      sp_member_merged = rbind(sp_member_merged, swed_member)
    }
  }
  write.csv(sp_member_merged, file = paste0("E:/Projects/BRCA/results/ST/SuperCell/",spID,"/", spID, "_metacell_members.csv"))
  cell_region = cell_region[cell_region$X %in% sp_member_merged$cellid,]
  merged_data <- sp_member_merged %>%
    left_join(cell_region, by = c("cellid" = "X"))
  merged_data$region[is.na(merged_data$region)] = "undefined"
  region_counts <- merged_data %>%
    group_by(meta_cell, region) %>%
    summarise(count = n(), .groups = 'drop')
  max_regions <- region_counts %>%
    group_by(meta_cell) %>%
    mutate(max_count = max(count)) %>%
    filter(count == max_count) %>%

    mutate(
      region_count = n(),
      final_region = ifelse(region_count > 1, 'undefined', region)
    ) %>%

    distinct(meta_cell, .keep_all = TRUE) %>%
    select(meta_cell, final_region)
  names(max_regions) <- c("meta_cell", "region")
  write.csv(max_regions, file = paste0("E:/Projects/BRCA/results/ST/SuperCell/", spID, "_metacell_region.csv"))
  metacell_region_all = rbind(metacell_region_all, max_regions)
}
write.csv(metacell_region_all, file = "E:/Projects/BRCA/results/ST/SuperCell/metacell_region_all.csv")

