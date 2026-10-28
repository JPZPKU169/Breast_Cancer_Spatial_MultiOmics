# Run WNN in a loop
sample_list_sm_matched = c(
  "P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P14","P15","P16","P26","P27","P30","P32","P33","P34","P35","P36","P37","P38","P39",
  "P40_1","P40_2","P42","P43_1","P43_2","P44","P45","P46","P47","P48","P49"
)

library(future)
plan(multisession, workers = 8)
options(future.globals.maxSize = 30 * 1024^3)

for(spID in sample_list_sm_matched[1]){
  print(paste0("processing ", spID))
  sm_obj = loadH5(FileName = paste0("E:/Projects/BRCA/Intermediate_Data/SM5u/matched_object_epi/", spID,"_match_st_epi.h5"),
                  assay = "Spatial", image_name = "slice1")
  spID_st = gsub("SM","ST", spID)
  st_obj = readRDS(file = paste0("E:/Projects/BRCA/Intermediate_Data/ST/object_unified/", spID_st, "_seg.rds"))
  st_obj = st_obj[,st_obj$celltype == "Epithelial"]
  
  epi_intersect = intersect(colnames(st_obj), colnames(sm_obj))
  print(paste0(as.character(length(epi_intersect)), " common cells are found"))
  sm_obj = sm_obj[,epi_intersect]
  st_obj = st_obj[,epi_intersect]
  if(all(colnames(st_obj@assays$Spatial) == colnames(sm_obj@assays$Spatial))){
    print("cell ids are in the same order")
  }
  st_obj@assays$SM = sm_obj@assays$Spatial
  # Spatial transcriptomics
  DefaultAssay(st_obj) = "Spatial"
  st_obj = NormalizeData(st_obj)
  st_obj = FindVariableFeatures(st_obj)
  st_obj = ScaleData(st_obj)
  st_obj = RunPCA(st_obj)
  # Spatial metabolomics
  DefaultAssay(st_obj) = "SM"
  st_obj = NormalizeData(st_obj) %>% FindVariableFeatures(nfeatures = 500) %>% ScaleData() %>% RunPCA(reduction.name = "sm.pca")
  
  # WNN
  st_obj <- FindMultiModalNeighbors(
    st_obj, reduction.list = list("pca", "sm.pca"), 
    dims.list = list(1:20, 1:20), modality.weight.name = "RNA.weight"
  )
  st_obj <- RunUMAP(st_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  st_obj <- FindClusters(st_obj, graph.name = "wsnn", algorithm = 3, resolution = 0.1, verbose = TRUE)
  st_obj <- FindClusters(st_obj, graph.name = "wsnn", algorithm = 3, resolution = 0.05, verbose = TRUE)
  
  saveRDS(st_obj, file = paste0("E:/Projects/BRCA/Intermediate_Data/ST/WNN_objects/", spID_st, "_wnn.rds"))
  
  
  output.dir = paste0("E:/Projects/BRCA/results/ST/WNN/", spID_st, "/")
  if(!dir.exists(output.dir)){
    dir.create(output.dir)
  }
  
  alternative_color_celltype = c("#f4a261","#2a9d8f","#a8dadc","#f5cac3","#ef476f",
                                 "#20854ECC","#9d4edd","#adc178","#E64B35CC","#5e60ce",
                                 "#B09C85CC","#DC0000CC","#8491B4CC","#F39B7FCC","#4DBBD5CC",
                                 "#A20056B2","#FFDC91CC","#5a189a","#3c6e71","#eaf4f4",
                                 "#ffd60a","#b7e4c7","#ff7b00","#001233","#9f86c0",
                                 "#ccdbfd","#41ead4","#495867","#004e98","#603808",
                                 "#ff5d8f","#70e000")
  alternative_color_celltype = c(alternative_color_celltype,alternative_color_celltype)
  
  cluster_name = names(table(st_obj$wsnn_res.0.1))[table(st_obj$wsnn_res.0.1) > 3]
  
  DimPlot(subset(st_obj, wsnn_res.0.1 %in% cluster_name), reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "wsnn_res.0.1") + 
    scale_color_manual(values = alternative_color_celltype)
  ggsave(paste0(output.dir, spID_st, "_wnn_umap_cluster_res.0.1.png"), width = 18, height = 15, units = "cm",
         dpi = 300)
  
  plt = SpatialDimPlot(subset(st_obj, wsnn_res.0.1 %in% cluster_name), group.by = "wsnn_res.0.1", image.alpha = 0,
                       pt.size.factor = seurat_spot_size[spID_st]) +
    scale_fill_manual(values = alternative_color_celltype) +
    theme(plot.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = "NA", size = 2))
  ggsave(paste0(output.dir, spID_st, "_wnn_cluster_spatial_res.0.1.png"), plot = plt,
         width = 20, height = 20, units = "cm", dpi = 700)
  
  cluster_name = names(table(st_obj$wsnn_res.0.05))[table(st_obj$wsnn_res.0.05) > 3]
  
  DimPlot(subset(st_obj, wsnn_res.0.05 %in% cluster_name), reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "wsnn_res.0.05") + 
    scale_color_manual(values = alternative_color_celltype)
  ggsave(paste0(output.dir, spID_st, "_wnn_umap_cluster_res.0.05.png"), width = 18, height = 15, units = "cm",
         dpi = 300)
  
  plt = SpatialDimPlot(subset(st_obj, wsnn_res.0.05 %in% cluster_name), group.by = "wsnn_res.0.05", image.alpha = 0,
                       pt.size.factor = seurat_spot_size[spID_st]) +
    scale_fill_manual(values = alternative_color_celltype) +
    theme(plot.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = "NA", size = 2))
  
  ggsave(paste0(output.dir, spID_st, "_wnn_cluster_spatial_res.0.05.png"), plot = plt,
         width = 20, height = 20, units = "cm", dpi = 700)
  
  st_obj@assays$Spatial@layers$scale.data = NULL
  st_obj@assays$SM@layers$scale.data = NULL
  saveH5(st_obj, FileName = paste0("E:/Projects/BRCA/Intermediate_Data/ST/WNN_objects/", spID_st, "_epi_wnn.h5"), 
         image_name = "slice1",  assay = "Spatial")
}


##### Define regions for P13 using transcriptomics alone #####
st13 = readRDS(file = paste0("E:/Projects/BRCA/Intermediate_Data/ST/object_unified/ST1088776_seg.rds"))
st13 = subset(st13, celltype == "Epithelial")

# Spatial transcriptomics
DefaultAssay(st13) = "Spatial"
st13 = NormalizeData(st13)
st13 = FindVariableFeatures(st13)
st13 = ScaleData(st13)
st13 = RunPCA(st13)
ElbowPlot(st13, ndims = 50)
st13 = FindNeighbors(st13, dims = 1:30)
st13 <- FindClusters(st13, resolution = 0.5)
st13 <- RunUMAP(st13, dims = 1:30)
DimPlot(st13, label = T)

st13@graphs
st13@assays$Spatial@layers$scale.data = NULL
saveH5(st13, FileName = paste0("E:/Projects/BRCA/Intermediate_Data/ST/WNN_objects/P13_epi_wnn.h5"), 
       image_name = "slice1",  assay = "Spatial")

