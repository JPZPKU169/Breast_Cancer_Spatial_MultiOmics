library(Seurat)

##### Integrate all myeloid cells #####
for(spID in sample_list_st){
  print(spID)
  object = readRDS(paste0("E:/Projects/BRCA/Intermediate_Data/ST/object_unified/",spID,"_seg.rds"))
  sp_mye = subset(object, celltype %in% c("DC","Macrophages","Mast"))
  
  if(!exists("mye_st")){
    mye_st = sp_mye
  }else{
    mye_st = merge(mye_st, sp_mye)
  }
}
layers_to_remove = names(mye_st@assays$Spatial@layers)
layers_to_remove = layers_to_remove[grepl("scale", layers_to_remove)]
for(sp_layer in layers_to_remove){
  mye_st@assays$Spatial@layers[[sp_layer]] = NULL
}

mye_st = JoinLayers(mye_st)
saveRDS(mye_st, file = "E:/Projects/BRCA/Intermediate_Data/ST/macro_objects/st_mye_before_integration.rds")

# Highly variable genes in myeloid cells
mye = readRDS(file = "E:/Projects/BRCA/Intermediate_Data/SC/SC_mye_object.rds")
mye = FindVariableFeatures(mye, nfeatures = 1500)

variable.features = VariableFeatures(mye)
variable.features = variable.features[variable.features %in% rownames(mye_st)]
saveRDS(variable.features, file = "E:/Projects/BRCA/Intermediate_Data/ST/macro_objects/mye_variable_features.rds")

# Server-side code
plan(multisession, workers = 16)
options(future.globals.maxSize = 80 * 1024^3)

mye = readRDS(file = "/io01dat03/PROJECT/zoujiping/brca/st/b_object_unified/st_mye_before_integration.rds")
variable.features = readRDS(file = "/io01dat03/PROJECT/zoujiping/brca/st/a_macro/mye_variable_features.rds")

mye <- ScaleData(mye, features = variable.features)
mye = RunPCA(mye, features = variable.features, npcs = 50)

mye[["Spatial"]] = split(mye[["Spatial"]], f = mye$sample)

mye = IntegrateLayers(
  object = mye, method = HarmonyIntegration,
  orig.reduction = "pca",
  features = variable.features,
  new.reduction = "harmony",theta = 2.5,
  verbose = TRUE
)

mye = JoinLayers(mye)

mye = FindNeighbors(mye, reduction = "harmony", dims = 1:30)
mye = FindClusters(mye, resolution = 0.8)
mye = RunUMAP(mye, reduction = "harmony", dims = 1:30)

saveRDS(mye, file = "/io01dat03/PROJECT/zoujiping/brca/st/b_object_unified/st_mye_integrated.rds")


