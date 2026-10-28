library(CellChat)
library(patchwork)
library(presto)
library(Seurat)

##### CellChat #####
table(brca$sub_celltype)
brca$cellchat_celltype = brca$sub_celltype
brca$cellchat_celltype[brca$cellchat_celltype %in% endo_celltype_levels] = "Endothelial"
brca$cellchat_celltype[brca$cellchat_celltype %in% c("CD4_Naive","CD4_Tcm","CD4_CCL5+_Tm","CD4_KLRB1+_Tem",
                                                     "CD4_CXCL13+_Tfh")] = "CD4T"
brca$cellchat_celltype[brca$cellchat_celltype %in% c("CD8_ZNF683+_CXCR6+_Trm","CD8_ITGA1+_Trm","CD8_ZNF683+_CXCR6-_Tm",
                                                     "CD8_GZMK+_Tem","CD8_Temra","CD8_Effector","CD8_ISG15+","CD8_Tex",
                                                     "CD8_KLRB1+_MAIT")] = "CD8T"
brca$cellchat_celltype[brca$cellchat_celltype %in% c("CD4_IKZF2+_Treg","CD4_TNFRSF9+_Treg")] = "Treg"
brca$cellchat_celltype[brca$cellchat_celltype %in% c("NK_CD56hi_CD16-","NK_CD56dim_CD16+","NK_TCF7+")] = "NK"
brca$cellchat_celltype[brca$cellchat_celltype %in% c("Early_PC_MS4A1low","PC_IGHG","PC_IGHA")] = "Plasma"

brca$cellchat_celltype = brca$sub_celltype
brca$cellchat_celltype[brca$sub_celltype == "Epithelial" & !(colnames(brca) %in% colnames(cancer))] = "non-malignant"

brca = brca[,!(brca$cellchat_celltype %in% c("Fibro_MKI67+","Macro_C8_MKI67+","T_MKI67+","non-malignant","Endo_MKI67+"))]
table(brca$cellchat_celltype)


data = GetAssayData(brca, assay = "RNA", layer = "data")
meta = brca@meta.data
colnames(meta)
colnames(meta)[7] = "samples"
meta$samples = factor(meta$samples)

# make object
brca_cellchat = createCellChat(object = data, meta = meta, group.by = "cellchat_celltype")

# select ligand-receptor interaction database
CellChatDB <- CellChatDB.human
CellChatDB.use = CellChatDB

brca_cellchat@DB <- CellChatDB.use
brca_cellchat = subsetData(brca_cellchat) # subset the expression data of signaling genes+
future::plan("multisession", workers = 4)
brca_cellchat = identifyOverExpressedGenes(brca_cellchat)
brca_cellchat <- identifyOverExpressedInteractions(brca_cellchat)
brca_cellchat <- projectData(brca_cellchat, PPI.human)
options(future.globals.maxSize = 50*1024^3)
brca_cellchat = computeCommunProb(brca_cellchat, raw.use = FALSE, type = "triMean") 
brca_cellchat = filterCommunication(brca_cellchat, min.cells = 10)

# extract inferred communication network
df.net.brca = subsetCommunication(brca_cellchat)

# Infer the cell-cell communication at a signaling pathway level
brca_cellchat <- computeCommunProbPathway(brca_cellchat)

levels(brca_cellchat@idents)
brca_cellchat <- aggregateNet(brca_cellchat)

