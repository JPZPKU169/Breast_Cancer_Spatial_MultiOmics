library(Seurat)
library(tidyverse)
library(ggplot2)

markers = c("EPCAM","AMBP","MUC1","KRT8","KRT19","FOXA1","GATA3","CD24","PROM1", # Luminal
            "MGP","SCGB2A2","SLPI","LTF","PTN","KIT","ALDH1A3","GABRP", # Luminal Progenitor
            "MUCL1","ERBB4","ANKRD30A","AZGP1","AGR2","STC2", # Luminal Mature
            "KRT5","KRT17","KRT14","MMP7","DST","MIR205HG","MT1X","OXTR","FST","MME","ITGA6", # Myoepithelial
            "ESR1","PGR","ERBB2", # Breast Cancer 
            "COL1A1","COL1A2","COL3A1", "DCN","ACTA2","SFRP2","PDGFRA","CDH11","ITGBL1", # Fibroblasts
            "MYH11","ACTG2","TAGLN", # SMC
            "MCAM","PDGFRB", # Pericytes
            "TRAC","TRBC2","CD3D","CD3E","CD8A","GZMB","GNLY", # T
            "NCAM1","NKG7", "KLRD1","KLRF1",
            'CD68','CD80','CD86','CD14','APOE','LYZ',"AIF1","CSF1R","LILRB4","CD163","LAMP2","LAMP3", #
            "PLD4","LILRA4","SPIB","IRF4", # pDC
            "CXCL8","G0S2","CSF3R","FCGR3B", # Neutrophil
            'MS4A1','CD79A','CD79B',"BLNK",'CD19','CXCR5','CD22',
            'IGHG1','IGKC','JCHAIN','TNFRSF17',"CD38","CD27",
            "TPSB2","TPSAB1","CPA3","MS4A2","HPGDS", # Mast
            "PNPLA2","CAV1","FABP4","PPARG","CEBPA","LEP", # Adipocyte
            "PECAM1","ENG","CD34","VWF","CDH5","SPARCL1" # Endothelial
)

sample_list_sc = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13","P14","P15","P17","P18","P19","P20",
                   "P21","P22","P23","P24","P25","P26","P27","P28","P29","P30","P31","P32","P37","P39","P40","P41","P42",
                   "P48","P49","P50","P51"
)


brca = NormalizeData(brca)
brca = FindVariableFeatures(brca)
brca = ScaleData(brca)
brca = RunPCA(brca)


immunogloblin_and_tcr = read.table(file = "D:/references/BioMart/Immunogloblin_and_TCR_genes.txt", sep = "\t", header = T)


##### Cluster directly without batch correction #####
brca = JoinLayers(brca)
brca = NormalizeData(brca)
brca = FindVariableFeatures(brca)
variable.features = VariableFeatures(brca)
if(length(grep("^RP([0-9]+-|[LS])", variable.features)) > 0){
  variable.features = variable.features[-grep("^RP([0-9]+-|[LS])", variable.features)]
}
if("MALAT1" %in% variable.features){
  variable.features = variable.features[!(variable.features == "MALAT1")]
}
if(any(variable.features %in% immunogloblin_and_tcr$Gene.name)){
  table(variable.features %in% immunogloblin_and_tcr)
  variable.features = variable.features[!(variable.features %in% immunogloblin_and_tcr$Gene.name)]
}
VariableFeatures(brca) = variable.features

brca = ScaleData(brca, vars.to.regress = c("nCount_RNA","percent.mt"))
brca = RunPCA(brca, npcs = 50)
ElbowPlot(brca, ndims = 50)

##### harmony #####
brca = IntegrateLayers(
  object = brca, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE, lambda = 3
)

ElbowPlot(brca, ndims = 50, reduction = "pca")

brca = FindNeighbors(brca, reduction = "harmony", dims = 1:35)
brca = FindClusters(brca, resolution = 1.2, cluster.name = "seurat_cluster_harmony")
brca = RunUMAP(brca, reduction = "harmony", dims = 1:35, reduction.name = "umap_harmony")
brca = RunTSNE(brca, reduction = "harmony", dims = 1:35, reduction.name = "tsne_harmony")
DimPlot(brca, reduction = "umap_harmony", label = T)
DimPlot(brca, reduction = "tsne_harmony", label = T)
VlnPlot(brca, markers, group.by = "seurat_cluster_harmony",stack = T, pt.size = 0)

clusteavr_markers = FindAllMarkers(brca)
top_cluster_markers = cluster_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
markers  = FindMarkers(brca, ident.1 = 39, min.pct = 0.1)


new.cluster.ids = c("T","T","Endothelial","Mysenchymal","Mysenchymal","Myeloid", # 0-5
                    "Epithelial","Epithelial","T","Epithelial","Mysenchymal", # 6-10
                    "B","Basal","Epithelial","Epithelial","Epithelial", # 11-15
                    "Epithelial","Epithelial","Epithelial","Epithelial","Epithelial", # 16-20
                    "Epithelial","T","T","Epithelial","T", # 21-25
                    "Epithelial","Epithelial","T","Plasma","Epithelial", # 26-30
                    "Epithelial","Epithelial","Myeloid","Myeloid","T"
)


Idents(brca) = brca@meta.data$seurat_clusters
length(new.cluster.ids)
table(new.cluster.ids)
head(Idents(brca))
names(new.cluster.ids) <- levels(brca)
brca <- RenameIdents(brca, new.cluster.ids)
DimPlot(brca, label = T, raster = FALSE)
brca$celltype = as.character(Idents(brca))

celltype_levels = c("Epithelial","Basal", "Mysenchymal", "Endothelial", "T","B","Plasma","Myeloid")
all(celltype_levels %in% unique(brca$celltype))
all(unique(brca$celltype) %in% celltype_levels)
brca$celltype = factor(brca$celltype, levels = celltype_levels)
DimPlot(brca, label = T, group.by = "celltype")


##### Plot the overall UMAP and DotPlot #####
library('ggsci')
library('ggplot2')
library("scales")
library(Scillus)
library(ggrepel)
library(tidyverse)
library(scales)
library(ggrastr)


int_umap = brca@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type = brca@meta.data$celltype) %>%
  cbind(SampleID = brca@meta.data$sample)
head(int_umap)


c("#f4a261","#2a9d8f","#adc178","#f5cac3","#5e60ce","#fad890","#7a529f","#ef476f",
  "#ae9a84","#71c9dd")

celltype_colors = c("#f4a261","#2a9d8f","#adc178","#f5cac3","#7a529f","#E64B35CC","#ef476f","#ae9a84","#71c9dd")

# Calculate the median coordinates of each cluster
int_umap <- int_umap[order(int_umap$cell_type,int_umap$umap_1),]
int_umap$cell_type <- factor(int_umap$cell_type, levels = celltype_levels)

celltype_med <- 
  int_umap %>% group_by(cell_type) %>% dplyr::summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

summary(int_umap)
lb = -14
rb = 17
db = -18
ub = 16.2

w = rb-lb
h = ub-db


ggplot(int_umap, aes(x = umap_1, y=umap_2, color = cell_type)) + 
  geom_point_rast(size = 0.5, alpha = 0.1,
                  shape = 16, stroke = 0, raster.dpi = 1000) +
  scale_color_manual(values = celltype_colors) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"))+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.1,'cm')) +
  xlim(lb,rb) +
  ylim(db,ub) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  geom_segment(data = data.frame(1),
               aes(x = lb + 0.05*w, y = db+0.05*h,   # Draw a horizontal line
                   xend = lb+0.05*w+0.18*w, yend = db+0.05*h),
               colour = "black", linewidth = 1 * 0.3528) + 
  geom_segment(data = data.frame(1),
               aes(x = lb + 0.05*w  , y = db+0.05*h, # Draw a vertical line
                   xend = lb+0.05*w , yend = db+0.05*h+0.18*h),
               colour = "black", size = 1 * 0.3528) + # Add the horizontal tSNE label
  annotate("text", x = lb + 0.05*w + 0.09*w, y = db+0.03*h, label = "UMAP_1",
           color="black",size = 4 * 0.3528, fontface="bold" ) + # Add the vertical tSNE label; fontface sets the font style
  annotate("text", x = lb + 0.02*w, y = db + 0.05*h+0.09*h, label = "UMAP_2",
           color="black",size = 4 * 0.3528, fontface="bold" ,angle=90)  +
  geom_label_repel(aes(label=cell_type),
                   color = "black",
                   label.size = NA,
                   label.padding =0.1,
                   na.rm = TRUE,
                   fill = alpha(c("white"),0.5),
                   size = 4 * 0.3528,
                   fontface="bold",data = celltype_med) + theme(legend.position = "none")

ggsave("E:/Projects/BRCA/results/SC/main_umap_2609.pdf",width=5,height = 5,units = "cm",
       dpi = 300)

# Plot the UMAP distribution of each sample
ggplot(int_umap, aes(x = umap_1, y=umap_2, color = sample_id)) + 
  geom_point_rast(size = 0.5, alpha = 0.1,
                  shape = 16, stroke = 0, raster.dpi = 1000) +
  #scale_color_manual(values = celltype_colors) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"))+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.1,'cm')
        ) +
  xlim(lb,rb) +
  ylim(db,ub) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  geom_segment(data = data.frame(1),
               aes(x = lb + 0.05*w, y = db+0.05*h,   # Draw a horizontal line
                   xend = lb+0.05*w+0.18*w, yend = db+0.05*h),
               colour = "black", linewidth = 1 * 0.3528) + 
  geom_segment(data = data.frame(1),
               aes(x = lb + 0.05*w  , y = db+0.05*h, # Draw a vertical line
                   xend = lb+0.05*w , yend = db+0.05*h+0.18*h),
               colour = "black", size = 1 * 0.3528) + # Add the horizontal tSNE label
  annotate("text", x = lb + 0.05*w + 0.09*w, y = db+0.03*h, label = "UMAP_1",
           color="black",size = 4 * 0.3528, fontface="bold" ) + # Add the vertical tSNE label; fontface sets the font style
  annotate("text", x = lb + 0.02*w, y = db + 0.05*h+0.09*h, label = "UMAP_2",
           color="black",size = 4 * 0.3528, fontface="bold" ,angle=90)  

ggsave("E:/Projects/BRCA/results/SC/main_sample_umap_2609.pdf",width=7,height = 5,units = "cm",
       dpi = 300)

#### dotplot
dotplot_markers <- c(
  "EPCAM","MUC1","KRT18","KRT19","ESR1","PGR","ERBB2", # Epithelial
  "KRT5","KRT17","KRT14","MME", # Basal
  "COL1A1","COL1A2","DCN","PDGFRA","CDH11", # Mysenchymal
  "PECAM1","ENG","CD34","VWF", # Endothelial
  "CD3D","CD3E","CD8A","GZMB","NKG7", # T cell
  "NCAM1","TYROBP","NCR1","KLRD1",
  'MS4A1','CD79A','CD79B',"BLNK", # B cell
  'IGHG1','IGKC','JCHAIN','TNFRSF17', # Plasma
  'CD68','CD80','CD86','CD14','APOE','LYZ',"AIF1" # Myeloid
)

dp = DotPlot(brca, features = dotplot_markers, group.by = "celltype") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

data = dp$data
theme_dotplot = theme(
  axis.text.y = element_text(size = 11),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = 11,angle = 45, hjust = 1),
  axis.ticks.x = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white"),
  axis.line = element_line(colour = "black", size = 0.5),
  legend.key = element_blank()
)
data$id = factor(data$id, levels = rev(celltype_levels))


ggplot(data, aes(x = features.plot, y = id)) + geom_point(aes(size = pct.exp, 
                                                              fill = avg.exp.scaled,
                                                              colour = avg.exp.scaled), shape = 21) +
  
  scale_fill_gradient2(low = "#2166AC", high = "#B2182B", 
                       limits = c(min(data$avg.exp.scaled), max(data$avg.exp.scaled)),
                       midpoint = mean(data$avg.exp.scaled)) +
  scale_color_gradient2(low = "#2166AC", high = "#B2182B", 
                        limits = c(min(data$avg.exp.scaled), max(data$avg.exp.scaled)),
                        midpoint = mean(data$avg.exp.scaled)) +
  scale_size_continuous(range = c(0.2, 2)) +
  guides(colour = FALSE,
         size = guide_legend(title = "Percent Expressed",),
         fill = guide_colorbar(title = "Average Expression")) +
  theme(
    axis.text.y = element_text(colour = "black",size = 4),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(colour = "black",size = 4,angle = 45, hjust = 1),
    axis.ticks.length = unit(1.5, "pt"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5 * 0.3528),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black", linewidth = 0.5 * 0.3528),
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 5),
    #legend.key = element_blank(),
    legend.key.size = unit(6, "pt"),
    legend.background = element_rect(size = 0.2)
  )
ggsave("E:/Projects/BRCA/results/SC/main_markers_dotplot.png",
       width = 12, height = 4.5, units = "cm")
ggsave("E:/Projects/BRCA/results/SC/main_markers_dotplot.pdf",
       width = 12, height = 4.5, units = "cm")
