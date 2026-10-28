##### Cluster B cells #####
b = subset(brca, celltype %in% c("B", "Plasma"))

b = NormalizeData(b)
b = FindVariableFeatures(b)
variable.features = VariableFeatures(b)
if(length(grep("^RP([0-9]+-|[LS])", variable.features)) > 0){
  variable.features = variable.features[-grep("^RP([0-9]+-|[LS])", variable.features)]
}
grep("^RP([0-9]+-|[LS])", variable.features)
"MALAT1" %in% variable.features

b = ScaleData(b, features = variable.features, vars.to.regress = c("nCount_RNA","percent.mt"))
b = RunPCA(b)
ElbowPlot(b, ndims = 50)

b[["RNA"]] = split(b[["RNA"]], f = b$sample)
b = IntegrateLayers(
  object = b, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = FALSE
)
b = JoinLayers(b)

b = FindNeighbors(b, reduction = "harmony", dims = 1:30)
b = FindClusters(b, resolution = 0.8)
b = RunUMAP(b, reduction = "harmony", dims = 1:30)
DimPlot(b, reduction = "umap", label = T)
DimPlot(b, group.by = "sample")


b_markers = FindAllMarkers(b, only.pos = T, logfc.threshold = 1.2, min.pct = 0.2)
top50_b_markers = b_markers %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 50) %>% ungroup()


qc.df = FetchData(b, vars = c("seurat_clusters", "nCount_RNA","nFeature_RNA","percent.mt"))
qc_median = qc.df %>% group_by(seurat_clusters) %>% dplyr::summarise(median_counts = median(nCount_RNA),
                                                                     median_Features = median(nFeature_RNA),
                                                                     median_mt = median(percent.mt))


b_cell_markers = c("TCL1A","YBX3","FCER2","IL4R","CD72","BACH2","IGHD","IGHM", # B naive TCL1A+
                   "NR4A1","CREM", # B naive NR4A2+
                   "ISG15","STAT1","IFI44L","IFI6","IFIT3", # IFN-response  
                   "HSPA1A","DNAJB1","HSPA1B",
                   "MT1X","MT2A","SLC30A1",
                   "TNF","EGR1",
                   "ZNF331","VPS37B","NR4A2","CD69",
                   "CD83","CCR7","NFKBID",
                   "CRIP2", # Bm marker
                   "CRIP1","S100A10","S100A4","ITGB1",
                   "CD27","TNFRSF13B","TXNIP","GPR183",
                   "DUSP4","FCRL5","ZEB2","ITGAX",
                   "ENO1","PSME2","NME1","APEX1","POLD2","POLE3","MYC", # pre-GC
                   "SUGCT","NEIL1","MEF2B",
                   "LMO2","GMDS","MARCKSL1",
                   "STMN1","MKI67","TOP2A",
                   "TXNDC5","MYDGF",
                   "MS4A1",
                   "XBP1","MZB1","JCHAIN",
                   "IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2")
VlnPlot(b, b_cell_markers, stack = T, pt.size = 0)
DotPlot(b, b_cell_markers) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### Annotate cell types #####
Idents(b) = b$seurat_clusters
new.cluster.ids <-c("B_Naive_TCL1A+","Bm_ITGB1+","Bm_GPR183+","Bm_GPR183+","B_MT1X+","ACB_EGR1+",
                    "PC_IGHG+","Bm_ITGB1+","Bm_ITGB1+","PC_IGHG+","B_RALGPS2+",
                    "AtM_DUSP4+","PC_early_MS4A1low","ACB_EGR1+","Bm_ITGB1+","Bm_GPR183+",
                    "Bm_GPR183+","Bm_ITGB1+","Bm_GPR183+","GCB","delete")
length(new.cluster.ids)
head(Idents(b))
names(new.cluster.ids) <- levels(b)
b <- RenameIdents(b, new.cluster.ids)
DimPlot(b, label = T)
b@meta.data$sub_celltype = as.character(Idents(b))

b = subset(b, sub_celltype != "delete")
ncol(b)

b_celltype_levels = c("B_Naive_TCL1A+","B_RALGPS2+","B_MT1X+","Bm_GPR183+","Bm_ITGB1+","ACB_EGR1+",
                      "AtM_DUSP4+","GCB","PC_early_MS4A1low","PC_IGHG+")
all(unique(b$sub_celltype) %in% b_celltype_levels)
all(b_celltype_levels %in% unique(b$sub_celltype))

b$sub_celltype = factor(b$sub_celltype, levels = b_celltype_levels)
levels(b$sub_celltype)
DimPlot(b, reduction = "umap", label = T, group.by = "sub_celltype")

b_celltype_label = paste0("c",1:length(b_celltype_levels))
names(b_celltype_label) = b_celltype_levels
b$sub_label = as.character(b_celltype_label[as.character(b$sub_celltype)])
b$sub_label = factor(b$sub_label, levels = paste0("c",1:length(b_celltype_levels)))


##### Plot the UMAP and DotPlot #####
library('ggsci')
library('ggplot2')
library("scales")
library(Scillus)
library(ggrepel)
library(tidyverse)
library(scales)
library(ggrastr)

int_umap = b@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_label = b@meta.data$sub_label) %>%
  cbind(cell_type = b@meta.data$sub_celltype) %>%
  cbind(SampleID = b@meta.data$sample)
head(int_umap)


celltype_colors <- c("#f4a261","#2a9d8f","#5e60ce","#f5cac3","#adc178","#fad890","#7a529f","#ef476f",
                     "#ae9a84","#71c9dd")

int_umap <- int_umap[order(int_umap$cell_type,int_umap$umap_1),]
int_umap$cell_label <- factor(int_umap$cell_label, levels = paste0("c",1:length(unique(int_umap$cell_label))))
levels(int_umap$cell_label)

celltype_med <- 
  int_umap %>% group_by(cell_label) %>% dplyr::summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

summary(int_umap)
lb = -14
rb = 7.5
db = -8
ub = 11

w = rb-lb
h = ub-db


ggplot(int_umap, aes(x = umap_1, y=umap_2, color = cell_label)) + geom_point_rast(size = 0.05, alpha = 0.7,
                                                                                  raster.dpi = 600) +
  scale_color_manual(values = celltype_colors,
                     breaks = paste0("c",1:length(unique(int_umap$cell_label))),
                     labels = paste0("c",1:length(unique(int_umap$cell_label)),": ",b_celltype_levels)) +
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
        legend.key.size = unit(0.5,'cm')) +
  xlim(lb,rb) +
  ylim(db,ub) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_segment(aes(x = lb + 0.05*w, y = db+0.05*h,   # Draw a horizontal line
                   xend = lb+0.05*w+0.18*w, yend = db+0.05*h),
               colour = "black", linewidth = 1 * 0.3528) + 
  geom_segment(aes(x = lb + 0.05*w  , y = db+0.05*h, # Draw a vertical line
                   xend = lb+0.05*w , yend = db+0.05*h+0.18*h),
               colour = "black", linewidth = 1 * 0.3528) + # Add the horizontal tSNE label
  annotate("text", x = lb + 0.05*w + 0.09*w, y = db+0.03*h, label = "UMAP_1",
           color="black",size = 4 * 0.3528, fontface="bold" ) + # Add the vertical tSNE label; fontface sets the font style
  annotate("text", x = lb + 0.02*w, y = db + 0.05*h+0.09*h, label = "UMAP_2",
           color="black",size = 4 * 0.3528, fontface="bold" ,angle=90)  +
  geom_label(aes(label=cell_label),
                   color = "black",
                   label.size = NA,
                   size = 4 * 0.3528,
                   label.padding = unit(0.1, "lines"),
                   na.rm = TRUE,
                   fill = alpha(c("white"),0.5),
                   fontface="bold",data = celltype_med) #+ theme(legend.position = "none")

ggsave("E:/Projects/BRCA/results/SC/B_Plasma/B_Plasma_umap.png",width=10,height = 6,units = "cm",
       dpi = 600)



### dotplot
dotplot_markers <- c(
  "MS4A1",
  "TCL1A","YBX3","FCER2","IL4R", # Naive
  "RALGPS2","PTPRJ",
  "MT1X","MT2A", # B_MT1X+
  "TXNIP", # Bm_TXNIP+
  "CCR7","GPR183", # Bm_GPR183+
  "ITGB1","CRIP1","S100A10","S100A4", # Bm_ITGB1+
  "EGR1","TNF","CD83","CD69", # ACB_EGR1+
  "DUSP4","FCRL5","ZEB2","ITGAX", # AtM
  "SUGCT","NEIL1","LMO2","GMDS", # GCB
  "XBP1","MZB1","JCHAIN" # Plasma
)


dp = DotPlot(b, features = dotplot_markers, group.by = "sub_celltype") + 
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
  legend.key = element_blank()  # Remove the legend key background
)
data$id = factor(data$id, levels = rev(b_celltype_levels))
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
  ) +
  guides(colour = FALSE,
         size = guide_legend(title = "Percent Expressed",),
         fill = guide_colorbar(title = "Average Expression"))

ggsave("E:/Projects/BRCA/results/SC/B_Plasma/B_Plasma_markers_dotplot.pdf",
       width = 12, height = 5, units = "cm")

