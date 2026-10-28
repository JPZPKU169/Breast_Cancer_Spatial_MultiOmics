##### Cluster myeloid cells #####
mye = subset(brca, celltype == "Myeloid")

mye = readRDS(file = "E:/Projects/BRCA/Intermediate_Data/SC/mye_filtered_Seurat_Object.rds")

immunogloblin_and_tcr = read.table(file = "D:/references/BioMart/Immunogloblin_and_TCR_genes.txt", sep = "\t", header = T)

mye = FindVariableFeatures(mye)
variable.features = VariableFeatures(mye)
table(variable.features %in% immunogloblin_and_tcr$Gene.name)
variable.features = variable.features[!(variable.features %in% immunogloblin_and_tcr$Gene.name)]
grep("^RP([0-9]+-|[LS])", variable.features)
variable.features = variable.features[-grep("^RP([0-9]+-|[LS])", variable.features)]
"MALAT1" %in% variable.features

mye = ScaleData(mye, features = variable.features, vars.to.regress = c("nCount_RNA","percent.mt"))
mye = RunPCA(mye, features = variable.features, npcs = 50)
ElbowPlot(mye, ndims = 50)

mye = FindNeighbors(mye, reduction = "pca", dims = 1:40)
mye = FindClusters(mye, resolution = 1)
mye = RunUMAP(mye, reduction = "pca", dims = 1:40)

### harmony
mye[["RNA"]] = split(mye[["RNA"]], f = mye$sample)
mye = IntegrateLayers(
  object = mye, method = HarmonyIntegration,
  new.reduction = "harmony",theta = 2.5,
  verbose = FALSE
)
mye = JoinLayers(mye)

mye = FindNeighbors(mye, reduction = "harmony", dims = 1:40)
mye = FindClusters(mye, resolution = 2)
mye = RunUMAP(mye, reduction = "harmony", dims = 1:40)
DimPlot(mye, reduction = "umap", label = T)
DimPlot(mye, group.by = "sample")

mye_markers = FindAllMarkers(mye, min.pct = 0.1, only.pos = T)
mye_markers = filter(mye_markers, p_val_adj < 0.01)
top50_deg_mye = mye_markers %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 50) %>% ungroup()


qc.df = FetchData(mye, vars = c("seurat_clusters", "nCount_RNA","nFeature_RNA","percent.mt"))
qc_median = qc.df %>% group_by(seurat_clusters) %>% dplyr::summarise(median_counts = median(nCount_RNA),
                                                                  median_Features = median(nFeature_RNA),
                                                                  median_mt = median(percent.mt))

##### Annotate cell types #####
Idents(mye) = mye$seurat_clusters
new.cluster.ids <-c("Macro_C6_SPP1+","cDC2_CD1C+","Macro_C6_SPP1+","Mast","Macro_C1_CX3CR1+","Mono_FCN1+",
                    "Macro_C4_FOLR2+","pDC","Macro_C5_PDE4C+","Mono_CXCL10+","Macro_C6_SPP1+",
                    "Macro_C2_CADM1+","Macro_C1_CX3CR1+","cDC2_CD1C+","Mono_FCN1+","Macro_C6_SPP1+",
                    "Macro_C4_FOLR2+","unknown","Macro_C4_FOLR2+","Mono_FCN1+","Macro_C7_MKI67+",
                    "cDC3_LAMP3+","Macro_C1_CX3CR1+","Macro_C4_FOLR2+","Macro_C3_CXCL10+","cDC1_CLEC9A+",
                    "delete","Macro_C2_CADM1+","Macro_C6_SPP1+","Macro_C2_CADM1+","cDC2_CD1C+",
                    "DC_CD207+","Macro_C6_SPP1+","Neutrophils")

length(new.cluster.ids)
head(Idents(mye))
names(new.cluster.ids) <- levels(mye)
mye <- RenameIdents(mye, new.cluster.ids)
DimPlot(mye, label = T)
mye@meta.data$sub_celltype = as.character(Idents(mye))
ncol(mye)
mye = subset(mye, sub_celltype != "delete")

mye_celltype_levels = c("cDC1_CLEC9A+","cDC2_CD1C+","cDC3_LAMP3+","DC_CD207+","pDC",
                        "Mono_FCN1+","Mono_CXCL10+",
                        "Macro_C1_CX3CR1+","Macro_C2_CADM1+","Macro_C3_CXCL10+","Macro_C4_FOLR2+",
                        "Macro_C5_PDE4C+","Macro_C6_SPP1+","Macro_C7_MKI67+",
                        "Mast","Neutrophils", "unknown")
all(unique(mye$sub_celltype) %in% mye_celltype_levels)
all(mye_celltype_levels %in% unique(mye$sub_celltype))

mye$sub_celltype = factor(mye$sub_celltype, levels = mye_celltype_levels)
levels(mye$sub_celltype)
DimPlot(mye, reduction = "umap", label = T, group.by = "sub_celltype")

mye_celltype_label = paste0("c",1:length(mye_celltype_levels))
names(mye_celltype_label) = mye_celltype_levels
mye$sub_label = as.character(mye_celltype_label[as.character(mye$sub_celltype)])
mye$sub_label = factor(mye$sub_label, levels = paste0("c",1:17))


##### Plot the UMAP and DotPlot #####
library('ggsci')
library('ggplot2')
library("scales")
library(Scillus)
library(ggrepel)
library(tidyverse)
library(scales)
library(ggrastr)

int_umap = mye@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_label = mye@meta.data$sub_label) %>%
  cbind(cell_type = mye@meta.data$sub_celltype) %>%
  cbind(SampleID = mye@meta.data$sample)
head(int_umap)

celltype_med <- 
  int_umap %>% group_by(cell_label) %>% dplyr::summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

summary(int_umap)
lb = -11
rb = 15.5
db = -7
ub = 13.5

w = rb-lb
h = ub-db

ggplot(int_umap, aes(x = umap_1, y=umap_2, color = cell_label)) + 
  geom_point_rast(size = 0.001, alpha = 0.5, raster.dpi = 600) +
  scale_color_manual(values = celltype_colors_ordered) +
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

ggsave("E:/Projects/BRCA/results/SC/Myeloid/Myeloid_umap_final_2.pdf", width=8,height = 6 ,units = "cm",
       dpi = 300)


ggplot(int_umap, aes(x = umap_1, y=umap_2, color = cell_label)) + 
  geom_point_rast(size = 0.001, alpha = 0.5, raster.dpi = 1000) +
  scale_color_manual(values = celltype_colors_ordered) +
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
           color="black",size = 4 * 0.3528, fontface="bold" ,angle=90)

ggsave("E:/Projects/BRCA/results/SC/Myeloid/Myeloid_umap_final_2_no_labelled.pdf", width=8,height = 6 ,units = "cm",
       dpi = 300)

# Plot the UMAP distribution of each sample
ggplot(int_umap, aes(x = umap_1, y=umap_2, color = SampleID)) + geom_point(size = 1, alpha = 0.7) +
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
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.6,'cm')
  ) +
  xlim(lb,rb) +
  ylim(db,ub) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_segment(aes(x = lb + 0.05*w, y = db+0.05*h,   # Draw a horizontal line
                   xend = lb+0.05*w+0.18*w, yend = db+0.05*h),
               colour = "black", linewidth = 1) + 
  geom_segment(aes(x = lb + 0.05*w  , y = db+0.05*h, # Draw a vertical line
                   xend = lb+0.05*w , yend = db+0.05*h+0.18*h),
               colour = "black", size = 1) + # Add the horizontal tSNE label
  annotate("text", x = lb + 0.05*w + 0.09*w, y = db+0.03*h, label = "UMAP_1",
           color="black",size = 4, fontface="bold" ) + # Add the vertical tSNE label; fontface sets the font style
  annotate("text", x = lb + 0.02*w, y = db + 0.05*h+0.09*h, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)  #+

ggsave("E:/Projects/BRCA/results/SC/Myeloid/Myeloid_sample_umap.pdf",width=24,height = 20,units = "cm",
       dpi = 300)

### dotplot
dotplot_markers = c("CLEC9A","XCR1","PTDSS1", # DC_CLEC9A
                    "CD1C","CLEC10A","FCER1A", # DC_CD1C
                    "LAMP3","FSCN1","CCL22", # DC_LAMP3
                    "CD207","CD1A","S100B", # DC_CD207
                    "PLD4","LILRA4","SPIB", # pDC
                    "FCN1","CD300E","CLEC12A","S100A8","LILRB2", # Monocytes
                    "CX3CR1","C3", # Macro_C1_CX3CR1+
                    "CADM1","TREM2", # Macro_C2_CADM1+
                    "CXCL10","CXCL9","CXCL11","ISG15", # Macro_C3_CXCL10+
                    "FOLR2","SLC40A1","MRC1","PLTP","LYVE1", # Macro_C4_FOLR2+
                    "PDE4C", "SNX31", # Macro_C5_PDE4C+
                    "SPP1","MMP9","RETN","CCL20", # Macro_6_SPP1
                    "MKI67","STMN1","TOP2A", # Macro_C9_MKI67
                    "TPSB2","CPA3","MS4A2", # Mast
                    "CXCL8","G0S2","CSF3R","FCGR3B" # Neutrophils
)

# Continue plotting
dp = DotPlot(mye, features = dotplot_markers, group.by = "sub_celltype") + 
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
data$id = factor(data$id, levels = rev(mye_celltype_levels))
data = filter(data, id != "unknown")
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

ggsave("E:/Projects/BRCA/results/SC/Myeloid/Myeloid_markers_dotplot_final_2.pdf",
       width = 14, height = 5.5, units = "cm")


##### Scoring #####
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

macrophage_color = c("#7a529f","#ef476f","#fe9f69","#f4d47f","#b7e4c7","#70d6ff")


macro = subset(mye, sub_celltype %in% c("Macro_C1_CX3CR1+","Macro_C2_CADM1+","Macro_C3_CXCL10+","Macro_C4_FOLR2+",
                                        "Macro_C5_PDE4C+","Macro_C6_SPP1+"))

gobp_sets = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gobp_sets = subset(gobp_sets, select = c("gs_name","gene_symbol")) %>% as.data.frame()

kegg_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_sets = subset(kegg_sets, select = c("gs_name","gene_symbol")) %>% as.data.frame()

cholesteral_efflux = filter(gobp_sets, gs_name == "GOBP_CHOLESTEROL_EFFLUX")
cholesteral_efflux = list(unique(cholesteral_efflux$gene_symbol))
steroid_biosynthesis = filter(kegg_sets, gs_name == "KEGG_STEROID_BIOSYNTHESIS")
steroid_biosynthesis = list(unique(steroid_biosynthesis$gene_symbol))

macro = AddModuleScore(macro, cholesteral_efflux, name = "cholesteral_efflux")
macro = AddModuleScore(macro, steroid_biosynthesis, name = "steroid_biosynthesis")

plot.df = FetchData(macro, c("sub_celltype","cholesteral_efflux1","steroid_biosynthesis1"))


plot.df$sub_celltype = factor(plot.df$sub_celltype, 
                              levels = c("Macro_C1_CX3CR1+","Macro_C2_CADM1+","Macro_C3_CXCL10+","Macro_C4_FOLR2+",
                                         "Macro_C5_PDE4C+","Macro_C6_SPP1+"))

compare_list = list(c("Macro_C1_CX3CR1+","Macro_C6_SPP1+"),c("Macro_C2_CADM1+","Macro_C6_SPP1+"),
                    c("Macro_C3_CXCL10+","Macro_C6_SPP1+"),
                    c("Macro_C4_FOLR2+","Macro_C6_SPP1+"),
                    c("Macro_C5_PDE4C+","Macro_C6_SPP1+"))


# cholesteral_efflux1
ggplot() +
  geom_violin(plot.df, mapping = aes(x = sub_celltype, y = cholesteral_efflux1, fill = sub_celltype), alpha = 0.7,
              linewidth = 0.5 * 0.3528) +
  geom_boxplot(plot.df, mapping = aes(x = sub_celltype, y = cholesteral_efflux1, fill = sub_celltype), width = 0.35,
               show.legend = FALSE, outlier.shape = NA, box.linewidth = 0.5 * 0.3528, median.linewidth = 0.8 * 0.3528,
               whisker.linewidth = 0.5 * 0.3528) +
  scale_fill_manual(values = macrophage_color) +
  geom_signif(plot.df, mapping = aes(x = sub_celltype, y = cholesteral_efflux1),
              comparisons = compare_list,
              step_increase = 0.07,
              map_signif_level = function(p_vec) {
                sapply(p_vec, function(p) {
                  if (p < 1e-3) {
                    sci_str <- formatC(p, format = "e", digits = 2)
                    sprintf("italic(P) == '%s'", sci_str)
                  } else {
                    sprintf("italic(P) == '%.3f'", p)
                  }
                })
              },
              tip_length = 0.02  * 0.3528,
              size =  0.1, textsize = 1.5,
              parse = TRUE) +
  ylab("GOBP_CHOLESTEROL_EFFLUX") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        #axis.line.x = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 4),
        axis.ticks.length = unit(1.5, "pt"),
        axis.ticks = element_line(colour = "black", linewidth = 0.5 * 0.3528),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", linewidth = 0.5 * 0.3528),
        legend.text = element_text(size = 4),
        legend.title = element_blank(),
        #legend.key = element_blank(),
        legend.key.size = unit(6, "pt"),
        legend.background = element_rect(size = 0.2))

ggsave(filename = "E:/Projects/BRCA/results/SC/Myeloid/cholesteral_efflux_score_final.pdf",
       width = 8, height = 4.5, units = "cm", dpi = 300)

# steroid_biosynthesis1
ggplot() +
  geom_violin(plot.df, mapping = aes(x = sub_celltype, y = steroid_biosynthesis1, fill = sub_celltype), alpha = 0.7,
              linewidth = 0.5 * 0.3528) +
  geom_boxplot(plot.df, mapping = aes(x = sub_celltype, y = steroid_biosynthesis1, fill = sub_celltype), width = 0.35,
               show.legend = FALSE, outlier.shape = NA, box.linewidth = 0.5 * 0.3528, median.linewidth = 0.8 * 0.3528,
               whisker.linewidth = 0.5 * 0.3528) +
  #scale_x_discrete(labels = custom_labels) +
  geom_signif(plot.df, mapping = aes(x = sub_celltype, y = steroid_biosynthesis1),
              comparisons = compare_list,
              step_increase = 0.07,
              tip_length = 0.02  * 0.3528,
              size =  0.1, textsize = 1.5,
              parse = TRUE) +
  scale_fill_manual(values = macrophage_color) +
  ylab("KEGG_STEROID_BIOSYNTHESIS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        #axis.line.x = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 4),
        axis.ticks.length = unit(1.5, "pt"),
        axis.ticks = element_line(colour = "black", linewidth = 0.5 * 0.3528),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", linewidth = 0.5 * 0.3528),
        legend.text = element_text(size = 4),
        legend.title = element_blank(),
        #legend.key = element_blank(),
        legend.key.size = unit(6, "pt"),
        legend.background = element_rect(size = 0.2))

ggsave(filename = "E:/Projects/BRCA/results/SC/Myeloid/KEGG_steroid_biosynthesis_score_final.pdf",
       width = 8, height = 4.5, units = "cm", dpi = 300)


##### DEG #####
macro = subset(mye, sub_celltype %in% c("Macro_C1_CX3CR1+","Macro_C2_CADM1+","Macro_C3_CXCL10+","Macro_C4_FOLR2+",
                                        "Macro_C5_PDE4C+","Macro_C6_SPP1+"))

macro$deg_label = "others"
macro$deg_label[macro$sub_celltype %in% c("Macro_C6_SPP1+")] = "LAM"
Idents(macro) = macro$deg_label
table(macro$deg_label)

lam_deg = FindMarkers(macro, ident.1 = "LAM", ident.2 = "others", logfc.threshold = 0)
lam_deg$gene = rownames(lam_deg)
write.csv(lam_deg, file = "E:/Projects/BRCA/results/SC/Myeloid/DEG_LAM_vs_other_macrophages_all.csv")
lam_deg = read.csv(file = "E:/Projects/BRCA/results/SC/Myeloid/DEG_LAM_vs_other_macrophages_all.csv", row.names = 1)

lam_deg_sig = filter(lam_deg, p_val_adj < 0.01)


lam_deg$cholesteral_associated = "no"
lam_deg$cholesteral_associated[rownames(lam_deg) %in% cholesteral_label$gene] = "yes"
lam_deg$cholesteral_associated = factor(lam_deg$cholesteral_associated, levels = c("yes","no"))
lam_deg$log_p_adj = -log10(lam_deg$p_val_adj)
lam_deg$log_p_adj[lam_deg$log_p_adj > 300] = 300

cholesteral_label = filter(lam_deg_cholesterol, avg_log2FC > 0.25 & p_val_adj < 0.01)
cholesteral_label$gene = rownames(cholesteral_label)

lam_deg_others = filter(lam_deg, cholesteral_associated == "no")
lam_deg_choles = filter(lam_deg, cholesteral_associated == "yes")

library(ggrepel)
ggplot() + 
  geom_point(lam_deg_others, mapping = aes(x = avg_log2FC, y = log_p_adj, color = cholesteral_associated), 
                      show.legend = FALSE, alpha = 0.5, size = 0.2) +
  geom_point(lam_deg_choles, mapping = aes(x = avg_log2FC, y = log_p_adj, color = cholesteral_associated), 
             show.legend = FALSE, alpha = 0.5, size = 0.2) +
  scale_color_manual(values = c("red","grey70")) +
  geom_vline(aes(xintercept=c(-0.25,0.25)),color="black",linetype="dashed", linewidth = 0.5 * 0.3528) +
  geom_text_repel(cholesteral_label, mapping = aes(x = avg_log2FC, y=log_p_adj,label=gene, fontface = "italic"),
                  box.padding = 0.1,segment.color = "black", max.overlaps = 200, force = 3, force_pull = 3,
                  segment.size = 0.1, size = 3 * 0.3528) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted P") +
  theme(
    axis.text.x = element_text(color = "black", hjust = 0.5, size = 4),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 5),
    axis.title.y = element_text(color = "black", size = 5),
    axis.ticks.length = unit(1.5, "pt"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5 * 0.3528),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black", linewidth = 0.5 * 0.3528),
    #plot.background = element_rect(color = "black",fill = "white"),
    legend.text = element_text(size = 4),
    legend.title = element_blank(),
    #legend.key = element_blank(),
    legend.key.size = unit(6, "pt"),
    legend.background = element_rect(size = 0.2)
    )

ggsave("E:/Projects/BRCA/results/SC/Myeloid/SPP1+_DEG_cholesteral_labelled_final.pdf", 
       width = 6, height = 4.7, dpi = 300, units = "cm")


##### Calculate each macrophage signature #####
macro_markers = FindAllMarkers(macro)

macro_markers_top100 = macro_markers %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 100) %>% ungroup()
write.csv(macro_markers_top100, 
          file = "E:/Projects/BRCA/results/SC/Myeloid/macro_signatures.csv")
macro_markers_top100 = read.csv("E:/Projects/BRCA/results/SC/Myeloid/macro_signatures.csv", row.names = 1)
