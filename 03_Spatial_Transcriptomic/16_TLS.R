tls_sample_list = c("P36","P38","P32","P27","P39","P47","P46","P33","P49","P34","P37","P43_2","P45","P5","P8","P26",  
                    "P4","P43_1","P42","P40_2","P11","P30","P10")

library(dbscan)
library(tidyverse)

spID = tls_sample_list[1]

##### Perform DBSCAN clustering on lymphocyte aggregate regions #####
for(spID in tls_sample_list[2:length(tls_sample_list)]){
  print(spID)
  
  sp_coords = read.csv(paste0("E:/Projects/BRCA/results/ST/Niche2607/", spID, "/coord_", spID, ".csv"))
  sp_niche_info = filter(niche_info, sample == spID)
  sp_niche_info = filter(sp_niche_info, niche %in% c("niche_6","niche_7"))
  sp_coords = filter(sp_coords, bin %in% sp_niche_info$bin)
  
  
  wanted_r = read.table(paste0("E:/Projects/BRCA/results/ST/Niche2607/", spID, "/wanted_r.txt"), sep = "\t")
  wanted_r = wanted_r$V1
  wanted_r = wanted_r / 80 * 30
  
  db_res <- dbscan(sp_coords[, c("x", "y")], eps = wanted_r, minPts = 4)
  sp_coords$dbscan_cluster = db_res$cluster

  cluster_count = data.frame(table(sp_coords$dbscan_cluster))
  cluster_count = filter(cluster_count, Freq > 500)
  cluster_count = cluster_count[order(cluster_count$Freq, decreasing = TRUE),]
  ordered_name = 1:nrow(cluster_count)
  names(ordered_name) = as.character(cluster_count$Var1)
  
  sp_coords = filter(sp_coords, dbscan_cluster %in% cluster_count$Var1)
  
  sp_coords$dbscan_cluster = ordered_name[as.character(sp_coords$dbscan_cluster)]
  sp_coords$dbscan_cluster = paste0("dbscan_cluster_",sp_coords$dbscan_cluster)
  table(sp_coords$dbscan_cluster)
  sp_coords$sample = spID
  sp_coords$celltype = sp_niche_info$celltype[match(sp_coords$bin, sp_niche_info$bin)]
  
  object = readRDS(file = paste0("E:/Projects/BRCA/Intermediate_Data/ST/object_unified/",spID,"_seg.rds"))
  object$tls_dbscan_cluster = "undefined"
  object$tls_dbscan_cluster[sp_coords$bin] = sp_coords$dbscan_cluster
  object$tls_dbscan_cluster = factor(object$tls_dbscan_cluster, levels = c(paste0("dbscan_cluster_", 1:nrow(cluster_count)),"undefined"))
  
  output.dir = paste0("E:/Projects/BRCA/results/ST/Niche2607/TLS/", spID)
  if(!(dir.exists(output.dir))){
    dir.create(output.dir)
  }
  
  alternative_color_celltype = c("#f4a261","#2a9d8f","#a8dadc","#f5cac3","#ef476f",
                                 "#4d9d71","#9d4edd","#adc178","#eb6f5d","#5e60ce",
                                 "#ae9a84","#DC0000CC","#8491B4CC","#f5af99","#71c9dd",
                                 "#A20056B2","#fad890","#5a189a","#3c6e71","#eaf4f4",
                                 "#ffd60a","#b7e4c7","#ff7b00","#001233","#9f86c0",
                                 "#ccdbfd","#41ead4","#495867","#004e98","#603808",
                                 "#ff5d8f","#70e000")
  
  
  sp_cluster_color = alternative_color_celltype[1:nrow(cluster_count)]
  sp_cluster_color = c(sp_cluster_color, "#d3d3d3")
  plt = SpatialDimPlot(object, group.by = "tls_dbscan_cluster", pt.size.factor = (seurat_spot_size[spID] + 1), image.alpha = 0) +
    scale_fill_manual(values = sp_cluster_color) + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  ggsave(filename = paste0(output.dir, "/TLS_dbscan_cluster_spatial.png"),
         plot = plt, width = 23, height = 17, units = "cm", dpi = 1000)
  
  write.csv(sp_coords, paste0(output.dir, "/TLS_dbscan_cluster_info.csv"))
  
  # Calculate the 98th percentile of TLS gene expression for each cluster
  tls_markers = c("CD3D", "MS4A1", "CR2", "FCER2", "BCL6", "MKI67")
  expr_df <- FetchData(object, vars = tls_markers, layer = "data")
  min_max_normalize <- function(x) {
    # Prevent division by zero when a gene has identical expression in all cells
    if (max(x) - min(x) == 0) {
      return(rep(0, length(x))) 
    }
    return((x - min(x)) / (max(x) - min(x) + 1e-12))
  }
  normalized_df <- as.data.frame(lapply(expr_df, min_max_normalize))
  rownames(normalized_df) = rownames(expr_df)
  normalized_df$cell_id = rownames(normalized_df)
  
  normalized_df = normalized_df[sp_coords$bin,]
  all(rownames(normalized_df) == sp_coords$bin)
  normalized_df$dbscan_cluster = sp_coords$dbscan_cluster
  
  cluster_gene_98th <- normalized_df %>%
    group_by(dbscan_cluster) %>%
    summarise(across(
      where(is.numeric), 
      ~ quantile(.x, probs = 0.98, na.rm = TRUE)
    ))
  cluster_gene_98th = cluster_gene_98th %>% add_column(.before = "dbscan_cluster", sample = spID)
  write.csv(cluster_gene_98th, file = paste0(output.dir, "/cluster_TLS_gene_98th_percentile.csv"))
}


##### Load expression data from all clusters and perform clustering #####
all_tls_expr = data.frame()
for(spID in tls_sample_list){
  sp_tls_expr = read.csv(paste0("E:/Projects/BRCA/results/ST/Niche2607/TLS/", spID, "/cluster_TLS_gene_98th_percentile.csv"),
                         row.names = 1)
  sp_tls_expr = sp_tls_expr %>% add_column(cluster_full = paste0(sp_tls_expr$sample, "_", sp_tls_expr$dbscan_cluster),
                                           .after = "dbscan_cluster")
  all_tls_expr = rbind(all_tls_expr, sp_tls_expr)
}

# Cluster using CR2 and FCER2
set.seed(50)
maturation_data = all_tls_expr[,c("CR2","FCER2")]
kmeans_result <- kmeans(maturation_data, centers = 2, nstart = 25)
all_tls_expr$maturation_cluster = kmeans_result$cluster

all_tls_expr$maturation_cluster[all_tls_expr$maturation_cluster == 1] = "fTLS"
all_tls_expr$maturation_cluster[all_tls_expr$maturation_cluster == 2] = "others"

write.csv(all_tls_expr, file = "E:/Projects/BRCA/results/ST/Niche2607/TLS/TLS_cluster_classification.csv")
all_tls_expr = read.csv("E:/Projects/BRCA/results/ST/Niche2607/TLS/TLS_cluster_classification.csv", row.names = 1)

##### B-cell state scores in lymphocyte aggregates and follicle-like TLS ######
for(spID in tls_sample_list){
  print(spID)
  object = readRDS(file = paste0("E:/Projects/BRCA/Intermediate_Data/ST/object_unified/",spID,"_seg.rds"))
  sp_tls_b_cell = filter(all_tls_b_cell, sample == spID)
  sp_b = subset(object, cells = sp_tls_b_cell$bin)
  
  sp_b$lymphocyte_region = "undefined"
  sp_b$lymphocyte_region[sp_tls_b_cell$bin] = sp_tls_b_cell$lymphocyte_region
  if(!exists("b_region")){
    b_region = sp_b
  }else{
    b_region = merge(b_region, sp_b)
  }
}

b_region = JoinLayers(b_region)

# GCB score
b = readRDS("E:/Projects/BRCA/Intermediate_Data/SC/SC_b_object.rds")

all_b_markers = FindAllMarkers(b)
gcb_markers = filter(all_b_markers, cluster == "GCB")
gcb_markers = gcb_markers$gene[1:100]

b_region = AddModuleScore(b_region, features = list(gcb_markers), name = "gcb_score")
gcb_score_plt = FetchData(b_region, c("lymphocyte_region", "gcb_score1"))
gcb_score_plt$lymphocyte_region = factor(gcb_score_plt$lymphocyte_region, levels = c("lymphocyte_region","fTLS"))
gcb_score_plt$cluster_full = all_tls_b_cell$cluster_full[match(rownames(gcb_score_plt), all_tls_b_cell$bin)]



cluster_score_98th <- gcb_score_plt %>%
  group_by(cluster_full) %>%
  summarise(across(
    where(is.numeric), 
    ~ quantile(.x, probs = 0.98, na.rm = TRUE)
  ))
cluster_score_98th$lymphocyte_region = gcb_score_plt$lymphocyte_region[match(cluster_score_98th$cluster_full,
                                                                             gcb_score_plt$cluster_full)]

library(ggsignif)

ggplot(cluster_score_98th, mapping = aes(x = lymphocyte_region, y = gcb_score1)) +
  geom_boxplot(mapping = aes(fill = lymphocyte_region), width = 0.6, alpha = 0.7, outlier.shape = NA,
               box.linewidth = 0.5 * 0.3528, median.linewidth = 0.8 * 0.3528,
               whisker.linewidth = 0.5 * 0.3528) +
  geom_jitter(mapping = aes(color = lymphocyte_region),
              show.legend = FALSE, alpha = 0.4, size = 0.2, width = 0.3) +
  scale_fill_manual(values = c("#fad890","#2a9789")) +
  scale_color_manual(values = c("#fad890","#2a9789")) +
  scale_x_discrete(labels = c("lymphocyte aggregates", "folicular-like TLS")) +
  ylab("98th percentile of GCB signature score") +
  geom_signif(comparisons = list(c("lymphocyte_region","fTLS")),
              #map_signif_level = function(p) sprintf("italic(P) == %.3f", p),
              tip_length = 0.02  * 0.3528,
              size =  0.1, textsize = 1.5,
              parse = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        #axis.line.x = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 5),
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
ggsave(filename = "E:/Projects/BRCA/results/ST/B/GCB_score_98th_percentile.pdf", width = 5, height = 5,
       units = "cm", dpi = 300)


# cell cycle score
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

hallmarker_sets = msigdbr(species = "Homo sapiens", category = "H")
hallmarker_sets = subset(hallmarker_sets, select = c("gs_name","gene_symbol")) %>% as.data.frame()

cell_cycle_genes = filter(hallmarker_sets, gs_name == "HALLMARK_G2M_CHECKPOINT")
cell_cycle_genes = cell_cycle_genes$gene_symbol

b_region = AddModuleScore(b_region, features = list(cell_cycle_genes), name = "cell_cycle_score")
cell_cycle_score_plt = FetchData(b_region, c("lymphocyte_region", "cell_cycle_score1"))
cell_cycle_score_plt$lymphocyte_region = factor(cell_cycle_score_plt$lymphocyte_region, levels = c("lymphocyte_region","fTLS"))

library(ggsignif)
ggplot(cell_cycle_score_plt, mapping = aes(x = lymphocyte_region, y = cell_cycle_score1)) +
  geom_violin(mapping = aes(fill = lymphocyte_region), alpha = 0.5, linewidth = 0.5 * 0.3528) +
  geom_boxplot(mapping = aes(fill = lymphocyte_region), width = 0.4, alpha = 0.5, outlier.shape = NA,
               box.linewidth = 0.5 * 0.3528, median.linewidth = 0.8 * 0.3528,
               whisker.linewidth = 0.5 * 0.3528) +
  scale_fill_manual(values = c("#fad890","#2a9789")) +
  scale_x_discrete(labels = c("lymphocyte aggregates", "folicular-like TLS")) +
  ylab("Cell cycle score in B cells") +
  geom_signif(comparisons = list(c("lymphocyte_region","fTLS")),
              #map_signif_level = function(p) sprintf("italic(P) == %.3f", p),
              tip_length = 0.02  * 0.3528,
              size =  0.1, textsize = 1.5,
              parse = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        #axis.line.x = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 5),
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
ggsave(filename = "E:/Projects/BRCA/results/ST/B/cell_cycle_score_fTLS_vs_lymphocyte_region.pdf", width = 5, height = 5,
       units = "cm", dpi = 300)
