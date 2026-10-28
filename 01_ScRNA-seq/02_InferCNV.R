library("infercnv")

##### Obtain the gene order file #####
library(AnnoProbe)
geneInfor <- annoGene(rownames(epi),"SYMBOL","human")
colnames(geneInfor)
geneInfor <- geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]
rownames(geneInfor) <- geneInfor$SYMBOL
geneInfor$SYMBOL <- NULL
saveRDS(geneInfor, file = "E:/Projects/BRCA/Intermediate_Data/SC/geneInfor_for_SC_infercnv")


all_ref_object = subset(brca, celltype %in% c("T","B","Myeloid"))
epi = subset(brca, celltype == "Epithelial")

sample_list_sc = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13","P14","P15","P17","P18","P19","P20",
                   "P21","P22","P23","P24","P25","P26","P27","P28","P29","P30","P31","P32","P37","P39","P40","P41","P42",
                   "P48","P49","P50","P51"
)
all(sample_list_sc %in% brca$sample)
all(brca$sample %in% sample_list_sc)

for(spID in sample_list_sc){
  
  ref_object = subset(all_ref_object, sample == spID)
  ref_group_to_remove = names(which(table(ref_object$celltype) < 10))
  ref_object = ref_object[,!(ref_object$celltype %in% ref_group_to_remove)]
  sp_tumor = subset(epi, sample == spID)
  anno_ref = data.frame(V1=colnames(ref_object), V2=paste(ref_object$celltype))
  anno_tumor = data.frame(V1=colnames(sp_tumor), V2="Epithelial Cells")
  anno_cnv = rbind(anno_ref, anno_tumor)
  rownames(anno_cnv) = anno_cnv$V1
  anno_cnv$V1 = NULL
  matrix = cbind(as.matrix(LayerData(ref_object, assay = "RNA", layer = "counts")), 
                 as.matrix(LayerData(sp_tumor, assay = "RNA", layer = "counts")))
  all_ref_celltype = c('T',"B","Myeloid")
  ref_celltype = all_ref_celltype[all_ref_celltype %in% anno_cnv$V2]
  cnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix,
                                 annotations_file = anno_cnv,
                                 delim = "\t",
                                 gene_order_file = geneInfor,
                                 ref_group_names = ref_celltype)
  saveRDS(cnv_obj, 
          file = paste0("E:/Projects/BRCA/results/SC/Epithelial/infercnv/input/", spID, "_infercnv_obj.Rds"))

}


##### Load inferCNV results and calculate CNV scores #####
for(sample in sample_list){
  dir.create(paste0("D:/Projects/Breast Cancer/results/SC/Epithelial/infercnv/",sample))
}
epi$cnv_score = 99999.99

library(scales)
library(cluster)

output_dir = "E:/Projects/BRCA/results/SC/Epithelial/infercnv/"
cancer_x_threshold_sd_multiplier <- 1.25
cancer_y_threshold_sd_multiplier <- 1.25
normal_x_threshold_sd_multiplier <- 1
normal_y_threshold_sd_multiplier <- 1.25

for(sample in sample_list_sc[1:length(sample_list_sc)]){
  print(paste0(sample, " processing"))
  obs = as.data.frame(t(read.table(file = paste0(output_dir, sample, "_infercnv/infercnv.observations.txt"))))
  rownames(obs) = gsub("\\.","-", rownames(obs))

  cnv_average = apply(obs, 2, mean)
  cnv_density_plot = density(cnv_average, bw="SJ")
  # Save a density plot of the mean copy number per gene 
  png(paste0(output_dir,sample, "_infercnv/average_CNV_density_plot.png"),width = 14, height = 12, units = "cm", res = 300)
  plot(cnv_density_plot, main = NA, xlab = "CNV values")
  dev.off()
  pdf(paste0(output_dir,sample, "_infercnv/average_CNV_density_plot.pdf"),width = 14, height = 12)
  plot(cnv_density_plot, main = NA, xlab = "CNV values")
  dev.off()
  
  scaled_df = as.data.frame(rescale(as.matrix(obs), c(-1,1)))
  CNA_values <- apply(scaled_df, 1, function(y) {  # Obtain CNV values for each cell
    return(mean(y^2))
  })
  CNA_value_df <- data.frame(
    row.names = names(CNA_values),
    CNA_value = CNA_values
  )
  # Identify cells with CNV scores in the top 5%
  CNA_order <- order(CNA_value_df$CNA_value, decreasing=T)
  ordered_CNA_values  <- data.frame(
    row.names = rownames(CNA_value_df)[CNA_order],
    CNA_value = CNA_value_df[CNA_order,]
  )
  top_cancer <- head(ordered_CNA_values, nrow(ordered_CNA_values)*0.05)
  
  # Obtain the mean CNV profile of cells in the top 5%
  top_cancer_CNV_average <- apply(obs[rownames(top_cancer),], 2, mean)
  # Calculate the correlation between each cell and the mean profile of the top CNV cells
  cancer_correlations <- apply(obs, 1, function(x) {
    if (length(unique(as.numeric(x))) == 1) {
      cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
                               cor.p.value="no_CNVs_recorded")
    } else {
      cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "kendall")
      cor_result <- data.frame(cor$estimate, cor$p.value)
    }
    return(cor_result)
  })
  correlation_df <- do.call("rbind", cancer_correlations)
  epithelial_metadata = cbind(CNA_value_df, correlation_df)
  
  # Plot the density of CNV scores
  CNA_density_plot <- density(epithelial_metadata$CNA_value, bw="SJ")
  pdf(paste0(output_dir,sample, "_infercnv/CNA_density_plot.pdf"))
  plot(CNA_density_plot, main=NA, xlab = "CNA value")
  dev.off()
  png(paste0(output_dir,sample, "_infercnv/CNA_density_plot.png"),width = 14, height = 12, units = "cm", res = 300)
  plot(CNA_density_plot, main=NA, xlab = "CNA value")
  dev.off()
  
  quad_df <- data.frame(
    row.names = rownames(epithelial_metadata),
    CNA_value = epithelial_metadata$CNA_value, 
    cor.estimate = epithelial_metadata$cor.estimate
  )

  
  CNA_mean <- mean(epithelial_metadata$CNA_value)
  CNA_std_dev <- sd(epithelial_metadata$CNA_value)
  cor_mean <- mean(epithelial_metadata$cor.estimate)
  cor_std_dev <- sd(epithelial_metadata$cor.estimate)

  x_int1 <- CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev)
  y_int1 <- cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev)
  
  x_int = x_int1
  y_int = y_int1
  
  epithelial_metadata$normal_cell_call <- "uncertain"
    
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
  ] <- "normal"
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate > y_int
  ] <- "cancer"
  p <- ggplot(epithelial_metadata, 
              aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
  p <- p + geom_point()
  p <- p + scale_color_manual(values=c("black", "#74add1"), 
                              labels=c("Cancer", "Normal"))
  p <- p + xlab("Infercnv level")
  p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
  p <- p + theme(legend.title = element_blank())
  p <- p + geom_vline(xintercept = x_int)
  p <- p + geom_hline(yintercept = y_int)
  png(paste0(output_dir,sample, 
             "_infercnv/normal_call_quad_plot_mean_of_scaled_squares.png"), 
      width = 12.5, height = 10, units = "cm", res = 300)
  print(p)
  dev.off()
  pdf(paste0(output_dir,sample, 
             "_infercnv/normal_call_quad_plot_mean_of_scaled_squares.pdf"),
      width = 7, height = 6)
  print(p)
  dev.off()
  write.table(epithelial_metadata,file = paste0(output_dir, sample, "_infercnv/cnv_info.csv"),
              sep = ",")
}

epi@meta.data$malignance = as.character(epi$malignance)
table(epi$malignance)
table(epi$sample, epi$malignance)

brca$malignance = "Normal"
tumor_cellid = colnames(epi)[epi$malignance == "cancer"]
brca@meta.data$malignance[colnames(brca) %in% tumor_cellid] = "Tumor"
table(brca$malignance)

