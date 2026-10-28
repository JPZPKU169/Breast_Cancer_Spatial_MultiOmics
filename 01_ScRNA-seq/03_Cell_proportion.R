library(ggplot2)
library(Seurat)
########### Stacked bar plot ##########
make_proportion_stackplot = function(object, group.by, sample_info, sample_level, celltype_info, celltype_level, remove_type="unknown", 
                                     colors){
  meta = object@meta.data[, c(celltype_info,sample_info,group.by)]
  meta = meta[!(meta[,celltype_info] %in% remove_type),]
  celltype_level = celltype_level[!(celltype_level %in% remove_type)]
  meta[,celltype_info] = factor(meta[,celltype_info], levels = celltype_level)
  meta[,sample_info] = factor(meta[,sample_info], levels = sample_level)
  meta = meta[order(meta[,sample_info]),]
  colors = colors[1:length(celltype_level)]
  
  sample_x = 0:(length(sample_level)-1)
  step = 0
  for(i in 2:length(sample_x)){
    if(meta[match(sample_level[i], meta[,sample_info]),group.by] != meta[match(sample_level[i-1], meta[,sample_info]),group.by]){
      step = step + 0.1
    }
    sample_x[i] = sample_x[i] + step
  }
  names(sample_x) = sample_level
  
  sample_celltype_prop <- as.data.frame(prop.table(table(meta[,celltype_info], meta[,sample_info]),2))
  colnames(sample_celltype_prop) = c("celltype" ,"sample", "Freq")
  sample_celltype_prop$x = sample_x[sample_celltype_prop$sample]
  
  ggplot(sample_celltype_prop, aes(x = x, y = Freq, fill = celltype)) +
    scale_fill_manual(values =colors) +
    geom_bar(stat = "identity", position = "stack") + ylab("Cell Proportion") +
    scale_x_continuous(limits = c(-0.5,max(sample_x)+1),breaks = sample_x, labels = sample_level) +
    theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 5),
          axis.line.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour = "black", size = 6),
          axis.ticks.length = unit(1.5, "pt"),
          axis.ticks = element_line(colour = "black", linewidth = 0.5 * 0.3528),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", linewidth = 0.5 * 0.3528),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.key.size = unit(6, "pt"),
          legend.background = element_rect(linewidth = 0.5 * 0.3528))
  
}


# All celltype
celltype_colors = c("#f4a261","#2a9d8f","#adc178","#f5cac3","#7a529f","#E64B35CC","#ef476f","#ae9a84","#71c9dd")
make_proportion_stackplot(brca, group.by = "receptor_status", sample_info = "patient_id", sample_level = sc_sample_level, 
                          celltype_info = "celltype", celltype_level = celltype_levels, remove_type="unknown", 
                          colors = celltype_colors)
ggsave(filename = "E:/Projects/BRCA/results/SC/proportion_stack_all_2609.png", width = 13, height = 4, units = "cm", dpi = 300)
ggsave(filename = "E:/Projects/BRCA/results/SC/proportion_stack_all_2609.pdf", width = 13, height = 4, units = "cm", dpi = 300)
