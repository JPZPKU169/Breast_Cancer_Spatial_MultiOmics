##### Summarize hotspot modules into metaprograms #####
sample_list_hotspot = c(
  "P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11_S1","P11_S2","P12","P13","P14","P15","P16","P26","P27","P30","P32","P33",
  "P34","P35","P36","P37","P38","P39","P40_1","P40_2","P42","P43_1","P43_2","P44","P45","P46","P47","P48","P49"
)

module_delete = c(
  "P43_7","P44_8","P35_9","P45_8","P37_9","P1_9","P1_10","P2_6","P2_10","P3_5","P3_10","P3_4","P30_9","P30_11","P4_7",
  "P4_11","P6_5","P7_1","P7_3","P47_9","P9_9","P10_3","P10_6","P38_8","P11_S1_1","P11_S2_4","P48_7","P49_5","P40_1_7",
  "P40_2_9","P40_2_10","P16_4","P27_7","P42_5","P32_7","P33_5","P33_9"
)

module_delete = data.frame(module_delete)
module_delete = module_delete %>% separate(col = "module_delete", into = c("sample", "part2"), sep = "_",extra = "merge")
ordered_samples_in_paper[gsub("SM","SP",module_delete$sample)]

all_program = list()
for(spID in sample_list_hotspot){
  program = read.csv(paste0("E:/Projects/BRCA/results/SM/hotspot/", spID, "/modules_metab_danb.csv"), check.names = FALSE)
  program$Gene = as.character(program$Gene)
  program <- program %>%
    group_by(Module) %>%
    dplyr::select(Module, Gene)
  
  program_list = split(program$Gene, program$Module)
  names(program_list) = paste0(spID, "_", names(program_list))
  
  programs_saved = names(program_list)
  programs_saved = programs_saved[!(programs_saved %in% module_delete)]
  program_list = program_list[programs_saved]
  all_program = c(all_program, program_list)
}

### Calculate pairwise Jaccard indices between programs
make_jaccard_matrix = function(program){
  jaccard = lapply(program, function(p){
    p_jaccard = lapply(program, FUN= function(q){
      return(length(intersect(p,q))/length(union(p,q)))
    })
    return(p_jaccard)
  })
  jaccard = do.call(rbind,jaccard)
  return(jaccard)
}

similarity_matrix = make_jaccard_matrix(all_program)
typeof(similarity_matrix)
similarity_matrix = as.data.frame(similarity_matrix)
similarity_matrix <- apply(similarity_matrix, 2, as.numeric)
rownames(similarity_matrix) = colnames(similarity_matrix)

nmf_hc = hclust(as.dist(1-similarity_matrix), method = "complete")
png(filename = "E:/Projects/BRCA/results/SM/hotspot/metaprograms/MP_cluster_tree.png", width = 120, height = 18, units = "cm",
    res = 300)
plot(nmf_hc)
dev.off()

clusters <- cutree(nmf_hc, h = 0.95)
clusters = sort(clusters)

clusters
ordered_clusters = as.data.frame(table(clusters))
ordered_clusters = ordered_clusters[order(ordered_clusters$Freq, decreasing = T),]
ordered_clusters$cluster_name = paste0("MP_",1:nrow(ordered_clusters))
cutoff = 9
ordered_clusters$cluster_name[ordered_clusters$Freq < cutoff] = "none" 

Clust_list = list()
Clust_genes = list()
library(readxl)
metab_anno = read_excel("E:/Projects/BRCA/results/SM/identification-5um.xlsx")

for(i in 1:(length(which(ordered_clusters$cluster_name != "none")))){
  id = ordered_clusters$clusters[ordered_clusters$cluster_name == paste0("MP_",i)]
  programs = names(clusters)[clusters == id]
  Clust_list[[paste0("MP_",i)]] = programs
  matched_genes = all_program[programs]
  matched_genes = unlist(matched_genes)
  matched_genes = data.frame(table(matched_genes))
  matched_genes = matched_genes[order(matched_genes$Freq, decreasing = T),]
  matched_genes = matched_genes$matched_genes[matched_genes$Freq > floor(0.5 * length(programs))]
  print(paste0("length of MP_",i, " is: ", length(matched_genes)))
  Clust_genes[[paste0("MP_",i)]] = matched_genes
  mp_info = data.frame(gene = as.character(matched_genes))
  mp_info$Compounds = metab_anno$Compounds[match(mp_info$gene, metab_anno$`Precursor (Da)`)]
  mp_info$HMDB = metab_anno$HMDB[match(mp_info$gene, metab_anno$`Precursor (Da)`)]
  mp_info$KEGG = metab_anno$KEGG[match(mp_info$gene, metab_anno$`Precursor (Da)`)]
  write.table(matched_genes, file = paste0("E:/Projects/BRCA/results/SM/hotspot/metaprograms/MP_",i,"_genes.csv"), 
              sep = ",", row.names = FALSE)
  write.table(mp_info, file = paste0("E:/Projects/BRCA/results/SM/hotspot/metaprograms/MP_",i,"_info.csv"), 
              sep = ",", row.names = FALSE)
}


##### Plot the hotspot metaprogram heatmap #####
Clust_list_filtered = Clust_list[1:9]
Clust_genes_filtered = Clust_genes[1:9]

# The metabolites originally assigned to MP_7 were incorporated into MP_6.
Clust_list_filtered$MP_6 = c(Clust_list_filtered$MP_6, Clust_list_filtered$MP_7)
Clust_list_filtered$MP_7 = NULL
names(Clust_list_filtered) = paste0("MP_", 1:8)

ordered_program_names = unlist(Clust_list_filtered)
similarity_matrix_ordered = similarity_matrix[ordered_program_names, ordered_program_names]

sample_label = colnames(similarity_matrix_ordered)

all(sample_label %in% names(color_patient))
all(names(color_patient) %in% sample_label)

receptor_label = clinical_subtype2[sample_label]
receptor_label = factor(receptor_label, levels = c("HR+", "HER2+", "HR+HER2+","TNBC"))
receptor_color = c("HR+"="#b0b0d4", "HER2+"="#e5b35f", "HR+HER2+"="#dc7c79", "TNBC"="#73c1e4")

sample_label = factor(sample_label, 
                      levels = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11_S1","P11_S2","P12","P13","P14","P15","P16",
                                 "P26","P27","P30","P32","P33","P34","P35","P36","P37","P38","P39",
                                 "P40_1","P40_2","P42","P43","P44","P45","P46","P47","P48","P49"
                                 ))

library(ComplexHeatmap)
library(circlize)

col_an = HeatmapAnnotation(sample = sample_label, 
                           `receptor status` = receptor_label,
                           col = list(sample = color_patient,
                                      `receptor status` = receptor_color), 
                           show_annotation_name = F,
                           annotation_legend_param = list(
                             sample = list(
                               at = levels(sample_label) # Use the explicitly specified sample order
                               # If the legend is too long, add an option such as ncol = 2 to arrange entries in multiple columns
                             )
                           ))

col_fun = colorRamp2(c(0,0.1,0.2,0.3,0.4,0.5),c("white","#fbd797","#ee7d5c","#ca3d71", "#7f2485","#3a2866"))

pdf("E:/Projects/BRCA/results/SM/hotspot/metaprograms/hotspot_MP_heatmap_2609.pdf",
    height = 9, width = 11)
Heatmap(similarity_matrix_ordered,
        border = "black",
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = FALSE,
        show_row_names = F,
        heatmap_legend_param = list(
          title = "jaccard index",
          at = c(0,0.3,0.5),
          legend_height = unit(3.5, "cm")
        ),
        col = col_fun,
        top_annotation = col_an)
dev.off()
