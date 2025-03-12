library(ReporterScore)
KO_htable <- load_KO_htable()
head(KO_htable, n = 20)
meta_skin <- meta_all[meta_all$Tissue=="Skin",]
###############
rownames(meta_skin) <- gsub("_TZO","",rownames(meta_skin) )
KO_abundance <- read.table("eggnog.KO.raw.txt",header = T,row.names = 1)
colnames(KO_abundance) <- gsub("_TZO.quant","",colnames(KO_abundance))
KOlist <- load_KOlist()
head(KOlist$pathway)

cat("Comparison: ", levels(factor(metadata$State.1)), "\n")
#> Comparison:  WT OE
# for microbiome!!!
reporter_res <- GRSA(KO_abundance, "State.1", metadata,
                     mode = "directed",
                     method = "t.test", perm = 999,
                     type = "pathway", feature = "ko",
)
plot_report_bar(reporter_res, rs_threshold = c(-2, 2), facet_level = TRUE)
plot_report_circle_packing(reporter_res, rs_threshold = c(-2.5, 2.5))
plot_KOs_in_pathway(reporter_res, map_id = "map00780")
plot_KOs_distribution(reporter_res, map_id = "map00780")
plot_KOs_network(reporter_res,
                 map_id = c("map00780", "map00785", "map00900"),
                 main = "", mark_module = TRUE
)
a <- reporter_res$modulelist
plot_KOs_heatmap(reporter_res,
                 map_id =a$id , only_sig = TRUE,
                 heatmap_param = list(cutree_rows = 1)
)

library(pheatmap)
library(pheatmap)
library(RColorBrewer)
b <- na.omit(KO_abundance)
b <- b[,-1]
b <- scale(b)

meta_skin2 <- data.frame(
  State = factor(meta_skin$State),
  Individual= factor(meta_skin$Individual),
  Gender=factor(meta_skin$Gender)
)
npg_colors <- pal_npg("nrc")(8)
vav <- pal_jama("default")(4)
va2 <- pal_d3("category10")(2)

annotation_colors = list(
  State = setNames(npg_colors, levels(meta_skin2$State)),
  Individual = setNames(vav, levels(meta_skin2$Individual)),
  Gender=setNames(va2, levels(meta_skin2$Gender))
)
rownames(meta_skin2) <- rownames(meta_skin)

p1 <- pheatmap(b,  
               color = colorRampPalette(c('blue','white','#FF4500'))(100),
               border_color = "black",  
               scale = "row", 
               cluster_rows =T, 
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","high"), 
               show_rownames = FALSE, 
               show_colnames = TRUE, 
               kmeans_k = 2000,
               #border_color = "white",
               clustering_method = "ward.D",
               annotation_col = meta_skin2,  
               annotation_colors = annotation_colors,
               annotation_names_row = FALSE,  
               annotation_names_col = FALSE, 
               fontsize = 8 
) 
