setwd("Result3/")
library(lme4)
library(cluster)
library(pctax)
library(data.table)
library(dplyr)
library(pcutils)
library(devtools)
library(ggplot2)
library(ggpmisc)
library(ape)
library(ggsci)
library(ggGenshin)
library(dplyr)

################metadata
meta_all <- read.csv("all/metadata_all.csv")
rownames(meta_all) <- meta_all$Sample
#######################
otutab <- read.table("all/species.tsv",header = T)
rownames(otutab) <- otutab$name
index <- otutab[,c(1:4)]
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
otutab <- otutab[,colnames(otutab) %in% meta_all$Sample ]
otutab_all <- otutab
otutab_spe<- otutab

######################
otutab <- otutab_spe
pcutils::trans(otutab, method = "log1") %>% head()
otutab <- otutab[rownames(otutab) %in% taxonomy2$Species,]
taxonomy2 <- taxonomy2[rownames(taxonomy2) %in% rownames(otutab),]
hebing(otutab, taxonomy2$Phylum, margin = 1, act = "sum") -> phylum
write.csv(phylum,"phylum.csv")
custom_colors <- c("#F9F9E6","#CCBEF3","#A7D8C9","#9CC9DE","#A1A1A1","#D8AEAC","#EFBE7E","#DEA0A4","#FFE28A","#A2A7BE")
meta_all <- meta_all[rownames(meta_all) %in% colnames(otutab),]
p1 <- stackplot(phylum, meta_all, group = "Tissue", topN = 10) +
  scale_fill_manual(values = custom_colors) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
hebing(otutab, taxonomy2$Class, margin = 1, act = "sum") -> Class
custom_colors <- c("#F9F9E6","#CCBEF3","#A7D8C9","#9CC9DE","#A1A1A1","#D8AEAC","#EFBE7E","#DEA0A4","#FFE28A","#A2A7BE")
p1 <- stackplot(Class, meta_all, group = "Tissue", topN = 10) +
  scale_fill_manual(values = custom_colors) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
stackplot(phylum, meta_all,
          group = "Tissue", style = "sample",
          group_order = F, flow = F,topN = 10,
) +
  scale_fill_manual(values = custom_colors)+
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_blank(),    
        axis.ticks.x = element_blank()) 
#################################
install.packages("umap")
library(umap)
umap_config <- umap.defaults
umap_config$n_neighbors <- 8  
umap_config$min_dist <- 0.3    
umap_config$n_components <- 2  
umap_config$metric <- "euclidean"  

result <- umap(t(otutab), config = umap_config)
umap_df <- as.data.frame(result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df <- merge(umap_df,meta_all,by="row.names")
rownames(umap_df ) <- umap_df$Row.names

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Tissue)) +
  geom_point(size = 2) +
  labs(title = "UMAP Projection") +
  theme_minimal()

#########tsne
library(Rtsne)
library(ggplot2)
library(ggExtra)
tsne_result <- Rtsne(t(otutab), 
                     dims = 2,           
                     perplexity = 20,    
                     theta = 0.5,       
                     max_iter = 500)    
tsne_df <- as.data.frame(tsne_result$Y)
rownames(tsne_df) <- colnames(otutab)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df  <- merge(tsne_df ,meta_all,by="row.names")
rownames(tsne_df) <- tsne_df$Row.names

p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Tissue, shape = Gender)) +
  geom_point(size = 2) +
  labs(title = "t-SNE Projection with Custom Shapes for Groups") +
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )
p
p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Tissue, shape = Gender)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_color_manual(values = c("#91D1C2", "#E64B35", "#4DBBD5")) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )
p_with_marginal <- ggMarginal(p, type = "density",groupColour = T, groupFill = TRUE)
p_with_marginal
ggMarginal(p, type = "density", fill = "blue")

###################################################
a <- rare_curve_sample(otutab_spe)
p <- plot(a)+  theme_minimal() +
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )
a <- rare_curve_species(otutab, mode = 1) 
p <- plot(a)+ 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),    
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )
plot(a)
group_box(meta_all["Age"], group = "Group", meta_all, alpha = T)
#########################################################################
otu_spe <- otutab

#########pctils oral
metadata_oral <- meta_all[meta_all$tissue=="Oral",]
otutab <- otutab_spe[,colnames(otutab_spe) %in% rownames(metadata_oral)]
a_diversity(otutab) -> a_res
plot(a_res, "Age", metadata_oral) +   theme_minimal() + 
  theme(
    panel.grid = element_blank(),    
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )
as.b_dist(dist1, group_df = metadata_oral["Age"]) -> b_dist1
plot(b_dist1, c_group = "intra", alpha = T) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),    
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
plot(b_dist1, mode = 2)+ 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )
metadata_skin <- meta_all[meta_all$tissue =="Skin",]

otutab <- otutab_spe[,colnames(otutab_spe) %in% rownames(metadata_skin)]
a_diversity(otutab) -> a_res
plot(a_res, "Age", metadata_skin) +   theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1)  #
  )
dist1 <- mat_dist(otutab, method = "bray")


as.b_dist(dist1, group_df = metadata_skin["Age"]) -> b_dist1
plot(b_dist1, mode = 2)+ 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )

metadata_Gut <- meta_all[meta_all$tissue=="Gut",]
otutab <- otutab_spe[,colnames(otutab_spe) %in% rownames(metadata_Gut)]
a_diversity(otutab) -> a_res
plot(a_res, "Age", metadata_Gut ) +   theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1)  #
  )
dist1 <- mat_dist(otutab, method = "bray")

as.b_dist(dist1, group_df = metadata_Gut["Age"]) -> b_dist1
plot(b_dist1, mode = 2)+  
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )



rownames(meta_gut) <- meta_gut$Sample
dist1 <- mat_dist(otutab, method = "bray")
as.b_dist(dist1, group_df = meta_gut["Age"]) -> b_dist1
plot(b_dist1, c_group = "intra", alpha = T)
plot(b_dist1, mode = 2) + theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )


