setwd("Result4/")
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
meta_all <- read.csv("data/space_all.csv")
meta_all$id <- gsub("_GUT","",meta_all$id)
rownames(meta_all) <- meta_all$id
#######################]
otutab <- read.table("data/species.tsv",header = T)
otutab <- as.data.frame(otutab)
rownames(otutab) <- otutab$name
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
colnames(otutab) <- gsub("GLDS.564_metagenomics_metaG_I4_","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
otutab <- otutab[,colnames(otutab) %in% meta_all$id ]
######################
pcutils::trans(otutab, method = "log1") %>% head()
otutab <- otutab[rownames(otutab) %in% taxonomy2$Species,]
taxonomy2 <- taxonomy2[rownames(taxonomy2) %in% rownames(otutab),]
hebing(otutab, taxonomy2$Phylum, margin = 1, act = "sum") -> phylum
custom_colors <- c("#F9F9E6","#CCBEF3","#A7D8C9","#9CC9DE","#A1A1A1","#D8AEAC","#EFBE7E","#DEA0A4","#FFE28A","#A2A7BE")

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

p1 <- stackplot(Class, meta_all, group = "State", topN = 10) +
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
        panel.border = element_rect(color = "black", fill = NA)) 

#########tsne
library(Rtsne)
library(ggplot2)
library(ggExtra)

tsne_result <- Rtsne(t(otutab), 
                     dims = 2,          
                     perplexity = 10,    
                     theta = 0.5,        
                     max_iter = 500)    
tsne_df <- as.data.frame(tsne_result$Y)
rownames(tsne_df) <- colnames(otutab)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df  <- merge(tsne_df ,meta_all,by="row.names")
rownames(tsne_df) <- tsne_df$Row.names
p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Tissue, shape = Individual)) +
  geom_point(size = 2) +
  labs(title = "t-SNE Projection with Custom Shapes for Groups") +
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )

p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Tissue, shape = Individual)) +
  geom_point(size = 5) +
  #labs(title = "t-SNE Projection with Custom Shapes and Colors for Groups") +
  

  scale_shape_manual(values = c(16, 17, 15, 18)) +
  

  scale_color_manual(values = c("#91D1C2", "#E64B35", "#4DBBD5")) +

  theme_minimal() + 
  theme(
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )

p_with_marginal <- ggMarginal(p, type = "density",groupColour = T, groupFill = TRUE)
ggMarginal(p, type = "density", fill = "blue")

a <- rare_curve_sample(otu_spe)
  
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
#########################################################################
otu_spe <- otutab
otutab <- read.table("data/all.tsv",header = T)
otutab <- as.data.frame(otutab)
rownames(otutab) <- otutab$name
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
colnames(otutab) <- gsub("GLDS.564_metagenomics_metaG_I4_","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
otu_all  <- otutab[,colnames(otutab) %in% meta_all$id ]
#########pctils oral
metadata_oral <- meta_all[meta_all$Tissue=="Oral",]
otutab <- otu_spe[,colnames(otu_spe) %in% rownames(metadata_oral)]
b_analyse(otutab, method = "pca") -> b_res
plot(b_res, "State.1", metadata_oral, bi = TRUE, rate = 0.9, sample_label = FALSE)
p <- pctax::plot(b_res, "State.1", metadata_oral, bi = T, rate = 0.5) + geom_text_repel(size = 1, max.overlaps = 50)
##########a
a_diversity(otutab) -> a_res
plot(a_res, "State.1", metadata_oral)+ 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),    
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
dist1 <- mat_dist(otutab, method = "bray")
as.b_dist(dist1, group_df = metadata_oral["State.1"]) -> b_dist1
plot(b_dist1, c_group = "intra", alpha = T) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),   
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )
plot(b_dist1, mode = 2)
###################
###########Desq2
otutab <- otu_all
meta_need <- meta_all[meta_all$Group==1 |meta_all$Group==3,]
otutab <- otu_all[,colnames(otu_all) %in% rownames(meta_need)]
######

meta_gut<- meta_need[meta_need$Tissue=="Gut",]
otutab_gut <- otutab[,colnames(otutab) %in% rownames(meta_gut)]
otutab_skin <- meta_gut
meta_skin<- meta_need[meta_need$Tissue=="Skin",]
otutab_skin <- otutab[,colnames(otutab) %in% rownames(meta_skin)]

meta_oral<- meta_need[meta_need$Tissue=="Oral",]
otutab_oral <- otutab[,colnames(otutab) %in% rownames(meta_oral)]
##################
library(DESeq2)
DES <- otutab_gut
col<- data.frame(row.names = colnames(DES),Group=meta_gut$State.1)
col$Group <- factor(col$Group)
levels(col$Group)
col$Group <- relevel(col$Group, ref = "Pre Flight")
dds <- DESeqDataSetFromMatrix(DES, col, design= ~ Group)
dds <- DESeq(dds)
des <- DES[rowMeans(DES)>2,] 
des <-des[complete.cases(des),]
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 8, parallel = FALSE) 
resgut <- results(dds1)
res1 <- data.frame(resgut, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
library(ggplot2)
library(ggrepel)
df4_label <- df4[df4$pvalue < 1.4e-6, ]
ggplot(df4, aes(logFC, -log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(pvalue), color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0, 1, 0.1),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.01, 0.7),
        legend.justification = c(0, 1)) +
  guides(col = guide_colourbar(title = "-Log10(p-value)"),
         size = "none") +
  geom_text_repel(data = df4_label, aes(label = df4_label[,1]), color = "black", size = 3) +
  xlab("logFC") +
  ylab("-Log10(p-value)")

################
DES <- otutab_oral

col<- data.frame(row.names = colnames(DES),Group=meta_oral$State.1)
col$Group <- factor(col$Group)
levels(col$Group)
col$Group <- relevel(col$Group, ref = "Pre Flight")
dds <- DESeqDataSetFromMatrix(DES, col, design= ~ Group)
dds <- DESeq(dds)
des <- DES[rowMeans(DES)>2,] 
des <-des[complete.cases(des),]
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 16, parallel = FALSE) 
resoral <- results(dds1)
res1 <- data.frame(resoral, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
library(ggplot2)
library(ggrepel)
df4_label <- df4[df4$pvalue < 1.4e-5, ]
ggplot(df4, aes(logFC , -log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(pvalue), color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0, 1, 0.1),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.01, 0.7),
        legend.justification = c(0, 1)) +
  guides(col = guide_colourbar(title = "-Log10(p-value)"),
         size = "none") +
  geom_text_repel(data = df4_label, aes(label = df4_label[,1]), color = "black", size = 3) +
  xlab("logFC") +
  ylab("-Log10(p-value)")
################
DES <- otutab_skin

col<- data.frame(row.names = colnames(DES),Group=meta_skin$State.1)
col$Group <- factor(col$Group)
levels(col$Group)
col$Group <- relevel(col$Group, ref = "Pre Flight")
dds <- DESeqDataSetFromMatrix(DES, col, design= ~ Group)
dds <- DESeq(dds)
des <- DES[rowMeans(DES)>2,] 
des <-des[complete.cases(des),]
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 8, parallel = FALSE) 
resskin <- results(dds1)
res1 <- data.frame(resskin, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
library(ggplot2)
library(ggrepel)
df4_label <- df4[df4$pvalue < 1e-7, ]
ggplot(df4, aes(logFC, -log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(pvalue), color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0, 1, 0.1),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.01, 0.7),
        legend.justification = c(0, 1)) +
  guides(col = guide_colourbar(title = "-Log10(p-value)"),
         size = "none") +
  geom_text_repel(data = df4_label, aes(label = df4_label[,1]), color = "black", size = 3) +
  xlab("logFC") +
  ylab("-Log10(p-value)")
##############

ress <- list(resoral,resgut,resskin)
names(ress) <- c("Oral", "Gut", "Skin")
library(volcano3D)
obj <- deseq_2x3_polar(ress,padj.method = "none",fc_cutoff = 0.01,pcutoff = 0.05,process = "two.side")

volcano3D(obj, type = 1, marker_size = 4,label_rows = NULL)
plot <- volcano3D(obj, type = 1, marker_size = 4, label_rows = NULL)
plot <- plot %>% layout(showlegend = FALSE)

plot <- volcano3D1(obj, type = 1, marker_size = 4)
######################
res1 <- data.frame(resskin, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
df4_skin <- df4
res1 <- data.frame(resgut, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
df4_gut <- df4

res1 <- data.frame(resoral, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
df4_oral <- df4
df4_oral$Tissue <- "oral"
df4_skin$Tissue <- "skin"
df4_gut$Tissue <- "gut"
df4_all <- rbind(df4_gut,df4_skin,df4_oral)
write.csv(df4_all,"space_tissue_deg.csv")
########
