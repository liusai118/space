setwd("Result3/")
###############Gut
meta <- read.csv("gut/meta.csv")
matrix <- read.table("gut/all.tsv",header = T)
rownames(matrix) <- matrix[,1]
m <- matrix[,-c(1:4)]
m <- m[-1,]
colnames(m) <- gsub(".taxonReads","",colnames(m))
m <- m[,colnames(m) %in% meta$ID]
m[is.na(m)] <- 0
meta <- meta[meta$ID %in%colnames(m),]
meta_yong <- meta %>%
  dplyr::filter(host.age <36)
meta_old <- meta %>%
  dplyr::filter(host.age > 57)
data_yong <- m[,colnames(m) %in% meta_yong$ID]
need <- c("Pseudomonadota","Bacillota","Euryarchaeota","Actinomycetota","Thermodesulfobacteriota")
data_old2 <- data_old[rownames(data_old) %in% need,]
data_yong2 <- data_yong[rownames(data_yong) %in% need,]
data_old <-  m[,colnames(m) %in% meta_old$ID]
dataolds <- colSums(data_old)
datayoung <- colSums(data_yong)
datao <- data_old
for (i in 1:ncol(data_old)){
  datao[,i] <- data_old[,i] /dataolds[i] * 100
}
datay <- data_yong
for (i in 1:ncol(data_yong)){
  datay[,i] <- data_yong[,i] /datayoung[i] * 100
}

col <- c(rep("young",27),rep("old",24))
colnames(data_old2) <- "old"
colnames(data_yong2) <- "young"
data <- cbind(datay,datao)
colnames(data) <- col
data <- data[rownames(data) %in% need,]

library(DESeq2)

DES <- cbind(data_yong,data_old)
need <- c("Megamonas funiformis","Megamonas","Selenomonadaceae","Selenomonadales","Negativicutes","Bacillota","Bacteria")
DES_need <- DES[rownames(DES) %in% need,]
sum <- rowSums(DES)
DES[,1] <- sum
col<- data.framDEScol<- data.frame( x = colnames(DES),condition = c(rep("young",27),rep("old",24)) )
dds <- DESeqDataSetFromMatrix(DES, col, design= ~ condition)
dds <- DESeq(dds)
des <- DES[rowMeans(DES)>2,] 
des <-des[complete.cases(des),]
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 27, parallel = FALSE) 
res_gut <- results(dds1)

res1 <- data.frame(res_gut, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
#############################################
df4_label <- df4[df4$pvalue < 1.4e-17, ]
library(ggplot2)
library(ggrepel)
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
  xlab("LogFC") +
  ylab("-Log10(p-value)")
###########all

meta_all <- read.csv("all/metadata_all.csv")
rownames(meta_all) <- meta_all$Sample

#######################
otutab <- read.table("all/all.tsv",header = T)
rownames(otutab) <- otutab$name
index <- otutab[,c(1:4)]
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
otutab <- otutab[,colnames(otutab) %in% meta_all$Sample ]
otutab_all <- otutab
#otutab <- otutab_all
metadata_oral <- meta_all[meta_all$tissue=="Oral",]
metadata_skin <- meta_all[meta_all$tissue=="Skin",]
meta_oral <- meta_all[meta_all$tissue=="Oral",]
meta_skin <- meta_all[meta_all$tissue=="Skin",]
young <- meta_skin[meta_skin$Age < 35,]
young$Group <- "young"
aging <- meta_skin[meta_skin$Age > 55,]
aging$Group <- "old"
meta_skin1 <- rbind(young,aging)
otutab_skin <- otutab_all[,colnames(otutab_all) %in% meta_skin1$Sample]
rownames(meta_skin1) <- meta_skin1$Sample
otutab_oral <- otutab_all[,colnames(otutab_all) %in% meta_oral1$Sample]
rownames(meta_skin1) <- meta_skin1$Sample
otutab <- otutab_skin
b_analyse(otutab, method = "pca") -> b_res
plot(b_res, "Group", meta_skin1, bi = T,sample_label = FALSE)
DES <- otutab_skin
col<- data.frame(row.names = colnames(DES),Group=meta_skin1$Group)
col$Group <- factor(col$Group)
levels(col$Group)
col$Group <- relevel(col$Group, ref = "young")
dds <- DESeqDataSetFromMatrix(DES, col, design= ~ Group)
dds <- DESeq(dds)
des <- DES[rowMeans(otutab_skin)>2,] 
des <-des[complete.cases(des),]
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 16, parallel = FALSE) 
resskin <- results(dds1)
res1 <- data.frame(resskin, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
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
  xlab("LogFC") +
  ylab("-Log10(p-value)")

##############

young <- meta_oral[meta_oral$Age < 31,]
young$Group <- "young"
aging <- meta_oral[meta_oral$Age > 54,]
aging$Group <- "old"

meta_oral1 <- rbind(young,aging)
meta_oral1 <- meta_oral1 [meta_oral1$Sample!= "SRR2228235",]
meta_oral1 <- meta_oral1 [meta_oral1$Sample!= "SRR2228548",]
otutab_oral <- otutab_all[,colnames(otutab_all) %in% meta_oral1$Sample]

rownames(meta_oral1) <- meta_oral1$Sample
otutab <- otutab_oral
b_analyse(otutab, method = "pca") -> b_res
plot(b_res, "Group", meta_oral1, bi = T,sample_label = FALSE)
plot(b_res, "Group", meta_oral1, mode = 3)
DES <- otutab_oral
meta_oral1 <- meta_oral1[rownames(meta_oral1) %in% colnames(DES),]
col<- data.frame(row.names = colnames(DES),Group=meta_oral1$Group)
# 确保 Group 是因子类型
col$Group <- factor(col$Group)
levels(col$Group)
col$Group <- relevel(col$Group, ref = "young")
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
df4_label <- df4[df4$pvalue < 1.4e-20, ]
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
young <- meta_skin[meta_skin$Age < 35,]
young$Group <- "young"
aging <- meta_skin[meta_skin$Age > 55,]
aging$Group <- "old"
meta_skin1 <- rbind(young,aging)
otutab_skin <- otutab_all[,colnames(otutab_all) %in% meta_skin1$Sample]
rownames(meta_skin1) <- meta_skin1$Sample
otutab_oral <- otutab_all[,colnames(otutab_all) %in% meta_oral1$Sample]
rownames(meta_skin1) <- meta_skin1$Sample
otutab <- otutab_skin
b_analyse(otutab, method = "pca") -> b_res
plot(b_res, "Group", meta_skin1, bi = T,sample_label = FALSE)
DES <- otutab_skin
col<- data.frame(row.names = colnames(DES),Group=meta_skin1$Group)
col$Group <- factor(col$Group)
levels(col$Group)
col$Group <- relevel(col$Group, ref = "young")
dds <- DESeqDataSetFromMatrix(DES, col, design= ~ Group)
dds <- DESeq(dds)
des <- DES[rowMeans(otutab_skin)>2,] 
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
  xlab("LogFC") +
  ylab("-Log10(p-value)")
##############

ress <- list(resoral,resgut,resskin)
names(ress) <- c("Oral", "Gut", "Skin")
library(volcano3D)
obj <- deseq_2x3_polar(ress,padj.method = "none",fc_cutoff = 0.01,pcutoff = 0.05,process = "two.side")
volcano3D(obj, type = 1, marker_size = 4,label_rows = NULL)
plot <- volcano3D(obj, type = 1, marker_size = 4, label_rows = NULL)
plot <- volcano3D1(obj, type = 1, marker_size = 4, label_rows = NULL)
plot <- plot %>% layout(showlegend = FALSE)
plot <- volcano3D1(obj, type = 1, marker_size = 4)
######################


res1 <- data.frame(resoral, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which( res1$padj < 0.05),]     
res1_down<- res1[which(res1$padj < 0.05),]   
res1_total <- rbind(res1_up,res1_down)
des <- otutab_oral[rowMeans(otutab_oral)>2,] 
df4<-merge(res1,des, by = "row.names", all = T)
df4 <- na.omit(df4)
df4$logFC <- df4$log2FoldChange * log10(2)
df4_oral <- df4
df4_oral <-df4_oral[,1:7]
df4_skin <-df4_skin[,1:7]
df4_gut  <- df4_gut[,1:7]
df4_oral$Tissue <- "oral"
df4_skin$Tissue <- "skin"
df4_gut$Tissue <- "gut"
df4_all <- rbind(df4_gut,df4_skin,df4_oral)
write.csv(df4_all,"aging_tissue_deg.csv")
