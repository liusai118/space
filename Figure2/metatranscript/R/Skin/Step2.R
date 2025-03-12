library(DESeq2)
otutab <- read.table("all_all_transcript.tsv",header = T)
otutab <- as.data.frame(otutab)
rownames(otutab) <- otutab$name
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
metadata <- meta_skin[meta_skin$Group=="1"|meta_skin$Group=="3",]
DES<- otutab[,colnames(otutab) %in% metadata$id ]
col<- data.frame(row.names = colnames(DES),Group=metadata$State.1)
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
df4_label <- df4[df4$pvalue < 1.4e-8, ]
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
  xlab("logFC") +
  ylab("-Log10(p-value)")

