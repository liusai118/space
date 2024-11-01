library(DESeq2)
DES <- transcript
metadata <- meta_skin[meta_skin$Group=="1"|meta_skin$Group=="3",]
metadata <- metadata[rownames(metadata) %in% colnames(DES),]
DES <- DES[,rownames(metadata)] 
col<- data.frame(row.names = colnames(DES),Group=metadata$State.1)
# 确保 Group 是因子类型
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
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  # 散点图:
  geom_point(aes(size = -log10(pvalue), color = -log10(pvalue))) +
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0, 1, 0.1),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1, 3)) +
  # 主题调整：
  theme_bw() +
  # 调整主题和图例位置：
  theme(panel.grid = element_blank(),
        legend.position = c(0.01, 0.7),
        legend.justification = c(0, 1)) +
  # 设置部分图例不显示：
  guides(col = guide_colourbar(title = "-Log10(p-value)"),
         size = "none") +
  # 添加标签，使用 geom_text_repel 确保标签不会重叠和超出边界：
  geom_text_repel(data = df4_label, aes(label = df4_label[,1]), color = "black", size = 3) +
  # 修改坐标轴：
  xlab("logFC") +
  ylab("-Log10(p-value)")

