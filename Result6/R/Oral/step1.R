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
################
meta_all <- read.csv("space_all.csv")
meta_all$id <- gsub("_GUT","",meta_all$id)
rownames(meta_all) <- meta_all$id
###############
otutab <- read.table("all.tsv",header = T)
otutab <- as.data.frame(otutab)
rownames(otutab) <- otutab$name
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
colnames(otutab) <- gsub("GLDS.564_metagenomics_metaG_I4_","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
meta_oral <- meta_all[meta_all$Tissue=="Oral",]
otutab <- otutab[,colnames(otutab) %in% meta_all$id ]
genomics <- otutab
###########################
otutab <- read.table("all_speceis_metatranscript.tsv",header = T)
otutab <- as.data.frame(otutab)
rownames(otutab) <- otutab$name
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
meta_oral <- meta_all[meta_all$Tissue=="Oral",]
otutab <- otutab[,colnames(otutab) %in% meta_oral$id ]
transcript <- otutab
pcutils::trans(otutab, method = "log1") %>% head()
otutab <- otutab[rownames(otutab) %in% taxonomy2$Species,]
taxonomy2 <- taxonomy2[rownames(taxonomy2) %in% rownames(otutab),]
hebing(otutab, taxonomy2$Species, margin = 1, act = "sum") -> phylum
write.csv(phylum,"phylum_spe.csv")
p1 <- stackplot(phylum, meta_oral, group = "State", topN = 10,style = "Sample") +
  #scale_fill_manual(values = custom_colors) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

stackplot(phylum, meta_oral,
          group = "Tissue", style = "sample",
          group_order = F, flow = F,topN = 10,
) +
  scale_fill_manual(values = pal_d3("category10")(10))+
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) 
p1 + scale_color_npg()  # 提供10种颜色
common <- intersect(rownames(genomics),rownames(transcript))
common_genomic <- genomics[common,]
total_reads1 <- colSums(common_genomic)
relative_common_genomic <- sweep(common_genomic, 2, total_reads1, FUN = "/")
common_getrans <- transcript[common,]
total_reads <- colSums(common_getrans)
relative_common_trans <- sweep(common_getrans, 2, total_reads, FUN = "/")
common2 <- relative_common_trans[,meta_oral$id] /relative_common_genomic[,meta_oral$id]
mean_gen <- rowMeans(common_genomic)
mean_tra <- rowMeans(common_getrans)
mean_gen2 <- rowMeans(relative_common_genomic)
mean_tra2 <- rowMeans(relative_common_trans)

data <- data.frame(gen=mean_gen,tra=mean_tra)
data2 <- data.frame(gen=log(mean_gen2),tra= log(mean_tra2))
data2 <- data2 %>%
  dplyr::filter(gen > -20)
data2 <- data2 %>%
  dplyr::filter(tra> -20)
###########


ggplot(data2, aes(x = gen, y = tra)) + 
  geom_point(color = "#E64B35", size = 1.5, alpha = 0.5) +  # 调整散点透明度
  geom_smooth(method = "lm", color = "black", linetype = 4, se = TRUE, level = 0.9999999999999999) +  # 添加虚线回归线，带95%置信区间
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
               formula = y ~ x, 
               parse = TRUE) +  # 添加回归方程、R²和p值
  ggtitle("Scatter Plot with Regression Line, p-value, and R²") + 
  xlab("log10(Relative abundance)") + 
  ylab("log10(Relative contribution to the MT)") + 
  theme_minimal() +  # 使用简洁主题
  theme(panel.grid = element_blank(),  # 删除背景网格线
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # 添加黑色边框


ggplot(data2, aes(x = gen, y = tra)) + 
  geom_point(color = "#E64B35", size = 1.5, alpha = 0.5) +  # 调整散点透明度
  geom_smooth(method = "lm", color = "black", linetype = 4, se = TRUE, level = 0.9999999999999999) +  # 添加虚线回归线，带95%置信区间
  theme_minimal() +  # 使用简洁主题
  theme(panel.grid = element_blank(),  # 删除背景网格线
        panel.border = element_rect(color = "black", fill = NA, size = 1,
                                    ),
        axis.text.x = element_blank(),  # 移除 x 轴刻度标签
        axis.text.y = element_blank(),  # 移除 y 轴刻度标签
        legend.position = "none")  # 移除图例)  # 添加黑色边框

ggplot(data2, aes(x = gen, y = tra)) + 
  geom_point(color = "#E64B35", size = 1.5, alpha = 0.5) +  # 调整散点透明度
  geom_smooth(method = "lm", color = "black", linetype = 0,,se = TRUE, level = 0.9999999999999999) +  # 添加虚线回归线，带95%置信区间
  theme_minimal() +  # 使用简洁主题
  theme(panel.grid = element_blank(),  # 删除背景网格线
        panel.border = element_rect(color = "white", fill = NA, size = 1,
        ),
        axis.text.x = element_blank(),  # 移除 x 轴刻度标签
        axis.text.y = element_blank(),  # 移除 y 轴刻度标签
        legend.position = "none")  +# 移除图例)  # 添加黑色边框
  labs(x = NULL, y = NULL)  # 移除横坐标和纵坐标的标题

##########
# 使用 dplyr 的 mutate_all 来替换非数值为 0
common2 <- common2 %>%
  mutate_all(~ ifelse(is.na(as.numeric(.)), 0, as.numeric(.)))
common2<- common2 %>%
  mutate_all(~ ifelse(is.infinite(.), 0, .))

rowsums <- rowSums(common2)
common21<-  common2[order(rowsums,decreasing = T),]
common21<- common21[1:20,]
common21<- t(common21)
common21 <- log(common21 +1)
data_long <- melt(common21, varnames = c("Sample", "Species"), value.name = "Expression")

# 使用ggplot2绘制小提琴图
# 绘制小提琴图，去除小提琴的边框
ggplot(data_long, aes(x = Expression, y = Species, fill = Species)) + 
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +  # 去除小提琴图的边框
  geom_boxplot(width = 0.1, position = position_dodge(0.9), fill = "white", 
               outlier.size = 1, lwd = 0.5) +  # 调整箱线图的线条粗细
  theme_minimal() +  # 使用简洁主题
  theme(
    panel.background = element_blank(),               # 去除背景填充
    panel.grid = element_blank(),                     # 去除背景网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    axis.line = element_line(color = "black")         # 添加X轴和Y轴线
  ) +
  theme(legend.position = "none") +  # 移除图例
  xlab("Relative Expression (MT/MG)") +  # 设置 X 轴标签
  ylab("") +  # 移除 Y 轴标签
  scale_fill_manual(values = pal_d3("category20")(20)) +  # 自动生成并设置颜色
  scale_x_log10() +  # 对 X 轴进行对数缩放
  ggtitle("Violin Plot of Species Expression Across Samples")  # 设置标题
# 绘制小提琴图
ggplot(data_long, aes(x = Expression, y = Species, fill = Species)) + 
  geom_violin(trim = FALSE, alpha = 0.7) +  # 绘制小提琴图
  geom_boxplot(width = 0.1, position = position_dodge(0.9), fill = "white", outlier.size = 1) +  # 添加箱线图
  theme_minimal() +  # 使用简洁主题
  theme(
    panel.background = element_blank(),               # 去除背景填充
    panel.grid = element_blank(),                     # 去除背景网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    axis.line = element_line(color = "black")         # 添加X轴和Y轴线
  ) +
  theme(legend.position = "none") +  # 移除图例
  xlab("Relative Expression (MT/MG)") +  # 设置 X 轴标签
  ylab("") +  # 移除 Y 轴标签
  scale_fill_manual(values = pal_d3("category20")(20)) +  # 自动生成并设置颜色
  scale_x_log10() +  # 对 X 轴进行对数缩放
  ggtitle("Violin Plot of Species Expression Across Samples")  # 设置标题
