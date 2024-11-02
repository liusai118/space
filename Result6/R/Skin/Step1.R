setwd("Result6/R/SKin/")
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
###############
meta_all <- read.csv("space_all.csv")
meta_all$id <- gsub("_GUT","",meta_all$id)
rownames(meta_all) <- meta_all$id
###############
otutab <- read.table("species.tsv",header = T)
otutab <- as.data.frame(otutab)
rownames(otutab) <- otutab$name
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
colnames(otutab) <- gsub("GLDS.564_metagenomics_metaG_I4_","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
meta_skin <- meta_all[meta_all$Tissue=="Skin",]
otutab <- otutab[,colnames(otutab) %in% meta_skin$id ]
genomics <- otutab
###########################
otutab <- read.table("skin_meta_transcript_all.tsv",header = T)
otutab <- as.data.frame(otutab)
rownames(otutab) <- otutab$name
otutab <-otutab[,-c(1:4)]
colnames(otutab) <- gsub(".cladeReads","",colnames(otutab))
otutab <- otutab %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

meta_skin <- meta_all[meta_all$Tissue=="Skin",]
rownames(meta_skin)  <- gsub("_TZO","", rownames(meta_skin) )
otutab <- otutab[,colnames(otutab) %in% rownames(meta_skin) ]
transcript <- otutab
meta_skin <- meta_skin[ rownames(meta_skin)  %in% colnames(transcript),]




pcutils::trans(otutab, method = "log1") %>% head()

otutab <- otutab[rownames(otutab) %in% taxonomy2$Species,]
taxonomy2 <- taxonomy2[rownames(taxonomy2) %in% rownames(otutab),]

hebing(otutab, taxonomy2$Species, margin = 1, act = "sum") -> Species
stackplot(Species, meta_skin,
          group = "Tissue", style = "sample",
          group_order = F, flow = F,topN = 10,
) +
  scale_fill_manual(values = pal_d3("category10")(10))+
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) 

common <- intersect(rownames(genomics),rownames(transcript))
common_genomic <- genomics[common,]
total_reads1 <- colSums(common_genomic)
relative_common_genomic <- sweep(common_genomic, 2, total_reads1, FUN = "/")
common_getrans <- transcript[common,]
total_reads <- colSums(common_getrans)
relative_common_trans <- sweep(common_getrans, 2, total_reads, FUN = "/")
common2 <- relative_common_trans[,rownames(meta_skin)] /relative_common_genomic[,meta_skin$id]


###############
common <- intersect(rownames(genomics),rownames(transcript))
common_genomic <- genomics[common,]
total_reads1 <- colSums(common_genomic)
relative_common_genomic <- sweep(common_genomic, 2, total_reads1, FUN = "/")
common_getrans <- transcript[common,]
total_reads <- colSums(common_getrans)
relative_common_trans <- sweep(common_getrans, 2, total_reads, FUN = "/")

meta_all <- read.csv("space_all.csv")
rownames(meta_all)  <- gsub("_TZO","",meta_all$id)
meta_skin <- meta_all[meta_all$Tissue=="Skin",]

meta_skin <- meta_skin[ rownames(meta_skin) %in% colnames(transcript),]
meta_skin <- meta_skin[rownames(meta_skin) %in% colnames(relative_common_trans),]
common2 <- relative_common_trans[,rownames(meta_skin)] /relative_common_genomic[,meta_skin$id]
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
  geom_point(color = "#4DBBD5", size = 1.5, alpha = 0.5) +  
  geom_smooth(method = "lm", color = "black", linetype = 4, se = TRUE, level = 0.9999999999999999) + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
               formula = y ~ x, 
               parse = TRUE) + 
  ggtitle("Scatter Plot with Regression Line, p-value, and R²") + 
  xlab("log10(Relative abundance)") + 
  ylab("log10(Relative contribution to the MT)") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, size = 1))  


