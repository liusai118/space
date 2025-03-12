setwd("Result6/R/Oral/")
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
hebing(otutab, taxonomy2$Species, margin = 1, act = "sum") -> Species
write.csv(Species,"phylum_spe.csv")

stackplot(Species, meta_oral,
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
  geom_point(color = "#E64B35", size = 1.5, alpha = 0.5) +  
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

##########

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
