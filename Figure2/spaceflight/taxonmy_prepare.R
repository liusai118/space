taxonomy1 <- read.csv("taxon.csv")
taxonomy2 <-taxonomy1[,8:2]
rownames(taxonomy2) <- taxonomy2$species
colnames(taxonomy2) <- c( "Kingdom" ,"Phylum"  ,"Class" ,  "Order" ,  "Family",  "Genus" ,   "Species")
library(microeco)
taxonomy3 <- taxonomy2[!taxonomy2$Kingdom=="0",]
taxonomy3 <- taxonomy3[!taxonomy3$Phylum=="0",]
taxonomy3 <- taxonomy3[!taxonomy3$Family=="0",]
taxonomy3 <- taxonomy3[!taxonomy3$Genus=="environmental samples",]
taxonomy3 <- taxonomy3[!taxonomy3$Species=="0",]
taxonomy3 <- taxonomy3[!taxonomy3$Order=="0",]
taxonomy3 <- taxonomy3[!taxonomy3$Class=="0",]
taxonomy3 <- taxonomy3[!taxonomy3$Genus=="0",]
taxonomy3 <- tidy_taxonomy(taxonomy3)
taxonomy3 <- taxonomy3[!taxonomy3$Genus=="g__",]
taxonomy3<- taxonomy3[rownames(taxonomy3) %in% rownames(otutab),]
