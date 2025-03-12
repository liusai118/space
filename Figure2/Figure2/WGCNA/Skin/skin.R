setwd("Result5/WGCNA/Skin")
library(WGCNA)
femData = read.csv("WGCNA_skin.csv",row.names = 1) 
femData<- log10(log10(femData+1) + 1)
datExpr0 = as.data.frame(t(femData)) 
gsg = goodSamplesGenes(datExpr0, verbose = 6);
gsg$allOK
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 9, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 9, minSize = 5)

table(clust)   
keepSamples = (clust==1)  
datExpr = datExpr0[keepSamples, ]  
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
traitData = read.csv("WGCNA_skin_meta.csv",row.names = 2)
traitData = traitData[,-1]
traitData = traitData[,-3]
dim(traitData)
names(traitData)
Samples = rownames(datExpr);
traitRows = match(Samples, rownames(traitData));
datTraits = traitData[traitRows, -1];
#rownames(datTraits) = traitData[traitRows, 2];


collectGarbage()



sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) 
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
type = "unsigned"
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
sft$powerEstimate
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 1, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM1",
                       verbose = 8)
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];


MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships"))

lantency = as.data.frame(datTraits$Age);
names(lantency) = "glycemia";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, lantency, use = "p"));#和体重性状的关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(lantency), sep="");
names(GSPvalue) = paste("p.GS.", names(lantency), sep="");

ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees =intramodularConnectivity(ADJ, moduleColors)
write.csv(Alldegrees, file = "intramodularConnectivity.csv")

pdf("GS vs. degree_weight3.pdf",width = 14,height = 8)

par(mfrow=c(4,5))

par(mar = c(5,5,3,3))
colorlevels=unique(moduleColors)
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]]; 
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees$kWithin[restrict1], 
                     geneTraitSignificance[restrict1,1], 
                     col=moduleColors[restrict1],
                     main=whichmodule, 
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}



dev.off()
colorlevels <- colorlevels[c(3,5,8,11)]

svg("GS vs. degree_weight3.svg",width = 3 ,height = 10)
par(mfrow=c(4,1))
par(mar = c(5,5,3,3))
for (i in c(1:length(colorlevels))) {whichmodule=colorlevels[[i]]; 
restrict1 = (moduleColors==whichmodule);
verboseScatterplot(Alldegrees$kWithin[restrict1], 
                   geneTraitSignificance[restrict1,1], 
                   col=moduleColors[restrict1],
                   main=whichmodule, 
                   xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)}
dev.off()
a1 <- unique(names(datExpr)[moduleColors %in% colorlevels])

#################
data <- t(femData[a1,])
data <- scale(data)
data <- merge(traitData,data,by="row.names")
data <- data[,-c(2,3)]
write.csv(data,"../../ML/data/skin_machine.csv")


