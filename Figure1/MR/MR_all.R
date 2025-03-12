library(ggplot2)
library(dplyr)
library(TwoSampleMR)
############
All <- rbind(res_all_gut,res_all_oral,res_all_skin)
All$id.exposure <- gsub("\\.\\.", ".", All$id.exposure)
All$type <- sub("\\..*$", "", All$id.exposure)
All$name <- sub("^[^.]+\\.([^.]+)\\..*$", "\\1", All$id.exposure)
All$id <- sub("^[^.]*\\.[^.]*\\.[^.]*\\.(\\d+).*", "\\1", All$id.exposure)
mydata_or <- generate_odds_ratios(All)
mydata_or <-subset(mydata_or,mydata_or$method=="Inverse variance weighted")
target<-subset(mydata_or,mydata_or$pval<0.05)
All$col1 <- c("#FFFFCC","#F2DFD0",rep("#F29A7F",3),rep("#FCEBE5",4),"#F2DFD0",rep("#F29A7F",4),rep("#FCEBE5",3),
              rep("#FFFCF5",20),"#FFFFCC",rep("#FCEBE5",5),rep("#FFFCF5",2) )
All <- all[,-12]
All <- all %>%
  mutate(col2 = case_when(
    p_col == "Postive effect(P<0.05)" ~ "#d55e00",
    p_col == "Negtive effect(P<0.05)" ~ "#0072b2",
    TRUE ~ NA_character_ 
  ))
###########

library(circlize)
data < read.csv("data/MR_all.csv")
sectors = data$id.exposure
circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = sectors, xlim = cbind(data$or_lci95, data$or_uci95))
circos.trackPlotRegion(ylim = c(0.9, 1.2), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  idx = which(data$id.exposure == sector.index)
  circos.rect(xleft = data$or_lci95[idx], xright = data$or_uci95[idx],
              ybottom = 0.88, ytop = 1.18, col = data$col1[idx], border = NA)
})
circos.trackPlotRegion(ylim = c(0.9, 1.2), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  idx = which(data$id.exposure == sector.index)
  circos.segments(
    x0 = 0.5, x1 = 0.5,
    y0 = data$or_lci95[idx],
    y1 = data$or_uci95[idx],
    col = data$COL[idx], lwd = 2
  )
})
