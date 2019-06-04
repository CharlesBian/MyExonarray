library(affy)
library(oligo)
library(genefilter)
library(pd.huex.1.0.st.v2)
setwd("/cluster/guest/yzhang/ExonArray/deal")
dir_cels = '/cluster/guest/yzhang/ExonArray/'
cel <- list.celfiles('/cluster/guest/yzhang/ExonArray/',full.name = TRUE)
cel_read <-read.celfiles(cel)

#for detect the assay
for(i in 1:length(sampleNames(cel_read))){
  name = paste("image",sampleNames(cel_read)[i],".jpg",sep="")
  image(cel_read[,i])
}

# normalize
fit1 <- fitProbeLevelModel(cel_read)
boxplot(fit1)
eset <- rma(cel_read)
EXP <-exprs(eset)

#filter the probe
xpa <- paCalls(cel_read)
head(xpa)
AP <-apply(xpa,1,function(x) any(x < 0.05))
xids <- as.numeric(names(AP[AP]))
head(xids)
pinfo <- getProbeInfo(cel_read)
head(pinfo)
fids <- pinfo[pinfo$fid %in% xids, 2]
head(fids)
data.exprs <- EXP[rownames(EXP) %in% fids,]

nrow(EXP)
nrow(data.exprs)
write.csv(data.exprs,"exonarray_exprs.csv")

MAplot(EXP)
write.exprs(eset)


























