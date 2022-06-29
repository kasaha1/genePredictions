# install Kang's basic functions package from the git-hub
if ("devtools" %in% installed.packages()[, "Package"]){cat("devtools is installed")}else(install.packages("devtools"))

devtools::install_github("kasaha1/kasaBasicFunctions")
library(kasaBasicFunctions)

#------------- Packages ----

packages <- c("tidyverse", "data.table")
kasa.instPak (packages)
packages_bioc <- c("pamr","impute")
kasa.instPak_bioc (packages_bioc)

# setwd & read dataset----
setwd("C:/Users/kangs/Documents/GitHub/genePredictions/")
## Z scoring (standardization of both training and prediction samples)
STD_method <- "Robust_Z_score" # "none", "STD", "Z_Score", "Robust_Z_score"
## you need to decided delta value

if(!file.exists("OutputPAMR")){dir.create("OutputPAMR")}


trainingDataset_r <- fread("rawData/signature.txt") %>% as.data.frame()
colnames(trainingDataset_r)[1] <- "gene"

trainingClass_r <- fread("rawData/trainingClass.txt") %>% as.data.frame()
colnames(trainingClass_r) <- c("sample","class")

predictDataset_r <- fread("rawData/testDataset.txt") %>% as.data.frame()
colnames(predictDataset_r)[1] <- "gene"

# data cleaning ----

# duplication removal
trainingDataset <-trainingDataset_r %>% kasa.duplicationRemovalBySD()
colnames(trainingDataset)[1] <- "gene"
predictDataset <- predictDataset_r %>% kasa.duplicationRemovalBySD()
colnames(predictDataset)[1] <- "gene"


# sample matching

tmp <- trainingDataset %>% kasa.transposeMatrix()
tmp.1 <- kasa.matchingRow(dataframeX = tmp,dataframeY = trainingClass_r, keycolX = "sample",keycolY = "sample")
trainingDataset <- tmp.1$dataframeX %>% kasa.transposeMatrix()
colnames(trainingDataset)[1] <- "gene"
trainingClass <- tmp.1$dataframeY
colnames(trainingClass) <- c("sample","class")

## unmatched samples
matchedSamples <- trainingClass$sample %>% as.vector()
unmatchedSamples.t <- tmp$sample[!(tmp$sample %in% matchedSamples)]
unmatchedSamples.c <- trainingClass$sample[!(trainingClass$sample%in% matchedSamples)]

# gene matching
tmp <- kasa.matchingRow(dataframeX = trainingDataset,dataframeY = predictDataset,keycolX = "gene",keycolY = "gene")
trainingDataset <- tmp$dataframeX
predictDataset <- tmp$dataframeY
colnames(trainingDataset)[1] <- "gene"
colnames(predictDataset)[1] <- "gene"

## unmatched gene
matchedgene <- trainingDataset$gene %>% as.vector()
unmatchedGene.t <- trainingDataset_r$gene[!(trainingDataset_r$gene %in% matchedgene)]
unmatchedGene.p <- predictDataset_r$gene[!(predictDataset_r$gene %in% matchedgene)]

# print unmatched
write.table(unmatchedSamples.t,file = "OutputPAMR/unmatchedSamples_training.txt",quote=F,row.names = F)
write.table(unmatchedSamples.c,file = "OutputPAMR/unmatchedSamples_class.txt",quote=F,row.names = F)
write.table(unmatchedGene.t,file = "OutputPAMR/unmatchedGene_training.txt",quote=F,row.names = F)
write.table(unmatchedGene.p,file = "OutputPAMR/unmatchedGene_prediction.txt",quote=F,row.names = F)

# Standardization by STD_method ---- "none", "STD", "Z_Score", "Robust_Z_score"
switch (
  STD_method,
  STD = {
    trainingDataset <-
      trainingDataset %>% kasa.geneMedianCentering() %>% kasa.geneStandardization()
    predictDataset <-
      predictDataset %>% kasa.geneMedianCentering() %>% kasa.geneStandardization()
  },
  Z_Score = {
    trainingDataset <- trainingDataset %>% kasa.geneZscoring()
    predictDataset <- predictDataset %>% kasa.geneZscoring()
  },
  Robust_Z_score = {
    trainingDataset <- trainingDataset %>% kasa.geneRobustZscoring()
    predictDataset <- predictDataset %>% kasa.geneRobustZscoring()
  }
)


# Data wrangling ----

trainingGeneMatrix <- trainingDataset[-1]
trainingClass <- trainingClass[-1]
trainingGenenames <- trainingDataset[1]

## check missing value of matrix
# kasa.dataCleaning(trainingGeneMatrix) # clear
## missing value is imputed using nearest neighbor averaging
trainingX <- trainingGeneMatrix %>% data.matrix() %>% impute.knn()
trainingX <- trainingX$data
kasa.dataCleaning(as.data.frame(trainingX))

colnames(trainingGenenames)[1] <- "genenames"
Genenames <- trainingGenenames %>% t %>% as.vector()

## transform 'class name' to numeric vector
levOfClass <- trainingClass %>% t %>% as.factor() %>% levels()
levOfClass_numeric <- levOfClass %>% as.factor() %>% as.numeric()
table.levels <- cbind(levOfClass,levOfClass_numeric) %>% as.data.frame()
trainingClass.m <- left_join(trainingClass,table.levels, by=c("class"="levOfClass"))
trainingY <- trainingClass.m$levOfClass_numeric %>% t %>% as.character() %>% as.numeric()


## merging dataset for training
mydata <- list(x=trainingX, y=trainingY, genenames = Genenames, geneid = c(1:length(Genenames)) )

## analysis start : Training ----
model <- pamr.train(mydata)
print(model)


## analysis start : CrossValidation ----

model.cv <- pamr.cv(fit = model, data = mydata)
print(model.cv)


pdf("OutputPAMR/CrossValidation_Plot.pdf",width = 15, height = 10,pointsize = 12)
pamr.plotcv(model.cv)
dev.off()


## analysis start : Threshold 0
Delta <- 0 

pdf("OutputPAMR/CrossValidation_Plot_0.pdf",width = 15, height = 10,pointsize = 12)
pamr.plotcvprob(model, mydata,threshold = Delta)
dev.off()

## analysis start : centroid
centroid_gene <- pamr.listgenes(model, mydata, Delta, genenames = TRUE) %>% as.data.frame
centroid_gene$id <- centroid_gene$id %>% as.numeric()
centroid_gene <- centroid_gene %>% arrange(id)
colnames(centroid_gene)[c(3:(2+length(levOfClass)))] <- levOfClass
write_delim(centroid_gene,file = "OutputPAMR/centroid_ByGenes.txt",delim = "\t")

## analysis start : prediction

## testDataset preparing
testDataset <- predictDataset
colnames(testDataset)[1] <- c("genenames")
testDataset.modi <- left_join(x=trainingGenenames,y=testDataset,by=c("genenames"))


kasa.dataCleaning(testDataset.modi)

# impute missing value
testX.p <- testDataset.modi[-1] %>% data.matrix() %>% impute.knn()
testX <- testX.p$data
kasa.dataCleaning(as.data.frame(testX))


res.class <- pamr.predict(fit = model,newx = testX,threshold = Delta,type = "class") %>% as.character() %>% as.numeric()
res.class.t <- levOfClass[res.class]

res.probability <- pamr.predict(fit = model,newx = testX,threshold = Delta,type = "posterior")
res.probability.m <- res.probability %>% round(digits = 3)
res.probability.t <- apply(res.probability, 1,max) %>% round(digits = 3)

res.ID <- colnames(testX)

res.table <- cbind(res.ID,res.class.t,res.probability.t,res.probability.m) %>% as.data.frame()
colnames(res.table) <- c("Sample","Class","posterior",levOfClass)
write_delim(x = res.table,file = "OutputPAMR/PredictionResult.txt",delim = "\t")
write_delim(x = trainingDataset,file = "OutputPAMR/trainingDataset_used.txt",delim = "\t")
write_delim(x = predictDataset,file = "OutputPAMR/testDataset_used.txt",delim = "\t")

