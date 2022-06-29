# Functions ----
BRBarrayInstall <- function ()
{
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")
  BiocManager::install(update = FALSE, ask = FALSE)
  
  # https://brb.nci.nih.gov/BRB-ArrayTools/ArrayToolsRPackages.html
  packages <- c("ROC","classpredict")
  new.pkg <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) {
    
    BiocManager::install(pkg = "ROC", update = FALSE, ask = FALSE)
    install.packages("https://brb.nci.nih.gov/BRB-ArrayTools/RPackagesAndManuals/classpredict_0.2.tar.gz", repos=NULL)
  }
  sapply(packages, require, character.only = TRUE)
}


# install Kang's basic functions package from the git-hub
if ("devtools" %in% installed.packages()[, "Package"]){cat("devtools is installed")}else(install.packages("devtools"))

devtools::install_github("kasaha1/kasaBasicFunctions")
library(kasaBasicFunctions)

#------------- Packages ----

packages <- c("tidyverse", "data.table")
kasa.instPak (packages)
BRBarrayInstall()

# setwd & read dataset----
setwd("C:/Users/kangs/Documents/GitHub/genePredictions/")
## Z scoring (standardization of both training and prediction samples)
STD_method <- "Robust_Z_score" # "none", "STD", "Z_Score", "Robust_Z_score"

if(!file.exists("OutputBRBarray")){dir.create("OutputBRBarray")}

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
write.table(unmatchedSamples.t,file = "OutputBRBarray/unmatchedSamples_training.txt",quote=F,row.names = F)
write.table(unmatchedSamples.c,file = "OutputBRBarray/unmatchedSamples_class.txt",quote=F,row.names = F)
write.table(unmatchedGene.t,file = "OutputBRBarray/unmatchedGene_training.txt",quote=F,row.names = F)
write.table(unmatchedGene.p,file = "OutputBRBarray/unmatchedGene_prediction.txt",quote=F,row.names = F)

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


# data preparing for class prediction ----
geneId_x_ <- trainingDataset$gene %>% as.data.frame()
colnames(geneId_x_) <- c("UniqueID")
x_x_ <- trainingDataset[-1]
filter_x_ <- rep(1,nrow(trainingDataset))
expdesign_x_ <- trainingClass

exprTrain_x_ <- x_x_
exprTest_x_ <- predictDataset[-1]


k <- count(expdesign_x_,class)
class_t.l1 <- k$n[1] %>% as.numeric()
class_t.l2 <- k$n[2] %>% as.numeric()

# projectPath <-"OutputBRBarray"
projectPath <- paste0(getwd(),"/OutputBRBarray")
outputName <- "classPrediction"
generateHTML <- TRUE
prevalence <- c(class_t.l1/(class_t.l1+class_t.l2),class_t.l2/(class_t.l1+class_t.l2))
names(prevalence) <- c(k$class[1], k$class[2])
geneId_m <- geneId_x_[c("UniqueID")]
cls_t <- expdesign_x_$class %>% as.vector()
resList <- classPredict(exprTrain = exprTrain_x_, exprTest = exprTest_x_, isPaired = FALSE, 
                        pairVar.train = NULL, pairVar.test = NULL, geneId_m,
                        cls = cls_t,
                        pmethod = c("ccp", "bcc", "dlda", "knn", "nc", "svm"), 
                        geneSelect = "igenes.univAlpha",
                        univAlpha = 0.001, univMcr = 0, foldDiff = 0, rvm = TRUE, filter = filter_x_, 
                        ngenePairs = 25, nfrvm = 10, cvMethod = 1, kfoldValue = 10, bccPrior = 1, 
                        bccThresh = 0.8, nperm = 0, svmCost = 1, svmWeight =1, fixseed = 1, 
                        prevalence = prevalence, projectPath = projectPath, 
                        outputName = outputName, generateHTML = generateHTML)
if (generateHTML)
  browseURL(file.path(projectPath, "Output", outputName,
                      paste0(outputName, ".html")))

# Output
write_delim(x = trainingDataset,file = paste0(projectPath,"/1_trainingDataset.txt"),delim = "\t")
write_delim(x = predictDataset,file = paste0(projectPath,"/3_predictDataset.txt"),delim = "\t")
write_delim(x = trainingClass,file = paste0(projectPath,"/2_trainingClass.txt"),delim = "\t")
write_delim(x=resList$predNewSamples,file = paste0(projectPath,"/4_PredictionResults.txt"),delim = "\t")
write_delim(x=resList$probNew,file = paste0(projectPath,"/5_probability_BCCP.txt"),delim = "\t")