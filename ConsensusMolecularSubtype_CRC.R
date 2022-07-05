# install Kang's basic functions package from the git-hub
if ("devtools" %in% installed.packages()[, "Package"]){cat("devtools is installed")}else(install.packages("devtools"))

devtools::install_github("kasaha1/kasaBasicFunctions")
library(kasaBasicFunctions)
devtools::install_github("Sage-Bionetworks/CMSclassifier")
library(CMSclassifier)

# install Kang's basic functions package from the git-hub
#------------- Packages ----
packages <- c("tidyverse", "data.table","HGNChelper")
kasa.instPak (packages)
#-----------------------------
#------------- Bioc_Packages ----
packages_bioc <- c("aliases2entrez","limma","impute")
kasa.instPak_bioc (packages_bioc)
#-----------------------------
HGNC <- update_symbols()
GeneName <- fread("rawData/CMS_reference_transcriptome.txt") %>% as.data.frame()
GeneName$GeneID <- GeneName$GeneID %>% as.character()
GeneName.list <- GeneName$GeneID %>% t() %>% as.vector()

inputdata <- fread("rawData/testDataset.txt") %>% as.data.frame()
colnames(inputdata)[1] <- "genes"
GeneList.rnaseq <- inputdata$genes %>% t() %>% as.vector()
ids <- convert_symbols(GeneList.rnaseq, HGNC,c=8) # Converts Hugo to enterzID 
inputdata_t <- inputdata
inputdata_t$genes <- ids$entrezID %>% as.character()

raw.data <- inputdata_t %>% filter(genes %in% GeneName.list) %>% kasa.duplicationRemovalBySD()
match <- inner_join(x = GeneName,y = raw.data, by=c("GeneID"="genes"))
datasetForCMS <- match[,-2] %>% kasa.transposeMatrix()
rownames(datasetForCMS) <- datasetForCMS[,1]
data.CMS <- datasetForCMS[,-1]

Rfcms <- CMSclassifier::classifyCMS(t(data.CMS),method="RF")[[3]]
SScms <- CMSclassifier::classifyCMS(t(data.CMS),method="SSP")[[3]]


result_SSCMS <- SScms %>% rownames_to_column()
result_RfCMS <- Rfcms %>% rownames_to_column()

dirOut <- "Output_CMS"
if(!file.exists(dirOut)){dir.create(dirOut)}
write_delim(result_SSCMS,file = paste0(dirOut,"/result_SSCMS.txt"),delim = "\t")
write_delim(result_RfCMS,file = paste0(dirOut,"/result_RfCMS.txt"),delim = "\t")
