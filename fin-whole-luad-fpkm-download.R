library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(maftools) # v‘2.2.10’

#2020.09.10
#packageVersion("TCGAbiolinks") - ‘2.14.1’
#packageVersion("tidyverse") - ‘1.3.0’
#packageVersion("SummarizedExperiment") - ‘1.16.1’

#2020.12.08
packageVersion('TCGAbiolinks') #‘2.18.0’
packageVersion("tidyverse") #‘1.3.0’
packageVersion("SummarizedExperiment") #‘1.20.0’
packageVersion("maftools") #‘2.6.0’

#> getwd()
# [1] "/Users/ferrenaa/Dropbox (EinsteinMed)/data/tam/aging/tcga_updated2020.09.23"

setwd('rawdata/')

#TCGA Biolinks has three steps
# 1 - query data, find it on the GDC portal
# 2 - actually download the data
# 3 - parse the data and read into R

#1 - query data
CancerProject <- "TCGA-LUAD"


query <- GDCquery(project = CancerProject,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM",
                  sample.type = "Primary Tumor")



#2 - download data. downloads to current working directory
GDCdownload(query)


#3 - parse and read in, as a RangedSummarizedExperiment-class object
rse <- GDCprepare(query = query)





#1 - get the gene expression matrix from the object
gem <- assays(rse)$`HTSeq - FPKM`

dim(gem); gem[1:5, 1:5] #print dimensions and a small sample of the (huge) gene exp mat





#2 - row ranges: stores the ensembl ID and the corresponding HUGO gene name
geneIDs <- as.data.frame(rowRanges(rse))

head(geneIDs) #take a quick look at the gene IDs object format





#3 - metadata, including sample IDs and some clinical / demographic data:
coldata <- colData(rse)

names(coldata) #print the fields included in this dataset
dim(coldata) #print dims of dataset: 533 rows, corresponding to 533 columns (samples) in gene exp mat





### parse everything ###

#get gem & coldata
gem <- assays(rse)$`HTSeq - FPKM`

coldata <- as.data.frame(colData(rse))

#read in purity estimates and add to coldata
cpe <- as.data.frame(readxl::read_excel('purityestimates.xlsx', skip = 3))[,1:7]
cpe <- cpe[cpe[,2] == 'LUAD',]

#remove samples not in my dataset...
cpe <- cpe[cpe$`Sample ID` %in% coldata$sample,]

#remove use cols and rename cols
cpe <- cpe[,c(1,7)]
colnames(cpe) <- c('sample', 'cpe')

#make sure cpe is numeric
cpe$cpe <- as.numeric(cpe$cpe)

#remove samples in coldata duplicated at "sample" level
# these samples appear to have multiple "portions"

coldata <- coldata[!duplicated(coldata$sample),]

#match cpe row order to coldata row order based on samples
coldata[1:5,1:3]
cpe[1:3,1:2]

cpe <- cpe[match(coldata$sample, cpe$sample),]

coldata$CPE <- cpe$cpe

#remove sample(s) without CPE estimate
coldata <- coldata[!is.na(coldata$CPE),]

rm(cpe)

#select the patients...
# must have age reported... use age_at_index
# must have kras mut (non-silent)
# if duplicate samples exist: check why...
# if due to multiple "vials" --> pick latest vial...

#get patients with age; use age_at_index; 
#some of these have age at diagnosis missing
#if look at non-NAs for age at diagnosis, there are no age at index missings...
cdx <- coldata[!is.na(coldata$age_at_index),]

#remove any dups: seem to be caused by multiple vials...
dup_patients <- cdx[duplicated(cdx$patient),"patient"]

message('\n\nMake sure the duplicates are due to multiple vials!!!\n\n')

#for each dup sample, get the highest vial number
samps_to_keep <- c()
for(dup_pat in dup_patients){
  
  #subset coldata for each duplicated patient
  cdx_dups <- cdx[cdx$patient == dup_pat,]
  
  #get vials; remap as numeric, select "highest" vial...
  # this paper selects B vials instead of A vials: https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30201-7
  vials <- substr(cdx_dups$sample, start = 16, stop = 16)
  
  vials <- plyr::mapvalues(factor(vials, levels = c('A','B','C')), 
                           from = c('A','B','C'), 
                           to = c(1,2,3))
  
  samp_retain <- cdx_dups$sample[which.max(vials)]
  
  samps_to_keep <- c(samps_to_keep, samp_retain )
  
  rm(cdx_dups, vials, samp_retain)
}

#samples to keep: add the high vial number dup patients
# plus the non-dup sample IDs

samps_to_keep <- c(samps_to_keep, cdx[!(cdx$patient %in% dup_patients), "sample"])


#final subset of the coldata and the GEM
coldata <- coldata[coldata$sample %in% samps_to_keep,]
gem <- gem[,colnames(gem) %in% rownames(coldata)]

rm(cdx, samps_to_keep)


table(colnames(gem) == rownames(coldata))




#final check 
dim(gem)


#try to convert the gem gene names...
geneIDs <- as.data.frame(rowRanges(rse))
geneIDs <- geneIDs[!duplicated(geneIDs$external_gene_name),]

gem <- gem[rownames(gem) %in% geneIDs$ensembl_gene_id,]
rownames(gem) <- geneIDs$external_gene_name



saveRDS(list(gem, coldata, geneIDs), 'luad-whole_gem_cd_list_FPKM.rds')
#saveRDS(list(gem, coldata, geneIDs), 'kras_gem_cd_list_FPKM.rds')
#saveRDS(list(gem, coldata, geneIDs), 'egfr_gem_cd_list_FPKM.rds')



