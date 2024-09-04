library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(biomaRt)
library(foreach)
library(doParallel)
library(TCGAbiolinks)
library(olsrr)
library(survival)
library(survminer)


set.seed(2020)

# set outdir:
# kras cohort first


outdir <- paste0('finalized_start-to-end-2021.02.20/E-clinicalvarscor/')

dir.create(outdir)

#kras
outdir <- paste0(outdir, '/wholeluad/')
dir.create(outdir)








#read in data
gemcdl <- readRDS('rawdata/luad-whole_gem_cd_list_FPKM.rds')
gem <- gemcdl[[1]]
md <- gemcdl[[2]]
rm(gemcdl)









# prep some of the confounding variables #

#site of resection
keep <- c('Upper lobe, lung', 'Lower lobe, lung')
md$site_of_resection_or_biopsy_recoded <- md$site_of_resection_or_biopsy
md[!(md$site_of_resection_or_biopsy_recoded %in% keep) ,"site_of_resection_or_biopsy_recoded"] <- 'Other lung'





#recode met and AJCC
md$ajcc_pathologic_t_recoded <- as.vector(md$ajcc_pathologic_t)
md$ajcc_pathologic_t_recoded[md$ajcc_pathologic_t_recoded %in% unique(md$ajcc_pathologic_t_recoded)[grepl('T1', unique(md$ajcc_pathologic_t_recoded))]] <- "T1"
md$ajcc_pathologic_t_recoded[md$ajcc_pathologic_t_recoded %in% unique(md$ajcc_pathologic_t_recoded)[grepl('T2', unique(md$ajcc_pathologic_t_recoded))]] <- "T2"

#remove Tx. only 3.
md[md$ajcc_pathologic_t_recoded=='TX', 'ajcc_pathologic_t_recoded'] <- NA

md$ajcc_pathologic_t_recoded <- factor(md$ajcc_pathologic_t_recoded)

md$ajcc_pathologic_m_recoded <- as.vector(md$ajcc_pathologic_m)
md$ajcc_pathologic_m_recoded[md$ajcc_pathologic_m_recoded %in% unique(md$ajcc_pathologic_m_recoded)[grepl('M1', unique(md$ajcc_pathologic_m_recoded))]] <- "M1"


#aajcc N, remove Nx (only 9), and also remove N3 (only 3)
# there is actually 1 NA in AJCC N score too... set it to NX first...?
md$ajcc_pathologic_n_recoded <- md$ajcc_pathologic_n
md[grepl('NX', md$ajcc_pathologic_n_recoded), 'ajcc_pathologic_n_recoded'] <- NA
md[grepl('N3', md$ajcc_pathologic_n_recoded), 'ajcc_pathologic_n_recoded'] <- NA


#try to grab overall mut burdern...
mut <- readRDS('rawdata/mut_maf.rds')

mut <- as.data.frame(maftools::getSampleSummary(mut))

#need to deal with different codes...
mut$Tumor_Sample_Barcode <- substr(mut$Tumor_Sample_Barcode, start = 1, stop = 12)

#kkeep only the ones in the MD (kras cohort)
mut <- mut[mut$Tumor_Sample_Barcode %in% substr(rownames(md),start = 1, stop = 12) ,]

#make sure they're sorted properly
mut <- mut[match(substr(rownames(md),start = 1, stop = 12), mut$Tumor_Sample_Barcode),]




#set up main clinical data
clindat <- data.frame(barcode = md$barcode,
                      Age = md$age_at_index,
                      # race = md$race #almost all white
                      purity = md$CPE,
                      pack_years_smoked = md$pack_years_smoked,
                      mutation_burden = mut$total,
                      
                      ajcc_pathologic_t = md$ajcc_pathologic_t_recoded ,
                      ajcc_pathologic_n = forcats::fct_infreq(md$ajcc_pathologic_n_recoded) ,
                      ajcc_pathologic_m = forcats::fct_infreq(md$ajcc_pathologic_m_recoded) ,
                      site_of_resection_or_biopsy_recoded = forcats::fct_infreq(md$site_of_resection_or_biopsy_recoded) ,
                      gender = forcats::fct_infreq(md$gender) 
                      
                      
)


rm(mut) #save space






### mutation burden ###

corres <- cor.test(clindat$Age, clindat$mutation_burden)


nicep <- ifelse(test = corres$p.value < 0.01, formatC(corres$p.value , format='e', digits=2), 
                no = round(corres$p.value , 2))

cp <- ggplot(clindat, aes(Age, mutation_burden))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = 'WHOLE LUAD',
       subtitle = paste0('N = ', nrow(clindat),
                         '\nPearson cor r = ', round(corres$estimate, 3),
                         '\nP (cor.test) = ', nicep))+
  theme_light()


cp


pdf(paste0(outdir, '/age_mutburden_cor.pdf'), height = 5, width = 5)
cp
dev.off()


