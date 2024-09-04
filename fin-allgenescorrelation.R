library(tidyverse)
library(biomaRt)
library(foreach)
library(doParallel)
library(TCGAbiolinks)
library(olsrr)


set.seed(2020)

#setwd('~/Documents/tam/aging/tcga_updated2020.09.23/') --> everything is within tcga_updated2020.09.23/, 
# keeping all navigation relative / local to here


### set parameters for run ###


#dir to write to
outdir <- "finalized_start-to-end-2021.02.20/all-genes-correlation/"
dir.create(outdir)







### define functions ###

#percent sign operators...
# https://stackoverflow.com/questions/12730629/what-do-the-op-operators-in-mean-for-example-in


scale=T
center=T


regressout <- function(gem, var_to_reg, scale=NULL, center=NULL){
  if(is.null(scale)){scale=T}
  if(is.null(center)){center=T}
  
  genes <- rownames(gem)
  
  if(!(scale == F & center == F )){
    
    gem <- as.data.frame(t(scale(t(gem), center=center, scale = scale)))
  }
  
  message('Initiating regression for matrix')
  total = nrow(gem) - 1
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  l <- vector(mode = 'list', length = nrow(gem))
  for(i in 1:nrow(gem)){
    x <- t(gem[i,])[,1]
    l[[i]] <- lm(x ~ var_to_reg)$residuals
    
    setTxtProgressBar(pb, i)
  }
  rgem <- as.data.frame(bind_rows(l))
  rownames(rgem) <- genes
  
  rgem
}












#read in data
gemcdl <- readRDS('rawdata/kras_gem_cd_list_FPKM.rds')
gem <- gemcdl[[1]]
md <- gemcdl[[2]]
geneids <- gemcdl[[3]]
rm(gemcdl)



# if(removeoutlier == T){
#   
#   gem <- gem[,colnames(gem) != 'TCGA-44-3917-01B-02R-A277-07']
#   md <- md[rownames(md) != 'TCGA-44-3917-01B-02R-A277-07',]
# }







#fpkm to TPM FPKMij / sum(FPKM,j)) * 1M
gemt <- apply(gem, 2, function(x){(x / sum(x))*10^6})

#quantile norm
gemt <- as.data.frame(preprocessCore::normalize.quantiles(gemt))
rownames(gemt) <- rownames(gem)
colnames(gemt) <- colnames(gem)



gem <- gemt
rm(gemt)






# filter for protein coding genes only
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
list <- getBM(attributes=c("external_gene_name","transcript_biotype", 'ensembl_gene_id'),
               filters = "ensembl_gene_id", values=geneids$ensembl_gene_id, mart=human)


coding <- list[list$transcript_biotype == 'protein_coding',]


# saveRDS(coding, 'rawdata/codinggenes.rds')
# coding <- readRDS('rawdata/codinggenes.rds')


#update 2021.10.26 --> some genes (like H2AFJ) are inconsisten in their symbol.
# need to use ensemblID primarily.

geneids <- geneids[match(coding$ensembl_gene_id, geneids$ensembl_gene_id),]

#subset, or not
gem <- gem[rownames(gem) %in% geneids$external_gene_name,]

rm(coding,list, geneids, human)








###   load up ####

#remove empty gene names (seems not to be a thing)
gem <- gem[rownames(gem) != '',]


#reduce gem size
gem <- gem[rowSums(gem) >= 10,]


#remove low CPE
md <- md[md$CPE >= 0.5,]
gem <- gem[,colnames(gem) %in% rownames(md)]

#remove missing CPE
md <- md[!is.na(md$CPE),]
gem <- gem[,colnames(gem) %in% rownames(md)]

#regress
gem <- regressout(gem = gem, var_to_reg = md$CPE, scale = T, center = T)








#get the age out
age <- md$age_at_index

#progress bar
total = nrow(gem) - 1
pb <- txtProgressBar(min = 0, max = total, style = 3)

reslist <- list()
for(genedex in c(1:nrow(gem)) ){
  
  #get gene name
  gene <- rownames(gem)[genedex]
  
  #get expression of the gene
  genevec <- as.vector(gem[genedex,])
  
  #make a gene df
  pdf <- data.frame('Age' = age,
                    'Gene' = t(genevec))
  
  #calculate cor
  cres <- cor.test(pdf[,1], pdf[,2])
  
  #save as a df
  resdf <- data.frame(genedex = genedex,
                      gene = gene,
                      cor = cres$estimate,
                      p = cres$p.value,
                      row.names = NULL, stringsAsFactors = F
  )
  
  #save df to list
  reslist[[genedex]] <- resdf
  
  
  setTxtProgressBar(pb, genedex)
  
}


#bind list as df and order by cor
res <- dplyr::bind_rows(reslist)

res <- res[order(res$cor, decreasing = T),]

#get FDR
res$fdr <- p.adjust(res$p, method = 'fdr')

beepr::beep()

res <- res[,2:5]


res_file <- paste0(outdir, '/all-genes-cor-resutls.csv')

write.csv(res, file = res_file,
          row.names = F)







#test top and bottom gene, select gene, plot, and run MVR

gene <- 'KLF4'
gene <- 'H2AFJ'

smoke = F


#set up, get age, gene exp vector, and cor/FDR
genevec <- t(gem[rownames(gem)==gene,])[,1]
generes <- res[res$gene == gene,]




#run without smoke

if(smoke == F){
  
  pdf <- data.frame(age= md$age_at_index, 
                    genevec = genevec)
  
  nicep <- ifelse(test = generes$fdr < 0.01, formatC(generes$fdr, format='e', digits=2), 
                  no = round(generes$fdr, 2))
  
  
  g <- ggplot(pdf, aes(age, genevec))+
    geom_point(size=3)+
    geom_smooth(method = 'lm', formula = y~x, se = T)+
    labs(title = gene,
         subtitle = paste0('N = ', ncol(gem),
                           '\nPearson cor r = ', round(generes$cor, 3),
                           '\nFDR (cor.test) = ', nicep))+
    ylab('Z-scaled Quantile-Normalized TPM')+ xlab('Age')+
    theme_linedraw() 
  
}


#jackknife plot with smoke status.
if(smoke==T){
  pdf <- data.frame(age= md$age_at_index, 
                    genevec = genevec,
                    pack_years_smoked = md$pack_years_smoked)
  

  
  pdf$smoking_status <- md$paper_Smoking.Status
  pdf[is.na(pdf$smoking_status), "smoking_status"] <- '[Not Available]'
  
  #levels
  # [1] "Current reformed smoker for > 15 years"     
  # [2] "Current reformed smoker for < or = 15 years"
  # [3] "Current smoker"                             
  # [4] "Lifelong Non-smoker"                        
  # [5] "[Not Available]" 
  
  pdf$smoking_status <- as.vector(pdf$smoking_status)
  
  pdf$smoker_status_cat <- pdf$smoking_status
  pdf$smoker_status_cat[pdf$smoker_status_cat %in% c('Current reformed smoker for > 15 years', 'Current reformed smoker for < or = 15 years', 'Current smoker')] <- 'History of smoking'
  
  #if categorizied is missing, but pack-years are recorded, then consider them as smoke history
  pdf[pdf$smoker_status_cat == '[Not Available]' & !is.na(pdf$pack_years_smoked),"smoker_status_cat"] <- 'History of smoking'
  
  nicep <- ifelse(test = generes$fdr < 0.01, formatC(generes$fdr, format='e', digits=2), 
                  no = round(generes$fdr, 2))
  
  pdf$smoker_status_cat <- factor(pdf$smoker_status_cat, levels = c('History of smoking', 'Lifelong Non-smoker', '[Not Available]'))
  
  #cor of smoke with module score
  cres_smoke <- cor.test(pdf$genevec, pdf$pack_years_smoked)
  nicep_smoke <- ifelse(test = cres_smoke$p.value < 0.01, formatC(cres_smoke$p.value, format='e', digits=2), 
                        no = round(cres_smoke$p.value, 2))
  
  
  
  g <- ggplot(pdf, aes(x = age, y = genevec, col = pack_years_smoked )) +
    geom_point(size = 3, aes(shape = smoker_status_cat))+
    scale_shape_manual(values = c(16,17,3), 'History of Smoking')+
    geom_smooth(method = 'lm', formula = y~x, se = T)+
    labs(title = gene,
         subtitle =   paste0('N = ', ncol(gem),
                             '\nPearson cor r = ', round(generes$cor, 3), ' ; Smoke cor with gene = ',round(cres_smoke$estimate,3),
                             '\nFDR (cor.test) = ', nicep, ' ; Smoke P = ',nicep_smoke))+
    ylab('Z-scaled Quantile-Normalized TPM')+ xlab('Age')+
    scale_color_gradient(low = '#B4B4B4', high = '#2D2D2D', na.value = 'purple3', 'Pack-Years Smoked\nPurple=NA')+
    theme_linedraw()
    

}


print(g)


#site of resection
keep <- c('Upper lobe, lung', 'Lower lobe, lung')
md$site_of_resection_or_biopsy_recoded <- md$site_of_resection_or_biopsy
md[!(md$site_of_resection_or_biopsy_recoded %in% keep) ,"site_of_resection_or_biopsy_recoded"] <- 'Other lung'





#recode met and AJCC
md$ajcc_pathologic_t_recoded <- as.vector(md$ajcc_pathologic_t)
md$ajcc_pathologic_t_recoded[md$ajcc_pathologic_t_recoded %in% unique(md$ajcc_pathologic_t_recoded)[grepl('T1', unique(md$ajcc_pathologic_t_recoded))]] <- "T1"
md$ajcc_pathologic_t_recoded[md$ajcc_pathologic_t_recoded %in% unique(md$ajcc_pathologic_t_recoded)[grepl('T2', unique(md$ajcc_pathologic_t_recoded))]] <- "T2"

md$ajcc_pathologic_t_recoded <- factor(md$ajcc_pathologic_t_recoded)

md$ajcc_pathologic_m_recoded <- as.vector(md$ajcc_pathologic_m)
md$ajcc_pathologic_m_recoded[md$ajcc_pathologic_m_recoded %in% unique(md$ajcc_pathologic_m_recoded)[grepl('M1', unique(md$ajcc_pathologic_m_recoded))]] <- "M1"



#try to grab overall mut burdern...
mut <- readRDS('rawdata/mut_maf.rds')

mut <- as.data.frame(maftools::getSampleSummary(mut))

#need to deal with different codes...
mut$Tumor_Sample_Barcode <- substr(mut$Tumor_Sample_Barcode, start = 1, stop = 12)

#kkeep only the ones in the MD (kras cohort)
mut <- mut[mut$Tumor_Sample_Barcode %in% substr(rownames(md),start = 1, stop = 12) ,]

#make sure they're sorted properly
mut <- mut[match(substr(rownames(md),start = 1, stop = 12), mut$Tumor_Sample_Barcode),]




# this is without standardized beta

pdf <- data.frame(Age = md$age_at_index,
                  genevec = genevec,
                  cpe = md$CPE,
                  pack_years_smoked = md$pack_years_smoked,
                  mutation_burden = mut$total,
                  
                  ajcc_pathologic_t = md$ajcc_pathologic_t_recoded ,
                  ajcc_pathologic_n = forcats::fct_infreq(md$ajcc_pathologic_n) ,
                  ajcc_pathologic_m = forcats::fct_infreq(md$ajcc_pathologic_m_recoded) ,
                  site_of_resection_or_biopsy_recoded = forcats::fct_infreq(md$site_of_resection_or_biopsy_recoded) ,
                  gender = forcats::fct_infreq(md$gender) 
                  
                  
)








#tab 2 bivar assoc
cor.test(pdf$genevec, pdf$Age)



#repeat it but try to get standardized betas
pdf_std <- data.frame(Age = scale(md$age_at_index),
                      genevec = genevec,
                      cpe = scale(md$CPE),
                      pack_years_smoked = scale(md$pack_years_smoked),
                      mutation_burden = scale(mut$total),
                      
                      ajcc_pathologic_t = md$ajcc_pathologic_t_recoded ,
                      ajcc_pathologic_n = forcats::fct_infreq(md$ajcc_pathologic_n) ,
                      ajcc_pathologic_m = forcats::fct_infreq(md$ajcc_pathologic_m_recoded) ,
                      site_of_resection_or_biopsy_recoded = forcats::fct_infreq(md$site_of_resection_or_biopsy_recoded) ,
                      gender = forcats::fct_infreq(md$gender) 
                      
                      
)




#omit NAs (mostly smoke...)
pdf <- na.omit(pdf)
pdf_std <- na.omit(pdf_std)


library(olsrr)









#full model
m <- lm(data = pdf, genevec ~ .)
summary(m)



#diagnostics
df <- data.frame(resid=m$residuals,
                 yhat = m$fitted.values)

#RVF plot
ggplot(df, aes(yhat, resid))+
  geom_point(size=2, shape=1, col='blue')+
  geom_hline(yintercept = 0, col='red', linetype='dashed')

#heteroscedasticity test
ols_test_breusch_pagan(m)

#normality of residuals
ggplot(df, aes(resid))+
  geom_histogram(aes(y=..density..), col='black', fill='steelblue')+
  stat_function(fun = dnorm, args = list(mean = mean(df$resid), sd = sd(df$resid)),
                col = 'red')

shapiro.test(df$resid)


#collinearity
ols_coll_diag(m)


#get standardized coefficients
m_std <- lm(data = pdf_std, genevec ~ .)
summary(m_std)


#prep to save; convert model results to df and add std be
m_save <- broom::tidy(m)
m_std <- broom::tidy(m_std)

m_save$estimate_standardized <- m_std$estimate


saveforapp=F
if(saveforapp == T){
  
  appdir <- 'shiny/shiny-tcga-kras/data/'
    
  #prep for app...
  # gem
  # md
  # res
  
  
  saveRDS(gem, paste0(appdir, '/gem.rds'))
  saveRDS(md, paste0(appdir, '/md.rds'))
  saveRDS(res, paste0(appdir, '/res.rds'))
  
}