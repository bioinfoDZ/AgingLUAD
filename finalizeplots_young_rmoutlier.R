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
outdir <- "finalized_start-to-end-2021.02.20/B_young_kras_correlation//update_final_rmoutlier"
dir.create(outdir)


#select deg list; aged is [[1]], young is [[2]]
deglist <- readRDS('finalized_start-to-end-2021.02.20/parsedmousedeglist.rds')
deg <- deglist[[2]]


### read in gene results ###
generes <- 'finalized_start-to-end-2021.02.20/B_young_kras_correlation///jackknifesummary.csv'
generes <- read.csv(generes)





#get input gene list; the object "rawgenelist" can be swapped out for any genelist.
rawgenelist <- deg$HGNC.symbol
genelist <- rawgenelist



#parallelization
# with 6, 200 genes takes few second, 1400 genes takes ~half hour...
numclust = 6




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


# genelist <- optimal
# numbins = 24
# numcontrolgenesperbin = 100

modulescore <- function(gem, genelist, numbins=NULL, numcontrolgenesperbin=NULL){
  if(is.null(numbins)){numbins=24}
  if(is.null(numcontrolgenesperbin)){numcontrolgenesperbin=100}
  
  testgenes <- genelist
  controlgenes <- rownames(gem)
  controlgenes <- controlgenes[!(controlgenes %in% testgenes)]
  
  #testmeans <- Matrix::rowMeans(gem[rownames(gem) %in% genes,])
  
  ### bin the genes by avg expression ###
  
  #get average gene expression
  avgs <- Matrix::rowMeans(gem)
  avgs <- avgs[order(avgs)]
  
  #cut; get the gene bins
  bins <- cut_number(avgs, n=numbins)
  
  bindf <- data.frame(genenames = names(avgs), 
                      bins = bins)
  
  #select control genes from same expression bins
  controlgenes <- c()
  for(genedex in 1:length(testgenes)){
    
    #get gene and bin
    gene <- testgenes[genedex]
    genebin <- bindf[bindf$genenames == gene,2]
    
    #select all genes from that bin, to sample from
    tmpcontrols <- bindf[bindf$bins == genebin,]
    
    #exclude the actual test genes
    tmpcontrols <- tmpcontrols[!(tmpcontrols$genenames %in% testgenes),]
    
    #if num controls exceeds bin size, take all genes in bin...
    numtotake <- ifelse(nrow(tmpcontrols) < numcontrolgenesperbin,
                        yes = nrow(tmpcontrols),
                        no = numcontrolgenesperbin)
    
    #sample the genes
    tmpcontrols <- sample(tmpcontrols$genenames, size = numtotake, replace = F)
    
    controlgenes <- unique(c(controlgenes, 
                             tmpcontrols
    ))
  }
  
  
  
  #get control gene mean expression for each sample
  samplemeans_controls <- Matrix::colMeans( gem[rownames(gem) %in% controlgenes,] )
  
  #get test genes mean expression for each samoke
  samplemeans_test <- Matrix::colMeans( gem[rownames(gem) %in% testgenes,] )
  
  #subtract them to get module score
  modulescore <- samplemeans_test - samplemeans_controls
  
  return(modulescore)
}














#read in data
gemcdl <- readRDS('rawdata/kras_gem_cd_list_FPKM.rds')
gem <- gemcdl[[1]]
md <- gemcdl[[2]]
rm(gemcdl)





gem <- gem[,colnames(gem) != 'TCGA-44-3917-01B-02R-A277-07']
md <- md[rownames(md) != 'TCGA-44-3917-01B-02R-A277-07',]








#fpkm to TPM FPKMij / sum(FPKM,j)) * 1M
gemt <- apply(gem, 2, function(x){(x / sum(x))*10^6})

#quantile norm
gemt <- as.data.frame(preprocessCore::normalize.quantiles(gemt))
rownames(gemt) <- rownames(gem)
colnames(gemt) <- colnames(gem)



gem <- gemt
rm(gemt)






#filter for protein coding genes only
# library(biomaRt)
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# list <- getBM(attributes=c("hgnc_symbol","transcript_biotype"), 
#                filters = "hgnc_symbol", values=rownames(gem), mart=human)
# 
# 
# coding <- list[list$transcript_biotype == 'protein_coding',]
#
# saveRDS(coding, 'rawdata/codinggenes.rds')
coding <- readRDS('rawdata/codinggenes.rds')

#subset, or not
gem <- gem[rownames(gem) %in% coding$hgnc_symbol,]

rm(coding)








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









#jackknifed module score corplot with raw gene list
finalgenelist <- generes[, "HGNC.symbol"]

m_opt <- scale(modulescore(gem, finalgenelist))

pdf <- data.frame(Age = md$age_at_index,
                  ModuleScore = m_opt)

cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))

g_jk <- ggplot(pdf, aes(x = Age, y = ModuleScore))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Raw gene list Module Score",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         ' ; Num final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3),
                         '\nP (cor.test) = ', nicep))+
  ylab('Z-Scaled Module Score')+
  theme_light() 



g_jk







#jackknife plot with smoke status.
pdf <- data.frame(ModuleScore = m_opt,
                  Age = md$age_at_index,
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

cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))

pdf$smoker_status_cat <- factor(pdf$smoker_status_cat, levels = c('History of smoking', 'Lifelong Non-smoker', '[Not Available]'))



pdf[pdf$smoker_status_cat=='Lifelong Non-smoker',"pack_years_smoked"] <- 0


#cor of smoke with module score
cres_smoke <- cor.test(pdf$ModuleScore, pdf$pack_years_smoked)
nicep_smoke <- ifelse(test = cres_smoke$p.value < 0.01, formatC(cres_smoke$p.value, format='e', digits=2), 
                      no = round(cres_smoke$p.value, 2))



g_jk_smoke <- ggplot(pdf, aes(x = Age, y = ModuleScore, col = pack_years_smoked )) +
  geom_point(size = 3)+
  scale_shape_manual(values = c(16,17,3), 'History of Smoking')+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Raw gene list Module Score + Smoking Status",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         ' ; Num final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3), ' ; Smoke cor with modulescore = ',round(cres_smoke$estimate,3),
                         '\nP (cor.test) = ', nicep, ' ; Smoke P = ',nicep_smoke))+
  ylab('Z-Scaled Module Score')+
  #scale_color_distiller(palette = 'Greys', 'Pack-years', direction = 1, na.value = 'Purple')+
  scale_color_gradient(low = '#B4B4B4', high = '#2D2D2D', na.value = 'purple3', 'Pack-Years Smoked\nPurple=NA')+
  #scale_color_gradientn(colors = rev(palfun(nrow(pdf))), limits = range(pdf$pack_years_smoked), 'Z-scored\nPack-Years')+
  #continuous_scale('color', scales::gradient_n_pal(pal), palette = 1)+
  #scale_color_viridis('Z-scaled\nPack-years', option = 'C', na.value = 'purple')+
  theme_light()


g_jk <- g_jk + theme_light() + scale_y_continuous(breaks = seq(-3,5))
g_jk_smoke <- g_jk_smoke + theme_light()+ scale_y_continuous(breaks = seq(-3,5))
#smoking plot: all lifelong nn-smkers have no pack-years.

pdffile <- paste0(outdir, '/plots_updated_noJK_rawgenelist.pdf')
pdf(pdffile, height=7, width=15)
g_jk + g_jk_smoke
dev.off()






#jackknifed module score corplot (no thresholding)
finalgenelist <- generes[generes$jackknifeclassification=='Signal', "HGNC.symbol"]

m_opt <- scale(modulescore(gem, finalgenelist))

pdf <- data.frame(Age = md$age_at_index,
                  ModuleScore = m_opt)

cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))

g_jk <- ggplot(pdf, aes(x = Age, y = ModuleScore))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Jackknifed Module Score",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         ' ; Num final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3),
                         '\nP (cor.test) = ', nicep))+
  ylab('Z-Scaled Module Score')+
  theme_light() 



g_jk







#jackknife plot with smoke status.
pdf <- data.frame(ModuleScore = m_opt,
                  Age = md$age_at_index,
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

cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))

pdf$smoker_status_cat <- factor(pdf$smoker_status_cat, levels = c('History of smoking', 'Lifelong Non-smoker', '[Not Available]'))



pdf[pdf$smoker_status_cat=='Lifelong Non-smoker',"pack_years_smoked"] <- 0


#cor of smoke with module score
cres_smoke <- cor.test(pdf$ModuleScore, pdf$pack_years_smoked)
nicep_smoke <- ifelse(test = cres_smoke$p.value < 0.01, formatC(cres_smoke$p.value, format='e', digits=2), 
                      no = round(cres_smoke$p.value, 2))



g_jk_smoke <- ggplot(pdf, aes(x = Age, y = ModuleScore, col = pack_years_smoked )) +
  geom_point(size = 3)+
  scale_shape_manual(values = c(16,17,3), 'History of Smoking')+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Jackknifed Module Score + Smoking Status",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         ' ; Num final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3), ' ; Smoke cor with modulescore = ',round(cres_smoke$estimate,3),
                         '\nP (cor.test) = ', nicep, ' ; Smoke P = ',nicep_smoke))+
  ylab('Z-Scaled Module Score')+
  #scale_color_distiller(palette = 'Greys', 'Pack-years', direction = 1, na.value = 'Purple')+
  scale_color_gradient(low = '#B4B4B4', high = '#2D2D2D', na.value = 'purple3', 'Pack-Years Smoked\nPurple=NA')+
  #scale_color_gradientn(colors = rev(palfun(nrow(pdf))), limits = range(pdf$pack_years_smoked), 'Z-scored\nPack-Years')+
  #continuous_scale('color', scales::gradient_n_pal(pal), palette = 1)+
  #scale_color_viridis('Z-scaled\nPack-years', option = 'C', na.value = 'purple')+
  theme_light()


g_jk <- g_jk + theme_light() + scale_y_continuous(breaks = seq(-3,5))
g_jk_smoke <- g_jk_smoke + theme_light()+ scale_y_continuous(breaks = seq(-3,5))
#smoking plot: all lifelong nn-smkers have no pack-years.

pdffile <- paste0(outdir, '/plots_updated_jackknifed.pdf')
pdf(pdffile, height=7, width=15)
g_jk + g_jk_smoke
dev.off()












#jackknifed module score WITH THRESHOLDING
finalgenelist <- generes[generes$thresholdedclassification=='Signal', "HGNC.symbol"]

m_opt <- scale(modulescore(gem, finalgenelist))

pdf <- data.frame(Age = md$age_at_index,
                  ModuleScore = m_opt)

cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))

g_jk <- ggplot(pdf, aes(x = Age, y = ModuleScore))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "THRESHOLDED Module Score",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         ' ; Num final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3),
                         '\nP (cor.test) = ', nicep))+
  ylab('Z-Scaled Module Score')+
  theme_light() 



g_jk







#jackknife plot with smoke status.
pdf <- data.frame(ModuleScore = m_opt,
                  Age = md$age_at_index,
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

cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))

pdf$smoker_status_cat <- factor(pdf$smoker_status_cat, levels = c('History of smoking', 'Lifelong Non-smoker', '[Not Available]'))



pdf[pdf$smoker_status_cat=='Lifelong Non-smoker',"pack_years_smoked"] <- 0


#cor of smoke with module score
cres_smoke <- cor.test(pdf$ModuleScore, pdf$pack_years_smoked)
nicep_smoke <- ifelse(test = cres_smoke$p.value < 0.01, formatC(cres_smoke$p.value, format='e', digits=2), 
                      no = round(cres_smoke$p.value, 2))



g_jk_smoke <- ggplot(pdf, aes(x = Age, y = ModuleScore, col = pack_years_smoked )) +
  geom_point(size = 3)+
  scale_shape_manual(values = c(16,17,3), 'History of Smoking')+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "THRESHOLDED Module Score + Smoking Status",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         ' ; Num final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3), ' ; Smoke cor with modulescore = ',round(cres_smoke$estimate,3),
                         '\nP (cor.test) = ', nicep, ' ; Smoke P = ',nicep_smoke))+
  ylab('Z-Scaled Module Score')+
  #scale_color_distiller(palette = 'Greys', 'Pack-years', direction = 1, na.value = 'Purple')+
  scale_color_gradient(low = '#B4B4B4', high = '#2D2D2D', na.value = 'purple3', 'Pack-Years Smoked\nPurple=NA')+
  #scale_color_gradientn(colors = rev(palfun(nrow(pdf))), limits = range(pdf$pack_years_smoked), 'Z-scored\nPack-Years')+
  #continuous_scale('color', scales::gradient_n_pal(pal), palette = 1)+
  #scale_color_viridis('Z-scaled\nPack-years', option = 'C', na.value = 'purple')+
  theme_light()


g_jk <- g_jk + theme_light() + scale_y_continuous(breaks = seq(-3,8))
g_jk_smoke <- g_jk_smoke + theme_light()+ scale_y_continuous(breaks = seq(-3,8))
#smoking plot: all lifelong nn-smkers have no pack-years.

pdffile <- paste0(outdir, '/plots_updated_THRESHOLDED.pdf')
pdf(pdffile, height=7, width=15)
g_jk + g_jk_smoke
dev.off()









### multivar regression with raw gene list



#raw gene list module score corplot
finalgenelist <- generes[, "HGNC.symbol"]

m_opt <- scale(modulescore(gem, finalgenelist))

md_save <- md


#recode some stuff...

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
                  ModuleScore = m_opt,
                  cpe = md$CPE,
                  pack_years_smoked = md$pack_years_smoked,
                  mutation_burden = mut$total,
                  
                  ajcc_pathologic_t = md$ajcc_pathologic_t_recoded ,
                  ajcc_pathologic_n = forcats::fct_infreq(md$ajcc_pathologic_n) ,
                  ajcc_pathologic_m = forcats::fct_infreq(md$ajcc_pathologic_m_recoded) ,
                  site_of_resection_or_biopsy_recoded = forcats::fct_infreq(md$site_of_resection_or_biopsy_recoded) ,
                  gender = forcats::fct_infreq(md$gender) 
                  
                  
)





#table 1 stuff
meansd <- function(x){
  message('Mean (SD) ', round(mean(x), 2), ' (', round(sd(x), 2), ')')
}

medianiqr <- function(x){
  message('Median [IQR] ', round(summary(x)[3], 2), 
          ' [', round(summary(x)[2], 2), ' - ',
          round(summary(x)[5],2), ']')
}

tabperc <- function(x){
  
  df <- data.frame(tab = as.numeric(table(x)),
                   perc = as.numeric(round( table(x) / length(na.omit(x))* 100, 2)) 
  )
  
  df[,3] <- paste0(df[,1], ' (', df[,2], ')')
  
  return(as.data.frame(df[,3], row.names = names(table(x))))
  
}


meansd(pdf$Age)
meansd(pdf$cpe)
medianiqr(pdf$pack_years_smoked)
medianiqr(pdf$mutation_burden)

tabperc(pdf$ajcc_pathologic_t)
tabperc(pdf$ajcc_pathologic_n)
tabperc(pdf$ajcc_pathologic_m)

tabperc(pdf$site_of_resection_or_biopsy_recoded)
tabperc(pdf$gender)
#tabperc(pdf$race_recoded)






#tab 2 bivar assoc
cor.test(pdf$ModuleScore, pdf$Age)



#repeat it but try to get standardized betas
pdf_std <- data.frame(Age = scale(md$age_at_index),
                      ModuleScore = m_opt,
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
m <- lm(data = pdf, ModuleScore ~ .)
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
m_std <- lm(data = pdf_std, ModuleScore ~ .)
summary(m_std)


#prep to save; convert model results to df and add std be
m_save <- broom::tidy(m)
m_std <- broom::tidy(m_std)

m_save$estimate_standardized <- m_std$estimate


mvrfile <- paste0(outdir, '/multivariableregressionresults_rawgenelist.csv')
m_save$term <- gsub(',', '', m_save$term)
write.csv(x= m_save, file = mvrfile, quote = F, row.names = F)











### multivar regression with nonthresholded module score



#jackknifed module score corplot (no thresholding)
finalgenelist <- generes[generes$jackknifeclassification=='Signal', "HGNC.symbol"]

m_opt <- scale(modulescore(gem, finalgenelist))

md_save <- md


#recode some stuff...

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
                  ModuleScore = m_opt,
                  cpe = md$CPE,
                  pack_years_smoked = md$pack_years_smoked,
                  mutation_burden = mut$total,
                  
                  ajcc_pathologic_t = md$ajcc_pathologic_t_recoded ,
                  ajcc_pathologic_n = forcats::fct_infreq(md$ajcc_pathologic_n) ,
                  ajcc_pathologic_m = forcats::fct_infreq(md$ajcc_pathologic_m_recoded) ,
                  site_of_resection_or_biopsy_recoded = forcats::fct_infreq(md$site_of_resection_or_biopsy_recoded) ,
                  gender = forcats::fct_infreq(md$gender) 
                  
                  
)





#table 1 stuff
meansd <- function(x){
  message('Mean (SD) ', round(mean(x), 2), ' (', round(sd(x), 2), ')')
}

medianiqr <- function(x){
  message('Median [IQR] ', round(summary(x)[3], 2), 
          ' [', round(summary(x)[2], 2), ' - ',
          round(summary(x)[5],2), ']')
}

tabperc <- function(x){
  
  df <- data.frame(tab = as.numeric(table(x)),
                   perc = as.numeric(round( table(x) / length(na.omit(x))* 100, 2)) 
  )
  
  df[,3] <- paste0(df[,1], ' (', df[,2], ')')
  
  return(as.data.frame(df[,3], row.names = names(table(x))))
  
}


meansd(pdf$Age)
meansd(pdf$cpe)
medianiqr(pdf$pack_years_smoked)
medianiqr(pdf$mutation_burden)

tabperc(pdf$ajcc_pathologic_t)
tabperc(pdf$ajcc_pathologic_n)
tabperc(pdf$ajcc_pathologic_m)

tabperc(pdf$site_of_resection_or_biopsy_recoded)
tabperc(pdf$gender)
#tabperc(pdf$race_recoded)






#tab 2 bivar assoc
cor.test(pdf$ModuleScore, pdf$Age)



#repeat it but try to get standardized betas
pdf_std <- data.frame(Age = scale(md$age_at_index),
                      ModuleScore = m_opt,
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
m <- lm(data = pdf, ModuleScore ~ .)
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
m_std <- lm(data = pdf_std, ModuleScore ~ .)
summary(m_std)


#prep to save; convert model results to df and add std be
m_save <- broom::tidy(m)
m_std <- broom::tidy(m_std)

m_save$estimate_standardized <- m_std$estimate


mvrfile <- paste0(outdir, '/multivariableregressionresults_jackknife.csv')
m_save$term <- gsub(',', '', m_save$term)
write.csv(x= m_save, file = mvrfile, quote = F, row.names = F)







## mvr with thresholded

#jackknifed module score corplot with thresholding
finalgenelist <- generes[generes$thresholdedclassification=='Signal', "HGNC.symbol"]

m_opt <- scale(modulescore(gem, finalgenelist))

md_save <- md


#recode some stuff...

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
                  ModuleScore = m_opt,
                  cpe = md$CPE,
                  pack_years_smoked = md$pack_years_smoked,
                  mutation_burden = mut$total,
                  
                  ajcc_pathologic_t = md$ajcc_pathologic_t_recoded ,
                  ajcc_pathologic_n = forcats::fct_infreq(md$ajcc_pathologic_n) ,
                  ajcc_pathologic_m = forcats::fct_infreq(md$ajcc_pathologic_m_recoded) ,
                  site_of_resection_or_biopsy_recoded = forcats::fct_infreq(md$site_of_resection_or_biopsy_recoded) ,
                  gender = forcats::fct_infreq(md$gender) 
                  
                  
)





#table 1 stuff
meansd <- function(x){
  message('Mean (SD) ', round(mean(x), 2), ' (', round(sd(x), 2), ')')
}

medianiqr <- function(x){
  message('Median [IQR] ', round(summary(x)[3], 2), 
          ' [', round(summary(x)[2], 2), ' - ',
          round(summary(x)[5],2), ']')
}

tabperc <- function(x){
  
  df <- data.frame(tab = as.numeric(table(x)),
                   perc = as.numeric(round( table(x) / length(na.omit(x))* 100, 2)) 
  )
  
  df[,3] <- paste0(df[,1], ' (', df[,2], ')')
  
  return(as.data.frame(df[,3], row.names = names(table(x))))
  
}


meansd(pdf$Age)
meansd(pdf$cpe)
medianiqr(pdf$pack_years_smoked)
medianiqr(pdf$mutation_burden)

tabperc(pdf$ajcc_pathologic_t)
tabperc(pdf$ajcc_pathologic_n)
tabperc(pdf$ajcc_pathologic_m)

tabperc(pdf$site_of_resection_or_biopsy_recoded)
tabperc(pdf$gender)
#tabperc(pdf$race_recoded)






#tab 2 bivar assoc
cor.test(pdf$ModuleScore, pdf$Age)



#repeat it but try to get standardized betas
pdf_std <- data.frame(Age = scale(md$age_at_index),
                      ModuleScore = m_opt,
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
m <- lm(data = pdf, ModuleScore ~ .)
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
m_std <- lm(data = pdf_std, ModuleScore ~ .)
summary(m_std)


#prep to save; convert model results to df and add std be
m_save <- broom::tidy(m)
m_std <- broom::tidy(m_std)

m_save$estimate_standardized <- m_std$estimate



mvrfile <- paste0(outdir, '/multivariableregressionresults_thresholded.csv')
m_save$term <- gsub(',', '', m_save$term)
write.csv(x= m_save, file = mvrfile, quote = F, row.names = F)






beepr::beep()

