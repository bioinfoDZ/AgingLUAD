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
outdir <- "finalized_start-to-end-2021.02.20/A_aged_kras_correlation/"
dir.create(outdir)


#select deg list; aged, young, etc
deglist <- readRDS('finalized_start-to-end-2021.02.20/parsedmousedeglist.rds')
deg <- deglist[[1]]




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





# #set parallelization...
cl <- parallel::makeForkCluster(numclust)
doParallel::registerDoParallel(cl)


#add module score for each genelist
#reslist <- list()
#for(listdex in c(1:length(genelist)) ) {
reslist <- foreach(listdex = c(1:length(genelist)) ) %dopar% {
  
  set.seed(2020)
  
  sublist <- genelist[-listdex]
  leaveoutgene <- genelist[listdex]
  
  pdf <- data.frame(Age = md$age_at_index,
                    ModuleScore = scale(modulescore(gem, sublist) ) )
  
  
  #get results
  cres <- cor.test(pdf[,1], pdf[,2])
  
  resdf <- data.frame(listdex = listdex,
                      leaveoutgene = leaveoutgene,
                      cor = cres$estimate,
                      p = cres$p.value,
                      row.names = NULL, stringsAsFactors = F
  )
  
  #return results, will be added to list
  resdf
  
}


doParallel::stopImplicitCluster()


# get whole gene list res

pdf <- data.frame(Age = md$age_at_index,
                  ModuleScore = scale(modulescore(gem, genelist)) )

cres <- cor.test(pdf[,1], pdf[,2])

resdf <- data.frame(listdex = 0,
                    leaveoutgene = "FULL GENE LIST",
                    cor = cres$estimate,
                    p = cres$p.value,
                    row.names = NULL, stringsAsFactors = F
)

#sort results and add in full gene list res
brd <- as.data.frame(data.table::rbindlist(reslist))
brd <- brd[order(brd$listdex),]
brd <- rbind(resdf, brd)

#repel label
top <- head( brd[order(brd$cor), "leaveoutgene"] )
top <- c(top, "FULL GENE LIST")
brd$repel <- NA
brd[brd$leaveoutgene %in% top,"repel"] <- brd[brd$leaveoutgene %in% top,"leaveoutgene"]

bottom <- tail( brd[order(brd$cor), "leaveoutgene"] )
brd[brd$leaveoutgene %in% bottom,"repel"] <- brd[brd$leaveoutgene %in% bottom,"leaveoutgene"]



jplot <- ggplot(brd, aes(x = listdex, y = cor))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = repel),
                           col = 'red',
                           box.padding = 2, force = 3, #nudge_y = -0.005
  )+
  geom_hline(yintercept = brd[1,'cor'], col = 'red', alpha = 0.5)+
  theme_linedraw()+
  xlab('Gene, sorted by magnitude of Mouse LogFC') +
  ylab('Signature Correlation with Age\nAfter Leaving out each gene')+
  labs(title = 'Jackknifed')







#uncomment to run all
#dirreslist[[dir]] <- reslist


#}


#dev.off()



#get final gene list, 
#leaving out the ones that make the correlation higher if excluded (noise genes)
# use absolute value to make direction agnostic; ie want away from null
#also remove the full gene list result from the df
brd_good <- brd[abs(brd$cor) <= abs(brd[1,'cor']),]
brd_good <- brd_good[-1,]


#save the jackknifing results
jkr <- brd[,2:4]
jkr$classification <- ifelse(abs(jkr$cor) <= abs(jkr[1,'cor']),
                             yes = 'Signal', no = 'Noise')
jkr[1,"classification"] <- NA




brd$select <- 'Remove'
brd[brd$leaveoutgene %in% brd_good$leaveoutgene, 'select'] <- 'Keep'

jplot_select <- ggplot(brd, aes(x = listdex, y = cor, col=select))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = repel),
                           col = 'black',
                           box.padding = 2, force = 3, #nudge_y = -0.005
  )+
  geom_hline(yintercept = brd[1,'cor'], col = 'red', alpha = 0.5)+
  theme_linedraw()+
  scale_color_manual(values = c('Red', 'grey'))+
  xlab('Gene, sorted by magnitude of Mouse LogFC') +
  ylab('Signature Correlation with Age\nAfter Leaving out each gene')+
  labs(title = 'Jackknifed')





#thresholded jacknifing results plot, showing 1SD below mean jk score
jkr_tmp <- jkr[-1,]
corsd <- sd(jkr$cor)

#to get threshold in direction agnostic manner, multiply SD by full gene list module score.
thres <- mean(jkr_tmp$cor) + (corsd * -sign(jkr[1,'cor']) )

brd$select_thres <- 'Remove'
brd[abs(brd$cor) <= abs(thres) , 'select_thres'] <- 'Keep'




jplot2 <- ggplot(brd, aes(x = listdex, y = cor, col=select_thres))+
  geom_point()+
  #geom_line()+
  ggrepel::geom_text_repel(aes(label = repel),
                           col = 'black',
                           box.padding = 2, force = 3, #nudge_y = -0.005
  )+
  geom_hline(yintercept = brd[1,'cor'], col = 'red', alpha = 0.5, size=1)+
  geom_hline(yintercept = mean(jkr_tmp$cor), col = 'blue',alpha = 0.5, size=1.5)+
  geom_hline(yintercept = mean(jkr_tmp$cor)+corsd, col = 'blue', linetype='dashed', alpha = 0.5, size=1.5)+
  geom_hline(yintercept = mean(jkr_tmp$cor)-corsd, col = 'blue', linetype='dashed', alpha = 0.5, size=1.5)+
  theme_linedraw()+
  scale_color_manual(values = c('Red', 'grey'))+
  xlab('Gene, sorted by magnitude of Mouse LogFC') +
  ylab('Signature Correlation with Age\nAfter Leaving out each gene')+
  labs(title = "Jackknifed + Thresholded")






jplot
jplot_select
jplot2



#jackknifed module score corplot (no thresholding)
finalgenelist <- brd_good$leaveoutgene

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
  theme_linedraw() 



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

#cor of smoke with module score
cres_smoke <- cor.test(pdf$ModuleScore, pdf$pack_years_smoked)
nicep_smoke <- ifelse(test = cres_smoke$p.value < 0.01, formatC(cres_smoke$p.value, format='e', digits=2), 
                      no = round(cres_smoke$p.value, 2))



g_jk_smoke <- ggplot(pdf, aes(x = Age, y = ModuleScore, col = pack_years_smoked )) +
  geom_point(size = 3, aes(shape = smoker_status_cat))+
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
  theme_linedraw()


g_jk
g_jk_smoke
#smoking plot: all lifelong nn-smkers have no pack-years.









#also get unoptimized approach
m_raw <- scale(modulescore(gem, genelist))

pdf <- data.frame(Age = md$age_at_index,
                  ModuleScore = m_raw)


cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))


g0 <- ggplot(pdf, aes(x = Age, y = ModuleScore))+
  geom_point(size=3)+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Raw Module Score",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         #  '\nNum final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3),
                         '\nP (cor.test) = ', nicep))+
  ylab('Z-Scaled Module Score')+
  theme_linedraw() 

g0




#jackknife plot with smoke status.
pdf <- data.frame(ModuleScore = m_raw,
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

#cor of smoke with module score
cres_smoke <- cor.test(pdf$ModuleScore, pdf$pack_years_smoked)
nicep_smoke <- ifelse(test = cres_smoke$p.value < 0.01, formatC(cres_smoke$p.value, format='e', digits=2), 
                      no = round(cres_smoke$p.value, 2))



g0_smoke <- ggplot(pdf, aes(x = Age, y = ModuleScore, col = pack_years_smoked )) +
  geom_point(size = 3, aes(shape = smoker_status_cat))+
  scale_shape_manual(values = c(16,17,3), 'History of Smoking')+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Raw Module Score + Smoking Status",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3), ' ; Smoke cor with modulescore = ',round(cres_smoke$estimate,3),
                         '\nP (cor.test) = ', nicep, ' ; Smoke P = ',nicep_smoke))+
  ylab('Z-Scaled Module Score')+
  #scale_color_distiller(palette = 'Greys', 'Pack-years', direction = 1, na.value = 'Purple')+
  scale_color_gradient(low = '#B4B4B4', high = '#2D2D2D', na.value = 'purple3', 'Pack-Years Smoked\nPurple=NA')+
  #scale_color_gradientn(colors = rev(palfun(nrow(pdf))), limits = range(pdf$pack_years_smoked), 'Z-scored\nPack-Years')+
  #continuous_scale('color', scales::gradient_n_pal(pal), palette = 1)+
  #scale_color_viridis('Z-scaled\nPack-years', option = 'C', na.value = 'purple')+
  theme_linedraw()


g0
g0_smoke
#smoking plot: all lifelong nn-smkers have no pack-years.






#to get threshold in direction agnostic manner, multiply SD by full gene list module score.
thres <- mean(jkr_tmp$cor) + (corsd * -sign(jkr[1,'cor']) )

brd$select_thres <- 'Remove'
brd[abs(brd$cor) <= abs(thres) , 'select_thres'] <- 'Keep'


#thresholded

#to get threshold in direction agnostic manner, multiply SD by full gene list module score.
thres <- mean(jkr_tmp$cor) + (corsd * -sign(jkr[1,'cor']) )


jkr_thresholded <- jkr_tmp[abs(jkr_tmp$cor) <= abs(thres) ,]




finalgenelist <- jkr_thresholded$leaveoutgene

m_opt <- scale(modulescore(gem, finalgenelist))

pdf <- data.frame(Age = md$age_at_index,
                  ModuleScore = m_opt)

cres <- cor.test(pdf[,1], pdf[,2])

nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                no = round(cres$p.value, 2))

g_thres <- ggplot(pdf, aes(x = Age, y = ModuleScore))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Jackknifed And Thresholded Module Score",
       subtitle = paste0('N = ', ncol(gem),
                         '\nNum genes precut= ', length(rawgenelist),
                         ' ; Num final genes = ', length(finalgenelist),
                         '\nPearson cor r = ', round(cres$estimate, 3),
                         '\nP (cor.test) = ', nicep))+
  ylab('Z-Scaled Module Score')+
  theme_linedraw() 


g_thres





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

#cor of smoke with module score
cres_smoke <- cor.test(pdf$ModuleScore, pdf$pack_years_smoked)
nicep_smoke <- ifelse(test = cres_smoke$p.value < 0.01, formatC(cres_smoke$p.value, format='e', digits=2), 
                      no = round(cres_smoke$p.value, 2))



g_thres_smoke <- ggplot(pdf, aes(x = Age, y = ModuleScore, col = pack_years_smoked )) +
  geom_point(size = 3, aes(shape = smoker_status_cat))+
  scale_shape_manual(values = c(16,17,3), 'History of Smoking')+
  geom_smooth(method = 'lm', formula = y~x, se = T)+
  labs(title = "Jackknifed & Threshollded Module Score + Smoke",
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
  theme_linedraw()






dev.off()


#plots
# 1) raw module score correlation plot with age bipot+linear regression
# 2) raw module score correlation plot with age bipot+linear regression, colored by smoke status
# 
# 3) jackknife result plot
# 4) jackknife result plot showing which genes are retained
# 5) jackknifed module score correlation plot
# 6) jackknifed module score correlation plot + smoke status
# 
# 7) jackknife result plot with 1std.dev threshold illustrated
# 8) jacknife+thresholding modue score correlation plot
# 9) jacknife+thresholding modue score correlation plot + smoke status

g0
g0_smoke

jplot
jplot_select
g_jk
g_jk_smoke

jplot2
g_thres
g_thres_smoke 



#save thresholded and unthresholded genes
jacknifesummary <- jkr_tmp[,c(1,4)]
jacknifesummary$thresholdedclassification <- 'Noise'
jacknifesummary[jacknifesummary$leaveoutgene %in% jkr_thresholded$leaveoutgene, "thresholdedclassification"] <- 'Signal'

colnames(jacknifesummary) <- c('HGNC.symbol', "jackknifeclassification", 'thresholdedclassification')




head( jacknifesummary )




### multivar regression with thresholded module score


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



beepr::beep()



### write out ###


# output files
# 1) jackknife results; thresholded jacknife results
# 2) mvr results (using Y = age and X =  theresholded module score, other X vars too)
# 3) pdf with plots.


jksfile <- paste0(outdir, '/jackknifesummary.csv')
write.csv(jacknifesummary, jksfile, quote = F, row.names = F)

mvrfile <- paste0(outdir, '/multivariableregressionresults.csv')
m_save$term <- gsub(',', '', m_save$term)
write.csv(x= m_save, file = mvrfile, quote = F, row.names = F)


pdffile <- paste0(outdir, '/plots.pdf')

pdf(pdffile, height=10, width=10)
g0
g0_smoke

jplot
jplot_select
g_jk
g_jk_smoke

jplot2
g_thres
g_thres_smoke
dev.off()


sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] olsrr_0.5.3         TCGAbiolinks_2.18.0 doParallel_1.0.16   iterators_1.0.13   
# [5] foreach_1.5.1       biomaRt_2.46.0      forcats_0.5.0       stringr_1.4.0      
# [9] dplyr_1.0.4         purrr_0.3.4         readr_1.4.0         tidyr_1.1.2        
# [13] tibble_3.0.6        ggplot2_3.3.3       tidyverse_1.3.0    
# 
# loaded via a namespace (and not attached):
#   [1] reticulate_1.18             R.utils_2.10.1              tidyselect_1.1.0           
# [4] RSQLite_2.2.1               AnnotationDbi_1.52.0        htmlwidgets_1.5.3          
# [7] grid_4.0.3                  BiocParallel_1.24.1         Rtsne_0.15                 
# [10] munsell_0.5.0               preprocessCore_1.52.0       units_0.6-7                
# [13] codetools_0.2-18            ica_1.0-2                   future_1.21.0              
# [16] miniUI_0.1.1.1              withr_2.4.1                 audio_0.1-7                
# [19] colorspace_2.0-0            Biobase_2.50.0              knitr_1.30                 
# [22] rstudioapi_0.13             Seurat_4.0.0                stats4_4.0.3               
# [25] ROCR_1.0-11                 tensor_1.5                  listenv_0.8.0              
# [28] labeling_0.4.2              MatrixGenerics_1.2.0        GenomeInfoDbData_1.2.4     
# [31] polyclip_1.10-0             farver_2.0.3                bit64_4.0.5                
# [34] downloader_0.4              parallelly_1.23.0           vctrs_0.3.6                
# [37] generics_0.1.0              lambda.r_1.2.4              xfun_0.19                  
# [40] BiocFileCache_1.14.0        R6_2.5.0                    GenomeInfoDb_1.26.2        
# [43] locfit_1.5-9.4              bitops_1.0-6                spatstat.utils_2.0-0       
# [46] DelayedArray_0.16.0         assertthat_0.2.1            promises_1.2.0.1           
# [49] scales_1.1.1                gtable_0.3.0                globals_0.14.0             
# [52] goftest_1.2-2               beepr_1.3                   rlang_0.4.10               
# [55] genefilter_1.72.0           splines_4.0.3               lazyeval_0.2.2             
# [58] broom_0.7.2                 modelr_0.1.8                reshape2_1.4.4             
# [61] abind_1.4-5                 backports_1.2.0             httpuv_1.5.5               
# [64] tools_4.0.3                 ellipsis_0.3.1              RColorBrewer_1.1-2         
# [67] BiocGenerics_0.36.0         ggridges_0.5.3              Rcpp_1.0.6                 
# [70] plyr_1.8.6                  progress_1.2.2              zlibbioc_1.36.0            
# [73] classInt_0.4-3              RCurl_1.98-1.2              prettyunits_1.1.1          
# [76] rpart_4.1-15                openssl_1.4.3               deldir_0.2-9               
# [79] pbapply_1.4-3               cowplot_1.1.1               S4Vectors_0.28.1           
# [82] zoo_1.8-8                   SeuratObject_4.0.0          SummarizedExperiment_1.20.0
# [85] haven_2.3.1                 ggrepel_0.9.1               cluster_2.1.0              
# [88] fs_1.5.0                    magrittr_2.0.1              data.table_1.13.6          
# [91] futile.options_1.0.1        scattermore_0.7             openxlsx_4.2.3             
# [94] reprex_0.3.0                lmtest_0.9-38               RANN_2.6.1                 
# [97] fitdistrplus_1.1-3          matrixStats_0.58.0          hms_0.5.3                  
# [100] patchwork_1.1.1             TCGAbiolinksGUI.data_1.10.0 mime_0.10                  
# [103] xtable_1.8-4                XML_3.99-0.5                VennDiagram_1.6.20         
# [106] rio_0.5.16                  readxl_1.3.1                IRanges_2.24.0             
# [109] gridExtra_2.3               compiler_4.0.3              KernSmooth_2.23-18         
# [112] crayon_1.4.1                R.oo_1.24.0                 htmltools_0.5.1.1          
# [115] mgcv_1.8-33                 later_1.1.0.1               ggVennDiagram_0.3          
# [118] geneplotter_1.68.0          lubridate_1.7.9.2           DBI_1.1.0                  
# [121] formatR_1.7                 dbplyr_2.0.0                MASS_7.3-53                
# [124] rappdirs_0.3.3              sf_0.9-6                    Matrix_1.2-18              
# [127] car_3.0-10                  cli_2.3.0                   R.methodsS3_1.8.1          
# [130] igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3            
# [133] foreign_0.8-80              plotly_4.9.3                xml2_1.3.2                 
# [136] annotate_1.68.0             XVector_0.30.0              maftools_2.6.0             
# [139] rvest_0.3.6                 digest_0.6.27               sctransform_0.3.2          
# [142] RcppAnnoy_0.0.18            spatstat.data_2.0-0         cellranger_1.1.0           
# [145] leiden_0.3.7                nortest_1.0-4               uwot_0.1.10                
# [148] curl_4.3                    shiny_1.6.0                 lifecycle_1.0.0            
# [151] nlme_3.1-151                jsonlite_1.7.2              carData_3.0-4              
# [154] futile.logger_1.4.3         viridisLite_0.3.0           askpass_1.1                
# [157] pillar_1.4.7                lattice_0.20-41             fastmap_1.1.0              
# [160] httr_1.4.2                  survival_3.2-7              glue_1.4.2                 
# [163] zip_2.1.1                   spatstat_1.64-1             png_0.1-7                  
# [166] bit_4.0.4                   class_7.3-17                stringi_1.5.3              
# [169] blob_1.2.1                  DESeq2_1.30.0               memoise_1.1.0              
# [172] e1071_1.7-4                 irlba_2.3.3                 future.apply_1.7.0 


