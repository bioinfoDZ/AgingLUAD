library(biomaRt)
library(tidyverse)



#genelist
mousedeg <- read.csv('finalized_start-to-end-2021.02.20/A-SPCCre_vs_Y-SPCCre.csv', header = T)

colnames(mousedeg)[1] <- 'MGI.symbol'

#gene selection
#padj < 0.05
genelist <- mousedeg[mousedeg$p_val_adj <= 0.05,]


deglist <- list(genelist[genelist$avg_logFC > 0,],
                genelist[genelist$avg_logFC < 0,] )

genelist <- list(genelist[genelist$avg_logFC > 0,'MGI.symbol'],
                 genelist[genelist$avg_logFC < 0,'MGI.symbol'] )

names(genelist) <- c('A-SPCCre_vs_Y-SPCCre', 'Y-SPCCre_vs_A-SPCCre')



#mouse to human gene names, function def


human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mousetohuman <- function(mousegenes){
  
  #based on this:
  #https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/

  
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = mousegenes , mart = mouse, 
                   attributesL = c("hgnc_symbol", "ensembl_gene_id"), 
                   martL = human, uniqueRows=T)
  
  genesV2
}



#convert to human
outdflist <- list()

for(subvec in c(1:length(genelist)) ){
  
  message(names(genelist)[subvec])
  
  mg <- mousetohuman(genelist[[subvec]])
  
  #need to figure out how to deal with duplicates; for now just remove.
  mg <- mg[!duplicated(mg$MGI.symbol),]
  mg <- mg[!duplicated(mg$HGNC.symbol),]
  
  #save outs...
  deg <- deglist[[subvec]]
  
  #sort by LFC
  deg <- deg[order(abs(deg$avg_logFC), decreasing = T),]

  #match human-gene ortholog df
  deg <- deg[deg$MGI.symbol %in% mg$MGI.symbol,]
  mg <- mg[match(deg$MGI.symbol, mg$MGI.symbol),]
  
  deg <- cbind(mg, deg[,-1])
  
  outdflist[[subvec]] <- deg
  
  
}


beepr::beep()


names(outdflist) <- names(genelist)

saveRDS(outdflist, 'finalized_start-to-end-2021.02.20/parsedmousedeglist.rds')



