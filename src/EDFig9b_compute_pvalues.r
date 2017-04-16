## This file is called from 'main_nci60_analysis.r'.
## This script is to generate 'Extended Data Fig. 9b', which visually show the 17 GOs
## were exclusively significantly enriched against the ferroptosis cluster.

## 

if(0){
  ## The following codes were run once to generate a object, 'slogp.fin.go'
  
  ## Enrichment Scores
  es.fin.go <- t(sapply(GOs.ge10[sig.gos],function(go){
    sapply(as.character(seq(53)),function(i){
      pog <- names(sort(corrs[used.genes.ge10,i],decreasing=T)) ## genes in descending order 
      pvals.ES.GSEA(pog,go,nperm=1)$ES
    })
  }))
  
  
  
  ## P-values with 10^3 iterations
  system.time(pval3 <- mclapply(GOs.ge10[sig.gos],function(go){ ## 340 secs
    sapply(as.character(seq(53)),function(i){
      pog <- names(sort(corrs[used.genes.ge10,i],decreasing=T)) ## genes in descending order
      pvals.ES.GSEA(pog,go,nperm=1000)$p.value
    })
  },mc.cores=(n.cores-1)))
  pval3 <- do.call(rbind,pval3)
  lab1 <- which(unlist(pval3)<1e-2,arr.ind=T)
  
  ## P-values with 10^4 iterations
  system.time(pval4 <- mclapply(seq(nrow(lab1)),function(j){ ## 385 secs
    cat("*")
    go <- GOs.ge10[[sig.gos[lab1[j,1]]]]
    col <- as.character(lab1[j,2])
    pog <- names(sort(corrs[used.genes.ge10,col],decreasing=T)) ## genes in descending order
    pvals.ES.GSEA(pog,go,nperm=10^4)$p.value
  },mc.cores=(n.cores-1)))
  pval4 <- unlist(pval4)
  lab2 <- lab1[which(unlist(pval4)<1e-3,arr.ind=T),]
  
  ## P-values with 10^5 iterations
  system.time(pval5 <- mclapply(seq(nrow(lab2)),function(j){
    cat("*")
    go <- GOs.ge10[[sig.gos[lab2[j,1]]]]
    col <- as.character(lab2[j,2])
    pog <- names(sort(corrs[used.genes.ge10,col],decreasing=T)) ## genes in descending order
    pvals.ES.GSEA(pog,go,nperm=10^5)$p.value
  },mc.cores=(n.cores-1)))
  pval5 <- unlist(pval5)
  lab3 <- lab2[which(unlist(pval5)<1e-4,arr.ind=T),]
  
  ## P-values with 10^6 iterations
  system.time(pval6 <- mclapply(seq(nrow(lab3)),function(j){
    cat("*")
    go <- GOs.ge10[[sig.gos[lab3[j,1]]]]
    col <- as.character(lab3[j,2])
    pog <- names(sort(corrs[used.genes.ge10,col],decreasing=T)) ## genes in descending order
    pvals.ES.GSEA(pog,go,nperm=10^6)$p.value
  },mc.cores=min(7,nrow(lab3))))
  pval6 <- unlist(pval6)
  lab4 <- lab3[which(unlist(pval6)<1e-5,arr.ind=T),]
  
  ##
  ## assemble computed P-values
  pval.fin.go <- pval3
  for(j in seq(nrow(lab1))){
    pval.fin.go[lab1[j,1],lab1[j,2]] <- pval4[[j]]
  }
  for(j in seq(nrow(lab2))){
    pval.fin.go[lab2[j,1],lab2[j,2]] <- pval5[[j]]
  }
  for(j in seq(nrow(lab3))){
    pval.fin.go[lab3[j,1],lab3[j,2]] <- pval6[[j]]
  }

  ## compute singed log P-values from esp.fin.go and es.
  slogp.fin.gos <- t(-log10(pval.fin.go)*sign(es.fin.go))
}

## Extended Data Fig. 9b
setwd(fig.dir)
pdf("EDFigure_9b.pdf",width=6,height=6)
##
o <- c("GO:0016491","GO:0009055","GO:0016614","GO:0016616","GO:0016651","GO:0006732",
       "GO:0006807","GO:0009308","GO:0005216","GO:0006812","GO:0030001","GO:0015672",
       "GO:0015674","GO:0002376","GO:0019221","GO:0007093","GO:0003729")
ordered2 <- ordered.sig.gos[o]
slogp.fin.gos <- slogp.fin.gos[,ordered2]
par(mar=c(10,5,4,1)+.1)
plot(range(0.5,17.5),range(slogp.fin.gos),type="n",xlab="",ylab="P-value",axes=F)
abline(v=1:17,col="grey80")
abline(h=0)
axis(2,seq(-6,6,2))
axis(1,at=1:17,labels=o,las=2,cex.axis=.7)
##
sapply(seq(17),function(i){
  points(jitter(rep(i,53),amount=.2),slogp.fin.gos[,i],pch=20,col=(seq(53)%in% 1)+1,cex=.3)
})
dev.off()

