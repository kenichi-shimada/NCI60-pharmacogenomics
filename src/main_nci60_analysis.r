#=================================================================================
# Title:        Model-based clustering of GI50 profiles of NCI60 cell lines
# Author:       Kenichi Shimada
# Date:         18 Sep 2014
# R-version:    2.15.2 (2012-10-26)
#=================================================================================

## This script should be run inside the 'nci60_analysis' directory.
## Each figure was already made and put in ./figs directory.

## working directory
dir <- getwd()
fig.dir <- paste(dir,"/figs",sep="")

## Load libraries -------------------------------------------------------------------------
## if not installed already, install from CRAN or Bioconductor
library(gplots) ## heatmap.2 function
library(RColorBrewer) ## gradation of colors in heatmap
library(mclust) ## model-based clustering
library(robustbase) ## lmrob for robust linear regression
library(org.Hs.eg.db) ## Gene ID <=> Gene SYMBOL
library(parallel) ## mclapply to use multiple cores 

## number of cores available in the machine
n.cores <- parallel:::detectCores() ## number of cores (this was originally run in an 8-cores machine)
## multicore is depricated, so if

## load r functions to use in this analysis
source("rfunctions.r")

## Load objects & trimming data?-----------------------------------------------------------
setwd(dir)
load("nci60-analysis.rda")
## this file contains the following R objects:

## gi: logGI50 values of 75446 (50839 unique) compounds in 59 cell lines

## gi.fin: logGI50 values of 21 (10 unique) compounds in 59 cell lines
## tissues: cell line name and tissue origins
## log.txn.ge80: RMA-transformed microarray data of 59 cell lines, selected genes within top 80% IQRs
## corrs: matrix containing the Spearman correlation between
## es.pvals: Each pathway's enrichment score and the p-value computed using GSEA.
## msigdb: MsigDB pathway data (version 3.1)
## slogp.fin.gos: signed log P-values of ferroptosis-associted GOs, to generate Extended Data Fig. 9b
####################################################################################################
## Computing GI50 profiles from public repository ##################################################
####################################################################################################
cmpds <- sub("\\..+","",rownames(gi)) ## 75446 compounds
uniq.cmpds <- unique(cmpds) ## 50839 unique compounds

## Compute unique GI50 profiles by taking the median of the replicates
system.time({ ## 222 secs
  med.gi <- t(sapply(uniq.cmpds,function(cmpd){
    ids <- which(cmpds %in% cmpd)
    if(length(ids)==1){
      x <- gi[ids,]
      med.x <- x-median(x,na.rm=T)
      return(med.x)
    }else{
      v <- apply(apply(gi[ids,],1,function(x)x <- x - median(x,na.rm=T)), ## median subtraction
                 1,function(y)median(y,na.rm=T)) ## median of replicates
      med.v <- v - median(v,na.rm=T)
      return(med.v)
    }
  }))
})

## Select GI50 profiles without missing values and whose IQR>0 : 6249 compounds
n.tested.cells <- apply(med.gi,1,function(x)sum(!is.na(x)))
iqrs <- apply(med.gi, 1, function(x)diff(quantile(x,c(.25,.75),na.rm=T))) #interquartile range (IQR)
used.gi <- med.gi[n.tested.cells == 59&iqrs>0,] ## 6249 compounds were selected

####################################################################################################
## Computing GI50 profiles of ferroptosis inducers #################################################
####################################################################################################
uniq.fin <- unique(sub("\\..+","",rownames(gi.fin))) ## 10 compounds

## computing a unique GI50 profile for each ferroptosis inducer
used.gi.fin <- t(sapply(uniq.fin,function(ur){
  thism <- apply(gi.fin[grep(paste("^",ur,"\\.",sep=""),rownames(gi.fin)),],
                 1,function(x)x-median(x,na.rm=T))
  x <- apply(thism,1,median,na.rm=T)
  return(x)
}))

## replace a missing value with an estimate
used.gi.fin["MEII","LC:EKVX"] <- mean(used.gi.fin[c("AE","PE","erastin"),"LC:EKVX"])

####################################################################################################
## Selecting cell line-selective compounds #########################################################
####################################################################################################
## IQRs of 6249 compounds and 10 ferroptosis inducers
used.iqrs <- apply(used.gi, 1, function(x)diff(quantile(x,c(.25,.75),na.rm=T))) ## 6249
used.iqrs.fin <- sort(apply(used.gi.fin, 1, function(x)diff(quantile(x,c(.25,.75),na.rm=T)))) ## 10

if(0){## Fig. 2a: barplot of cell-line-selectivity with in ferroptosis inducers
  setwd(fig.dir)
  pdf("figure_2a.pdf",width=3,height=4)
  mir <- used.iqrs.fin##;names(mir)[8] <- "RSL-CIL56"
  fi.class <- list(gsh=c("erastin","AE","PE","MEII"),
                   gpx4=c("RSL3","1S3R","ML162","ML210"),
                   cil56=c("CIL56","FIN56"))
  fi.mir <- lapply(fi.class[1:3],function(x)mir[x])
  cols <- 1##brewer.pal(3,"Set2")
  par(mar=c(7,5,1,1)+.1)
  plot(c(1,3),range(c(0,unlist(fi.mir))),type="n",xlim=c(.8,4),
       axes=F,xlab="",ylab="Cell-line selectivity")
  abline(v=1:3,col="grey80")
  par(xpd=T)
  set.seed(123)
  sapply(1:3,function(i){
    mir1 <- fi.mir[[i]]
    xs <- jitter(rep(i,length(mir1)),amount=.15)
    points(xs,mir1,col=1,pch=20,cex=1)
    ##text(xs,mir1,names(mir1),pos=4)
  })
  par(xpd=F)
  axis(1,at=1:3,labels=c("GSH depletors","GPX4 inhibitors","CIL56 analogs"),las=2)
  axis(2)
  abline(h=min(unlist(fi.mir[1:2])),col=2,lty=2)
  ## boxplot
  dev.off()
}

if(0){## Fig. 4a: boxplot of cell-line selectivity (IQR of each GI50 profile)
  setwd(fig.dir)
  pdf("figure_4a.pdf",width=3,height=4)
  par(mar=c(7,4,1,7)+.1)
  mir <- used.iqrs.fin##;names(mir)[8] <- "RSL-CIL56"
  mgi50 <-  list(mir,used.iqrs)
  a <- boxplot(mgi50,border=NA,cex=0,axes=F)
  axis(2)
  axis(1,at=1:2,labels=c("Others","Ferroptosis"),las=2)
  set.seed(1245)
  points(jitter(rep(1,length(mgi50[[2]])),amount=.35),mgi50[[2]],pch=20,cex=.15,col="grey70")
  boxplot(mgi50[2:1],cex=0,axes=F,add=T,border="grey40")
  set.seed(543)
  cols <- brewer.pal(6,"Set1")
  iqrs.fins.cols <- "grey20"##cols[c(2,2,2,2,3,3,3,1,3)]
  abline(h=min(mgi50[[1]][-1]),col=2,lty=2)
  points(jitter(rep(2,length(mgi50[[1]])),amount=.2),mgi50[[1]],pch=20,cex=.7,col=iqrs.fins.cols)
  dev.off()
}

## Selecting cell line-selective compounds
has.large.iqr <- used.iqrs > min(used.iqrs.fin[-1]) ## 2555 compounds
sele.gi <- rbind(used.gi[has.large.iqr,],used.gi.fin) ## GI50 values of 2565 cell-line-selective lethals
####################################################################################################
## Model-based clustering of GI50 profiles of cell line-selectie comopunds #########################
####################################################################################################
## find the optimal model and # of clusters
n.clusts.range <-  30:70
system.time({## 244 secs
  bics <- mclustBIC(sele.gi,G=n.clusts.range,modelNames=c("EII","VII"),
                    prior=priorControl())
})
max.ind <- which(bics==max(bics),arr.ind=T)
model <- colnames(bics)[max.ind[2]] ## optimal model: VII
g <- as.numeric(rownames(bics)[max.ind[1]]) ## optimal # of clusters: 53
##
if(0){
  ## Figure choose the optimal number of clusters
  pdf("choose_n_clusters_with_BIC.pdf")
  plot(n.clusts.range,bics[,model],type="b",
       xlab="# of clusters",ylab="Bayesian Information Criterion",
       main=paste("model:",model))
  abline(v=g,col=2)
  dev.off()
}
system.time(mc <- Mclust(sele.gi,G=g,modelNames=model,prior=priorControl())) ## 47 secs

## Discover the ferroptosis cluster
mem <- mc$classification ## compounds' memberships
mem.fins <- mem[uniq.fin]
fin.clust <- names(sort(table(mem.fins),decreasing=T))[1] ## ferroptosis cluster is "1"-st cluster

## GI50 profile of each cluster
gi.clust <- sapply(seq(g),function(m)apply(sele.gi[mem==m,],2,median)) ## gi50 profiles of the 53 clusters
colnames(gi.clust) <- as.character(seq(g))

if(0){
  ## Fig 4b heatmap
  setwd(fig.dir)
  pdf("figure_4b_heatmap.pdf",width=7,height=5) ##fig.4b
  colors <- c("white","red")[seq(as.numeric(g)) %in% fin.clust+1]
  xx <- myheatmap(gi.clust,cexRow=.5,cexCol=.5,keysize=1.5,RowSideColors=colors)
  dev.off()
  ## Fig 4b barplot
  pdf("figure_4b_barplot.pdf",width=2,height=5) ##fig.4b
  tab.fins <- table(mem,rep(0:1,c(2555,10)));rownames(tab.fins) <- as.character(1:53)
  sort.tab.fins <- tab.fins[xx$rowInd,]
  barplot(t(sort.tab.fins[,2:1]),beside=F,las=2,col=c(2,"lightblue"),horiz=T,cex.names=.3,border=NA)
  dev.off()
  ##
  pdf("figure_4c.pdf",width=5,height=4) ##fig.4c
  gi50.within.fi.clust <- t(sele.gi[mem==fin.clust,])
  myheatmap(gi50.within.fi.clust[xx$colInd,],cexRow=.7,cexCol=.5,keysize=2,
            Colv=FALSE,dendrogram="row")
  dev.off()
}

##
if(0){
  setwd(dir)
  write.csv(gi,"logGI50_all.csv")
  write.csv(used.gi.fin,"logGI50_FIN.csv")
  write.csv(cbind(mem,sele.gi),"cell-line.selective.csv")
}
####################################################################################################
## Transcriptome analysis of the ferroptosis cluster ###############################################
####################################################################################################
## transcriptome data are available at 
## https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data
## used GENELOGIC's U133 array data (RMA-normalized)
## Order of cell lines should be the sames between GI50 matrix and transcription matrix ---

t.gi50 <- t(gi.clust) ##transpose GI50 profiles
colnames(t.gi50) <- sub(":",".",colnames(t.gi50))

txn <- log.txn.ge80[,colnames(t.gi50)] ## transcriptional expression
identical(colnames(t.gi50),colnames(txn)) ## TRUE
GOs <- msigdb[grep("^GO",names(msigdb))] ## 1454 GOs

## genes that are in both transcriptome (txn) and pathway database (msigdb) are tested
txns <- rownames(txn)
genes <- unique(unlist(msigdb)) ## all the genes annotated in msigdb
used.genes <- genes[genes %in% txns] ## 10850 genes' expression is measured
used.msigdb <- lapply(msigdb,function(x)x[x %in% used.genes]) ## removed genes not measured from msigdb
msigdb.ge10 <- used.msigdb[sapply(used.msigdb,length)>=10] ##2763 gene sets contain 10 genes or more
used.genes.ge10 <- unique(unlist(msigdb.ge10)) ## 10850 genes in the above gene sets
GOs.ge10 <- msigdb.ge10[grep("^GO",names(msigdb.ge10))] ## 1222 gene sets are Gene Ontology terms

## Spearman correlation btwn each gene expression and each GI50 profile among 59 cell lines
## 10850 genes x 53 compound clusters ## 152 secs with n.cores=8
system.time(corrs <- mclapply(rownames(t.gi50),function(cmpd){
  sapply(used.genes.ge10,function(gid){
    cor(txn[gid,],t.gi50[cmpd,],method="spearman",use="complete.obs")
  })
},mc.cores=(n.cores-1)))
corrs <- do.call(cbind,corrs)
colnames(corrs) <- rownames(t.gi50)
rownames(corrs) <- used.genes.ge10
pos.ordered.genes <- names(sort(corrs[used.genes.ge10,fin.clust],decreasing=T))

## scatterplots (txn.exppression vs GI50) of the largest and smallest correlations
if(0){
  ## Extended Data Fig. 9a
  setwd(fig.dir)
  pdf("EDfigure_9a_hist.pdf",width=5,height=5)
  hist(corrs[,fin.clust],pch=20,cex=.1,axes=T,ylab="Spearman Coef.(mRNA vs GI50)",breaks=30,col="grey80")
  abline(v=0)
  dev.off()
  
  max.corr.txn <- head(pos.ordered.genes,1) ## GeneID: 1728
  max.sym <- mget(max.corr.txn,org.Hs.egSYMBOL) ## SYMBOL: NQO1
  max.corr <- cor(txn[max.corr.txn,],gi.clust[,fin.clust],method="spearman")
  
  min.corr.txn <- tail(pos.ordered.genes,1) ## GeneID: 84142
  min.sym <- mget(min.corr.txn,org.Hs.egSYMBOL) ## SYMBOL: FAM175A
  min.corr <-  cor(txn[min.corr.txn,],gi.clust[,fin.clust],method="spearman")

  pdf("EDfigure_9a_max_min_genes.pdf",width=4,height=4)
  par(mar=c(5,5,3,3)+.1)
  plot(txn[max.corr.txn,],gi.clust[,fin.clust],pch=20,
       xlab=paste(max.sym,"mRNA"),ylab="logGI50",main="",cex=.7)
  coef<- summary(lmrob(gi.clust[,fin.clust]~txn[max.corr.txn,]))$coef[1:2]
  abline(coef,col=2)
  plot(txn[min.corr.txn,],gi.clust[,fin.clust],pch=20,
       xlab=paste(min.sym,"mRNA"),ylab="logGI50",cex=.7)
  coef<- summary(lmrob(gi.clust[,fin.clust]~txn[min.corr.txn,]))$coef[1:2]
  abline(coef,col=2)
  dev.off()
}


##[1st-criterion: significance] Computing the significance of the enrichment score with 1000 permutations
system.time(es.pvals <- mclapply(GOs.ge10,function(go){ ## 418 secs
  pvals.ES.GSEA(pos.ordered.genes,go,nperm=1000) ## enrichment score and P-values for each GO term
},mc.cores=(n.cores-1))) 
slogps <- sapply(es.pvals,function(x){ ## signed log P-values
  p <- x$p.val
  slogp <- -log10(p)*sign(x$ES)
  return(slogp)
})

##[2nd-criterion: exclusiveness] Computing enrichment scores of each GO term against 53 clusters
grps <- colnames(corrs)
es.grps <- sapply(grps,function(grp){
  pos <- names(sort(corrs[,grp],decreasing=T))
  ES.grp <- sapply(GOs.ge10,function(go){
    set.seed(12345)
    GSEA.EnrichmentScore3(pos,go)
  })
})

fi.gos <- apply(es.grps,1,function(x){
  x <- unlist(x)
  o <- order(x,decreasing=T)
  if (o[1]==fin.clust){
    return("pos")
  }else if(o[length(o)]==fin.clust){
    return("neg")
  }else{
    return("NA")
  }
})

## summary: GOs satisfying both criteria.
pos.corr.gos <- names(which(fi.gos=="pos"&slogps==Inf)) ## 13 GOs (15??)
neg.corr.gos <- names(which(fi.gos=="neg"&slogps==-Inf)) ## 4 GOs
##
sig.gos <- c(pos.corr.gos,neg.corr.gos)
names(sig.gos) <- c("GO:0015674","GO:0006732","GO:0006807","GO:0030001","GO:0009308","GO:0006812",
                    "GO:0015672","GO:0016616","GO:0009055","GO:0016491","GO:0016651","GO:0005216",
                    "GO:0016614","GO:0007093","GO:0002376","GO:0019221","GO:0003729")
names(pos.corr.gos) <- names(sig.gos[1:13])
names(neg.corr.gos) <- names(sig.gos[14:17])

####################################################################################################
## Analysis of Overlapped genes between GOs ########################################################
####################################################################################################
## analysis of overlapped genes between GO Terms: discover hierarchical structures among GO Terms
pos.overlap <- sapply(pos.corr.gos,function(i){
  sapply(pos.corr.gos,function(j){
    sum(GOs.ge10[[i]] %in% GOs.ge10[[j]])/min(sapply(GOs.ge10[c(i,j)],length))
  })
})
pos.goids <- names(pos.corr.gos)
colnames(pos.overlap) <- rownames(pos.overlap) <- pos.goids
##
neg.overlap <- sapply(neg.corr.gos,function(i){
  sapply(neg.corr.gos,function(j){
    sum(GOs.ge10[[i]] %in% GOs.ge10[[j]])/min(sapply(GOs.ge10[c(i,j)],length))
  })
})
neg.goids <- names(neg.corr.gos)
colnames(neg.overlap) <- rownames(neg.overlap) <- neg.goids

##
if(0){
  setwd(fig.dir)
  pdf("overlap-between-GOs.pdf",width=5,height=5)
  cols <- colorRampPalette(c("white","darkblue"))(10)
  n1 <- heatmap.2(pos.overlap,col=cols,trace="none",
                  margin=c(10,10),cexRow=.7,cexCol=.7,
                  hclustfun=hclustfun,
                  distfun=distfun) ## overlap of genes between positively correlated GOs.
  n2 <- heatmap.2(neg.overlap,col=cols,trace="none",
                  margin=c(10,10),cexRow=.7,cexCol=.7,
                  hclustfun=hclustfun,
                  distfun=distfun) ## overlap of genes between negatively correlated GOs.
  dev.off()
}
if(0){
  ## Extended Data Fig.9c Significant and Exclusive GO Terms against the Ferroptosis Cluster.
  setwd(fig.dir)
  ordered.sig.gos <- c(pos.corr.gos[rev(n1$rowInd)],neg.corr.gos[rev(n2$rowInd)])
  write.csv(cbind(Name=ordered.sig.gos,
                  Correlation=rep(c("positive","negative"),c(13,4)),
                  GO.ID=names(ordered.sig.gos),
                  N.Genes=sapply(GOs.ge10[ordered.sig.gos],length),
                  GO.Class=c(1,1,1,1,1,2,3,3,4,4,4,4,4,5:8)),## GO.Class was from overlap analysis
            "EDFigure_9c.csv",row.names=FALSE)
}

##
if(0){
  ## Extended Data Figure 9b: ferroptosis-associated GOs against all the clusters.
  setwd(dir)
  source("EDFig9b_compute_pvalues.r")
}

## NAD(P)(H)-dependent gene set (after Figure 4d)
nadph.genes <- unlist(GOs.ge10[ c("GO_OXIDOREDUCTASE_ACTIVITY_GO_0016616",
                                  "GO_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NADH_OR_NADPH")])

## Running enrichment scores from GSEA
if(0){
  r1 <- GSEA.EnrichmentScore(pos.ordered.genes,nadph.genes,weighted.score.type=0,correl.vector=NULL)
  ## Extended Data Fig. 9d: Visualizing Enrichment of NAD(P)(H)-dependent gene set
  setwd(fig.dir)
  pdf("EDFigure_9d.pdf",width=5,height=5)
  plot(c(-0.08,r1$RES),type="n",xlab="",ylab="",axes=F,
       main="Gene Set (GO:0016616 + GO:0016651)")
  lines(rep(r1$arg.ES,2),c(0,r1$ES),col="grey60",lty=2)
  lines(c(par()$usr[1],r1$arg.ES),rep(r1$ES,2),col="grey60",lty=2)
  lines(c(0,r1$RES),col=3,lwd=2)
  mtext("Running ES",2,line=3)
  mtext("Gene Index",1,line=3)
  ind <- which(r1$indicator==1)+1
  x <- sapply(seq(ind),function(i)lines(rep(ind[i],2),c(0,-.04),col=2))
  cor.fi <- sort(corrs[pos.ordered.genes,fin.clust],decreasing=T)
  rn <- range(cor.fi,na.rm=T)
  rng <- round(rn/max(abs(rn))*256)
  cols <- colorRampPalette(c("blue","white","red"))(513)[match(seq(rng[1],rng[2]),-256:256)]
  scale <- seq(rn[1],rn[2],length.out=length(cols))
  x <- sapply(seq(r1$RES),function(i){
    j <- which.min(abs(scale-cor.fi[i]))
    lines(rep(i,2)+1,c(-.04,-.08),col=cols[j])
  })
  abline(h=c(0,-.04))
  box()
  axis(2,at=seq(0,.4,.1),las=1)
  axis(1,at=c(0,5000,10000))
  dev.off()
  
  ## Extended Data Fig. 9e
  border <- r1$arg.ES
  ranks <- which(r1$indicator==1)
  nadph.eids <- pos.ordered.genes[ranks]
  nadph.symbols <- unlist(mget(nadph.eids,org.Hs.egSYMBOL))
  nadph.tab <- cbind(rank=ranks,
                     gene.symbol=nadph.symbols,
                     entrez.id=nadph.eids)
  ind.16 <- which(nadph.eids %in% GOs.ge10[["GO_OXIDOREDUCTASE_ACTIVITY_GO_0016616"]])
  ind.51 <- which(nadph.eids %in% GOs.ge10[["GO_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NADH_OR_NADPH"]])
  write.csv(nadph.tab[ind.16,],"EDFigure_9e_genes_go_0016616.csv",row.names=F)
  write.csv(nadph.tab[ind.51,],"EDFIgure_9e_genes_go_0016651.csv",row.names=F)
}


es.nadph <- sapply(grps,function(grp){ ## enrichment score of NAD(P)(H)-dependent genes against 53 clusters.
  pos <- names(sort(corrs[used.genes.ge10,grp],decreasing=T))
  set.seed(12345)
  ES.grp <- GSEA.EnrichmentScore3(pos,nadph.genes)
})
     
## p-values of NAD(P)(H)-associated redox genes
system.time(p0 <- sapply(grps,function(grp){ # 63 secs
  pos <- names(sort(corrs[used.genes.ge10,grp],decreasing=T))
  pvals.ES.GSEA(pos,nadph.genes,nperm=1000)$p.value
}))
##
ind.p0 <- which(p0<0.01)
system.time(p1 <- sapply(ind.p0,function(i){ # 151 secs (53 secs)
  pos <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
  pvals.ES.GSEA(pos,nadph.genes,nperm=1e4)$p.value
}))
##
ind.p1 <- names(p1)[which(p1<1e-3)] ## only "1"
if(0){
  system.time({ ## takes 4.6 hours with 7 cores
    nperm <- 1e8
    pos <- names(sort(corrs[used.genes.ge10,"1"],decreasing=T))
    original.ES <- GSEA.EnrichmentScore3(pos,nadph.genes)
    perm.ES <- unlist(mclapply(seq(nperm),function(i){
      set.seed(i)
      GSEA.EnrichmentScore3(sample(pos),nadph.genes)
    },mc.cores=7,mc.set.seed=FALSE))
    p.value <- min(sum(perm.ES > original.ES),sum(perm.ES < original.ES))/nperm
    p2 <- p.value ## p2 <- 3e-08
  })
}else{
  p2 <- c(`1`=3e-08)
}

## Figure 4e: NAD(P)(H)-dependent genes' exclusive correlation against the ferroptosis cluster
pvalues <- p0
pvalues[names(p1)] <- p1
pvalues[names(p2)] <- p2
slogp.nadph <- -log10(pvalues)*sign(es.nadph)
##
if(0){
  setwd(fig.dir)
  pdf("Figure_4e.pdf",width=1.8,height=4)
  set.seed(1)
  plot(jitter(rep(1,53),amount=.2),slogp.nadph,pch=20,xlim=c(0,2),type="n",axes=F,
       xlab="",ylab="",ylim=c(-2,8))
  abline(h=0)
  abline(v=1,col="grey80")
  axis(2,at=seq(-2,8,2))
  mtext("P-value",2,3)
  points(jitter(rep(1,53),amount=.3),slogp.nadph,pch=20,col=(seq(53) %in% fin.clust)+1,cex=.5)
  dev.off()
}


##
sessionInfo()
##
## R version 2.15.2 (2012-10-26)
## Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] org.Hs.eg.db_2.8.0   RSQLite_0.11.2       DBI_0.2-5           
##  [4] AnnotationDbi_1.20.7 Biobase_2.18.0       BiocGenerics_0.4.0  
##  [7] robustbase_0.9-7     mclust_4.0           gplots_2.11.0       
## [10] MASS_7.3-23          KernSmooth_2.23-10   caTools_1.14        
## [13] gdata_2.12.0         gtools_2.7.1         RColorBrewer_1.0-5  
## [16] kslib_2.2           
## 
## loaded via a namespace (and not attached):
## [1] IRanges_1.16.6  bitops_1.0-4.2  compiler_2.15.2 stats4_2.15.2  
## [5] tools_2.15.2   
##
