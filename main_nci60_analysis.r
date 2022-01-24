#=================================================================================
# Title:        Model-based clustering of GI50 profiles of NCI60 cell lines
# Author:       Kenichi Shimada
# Date:         18 Sep 2014
# R-version:    2.15.2 (2015-02-07)
#=================================================================================

## Note on 01/21/2022:
## 
## This script should be run inside the root directory of the nci60_analysis.
## If redo the analysis after downloading the github repo
## (https://github.com/kenichi-shimada/NCI60-pharmacogenomics), set the root dir of the repo to 'dir':
##
## dir <- "/root/of/github/repo/ # specified in .Rprofile
setwd(dir)

## working directory
fig.dir <- paste(dir,"/figs",sep="") ## where you save figures under wd

## Load libraries -------------------------------------------------------------------------
library(gplots) ## heatmap.2 function
library(RColorBrewer) ## gradation of colors in heatmap
library(mclust) ## model-based clustering
library(robustbase) ## lmrob for robust linear regression
library(org.Hs.eg.db) ## Gene ID <=> Gene SYMBOL
library(parallel) ## mclapply to use multiple cores 
library(qvalue)

## number of cores available in the machine
n.cores <- parallel:::detectCores() 

## load r functions to use in this analysis
source("src/rutils.r")

## Load objects & trimming data
x <- load("rda/nci60-analysis.rda")

## 'nci60-analysis.rda' contains the following R objects:

## gi: logGI50 values of 75446 (50839 unique) compounds in 59 cell lines
## gi.fin: logGI50 values of 21 (10 unique) ferroptosis inducers in 59 cell lines
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
system.time({ ## 172 secs (run on 1/23/22)
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
save(med.gi,file="rda/med.gi.rda")

if(0){
  load("rda/med.gi.rda") # the rda file is available on github
}

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

## impute a missing value
used.gi.fin["MEII","LC:EKVX"] <- mean(used.gi.fin[c("AE","PE","erastin"),"LC:EKVX"])

source("src/fig3_fin_nci60.r")
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
  mgi50 <-  list(used.iqrs,mir)
  a <- boxplot(mgi50,border=NA,cex=0,axes=F)
  axis(2)
  axis(1,at=1:2,labels=c("Others","Ferroptosis"),las=2)
  set.seed(543)
  boxplot(mgi50,cex=0,axes=F,add=T,border="grey40")
  points(jitter(rep(1,length(mgi50[[1]])),amount=.35),mgi50[[1]],pch=20,cex=.15,col="grey20")
  cols <- brewer.pal(6,"Set1")
  set.seed(1245)
  abline(h=min(mgi50[[2]][-1]),col=2,lty=2)
  points(jitter(rep(2,length(mgi50[[2]])),amount=.2),mgi50[[2]],pch=20,cex=.7,col="grey20")
  dev.off()
}

## Selecting cell line-selective compounds
has.large.iqr <- used.iqrs > min(used.iqrs.fin[-1]) ## 2555 compounds
pval <- t.test(used.iqrs, used.iqrs.fin)$p.value ## 0.00869
sele.gi <- rbind(used.gi[has.large.iqr,],used.gi.fin) ## GI50 values of 2565 cell-line-selective lethals
####################################################################################################
## Model-based clustering of GI50 profiles of cell line-selectie comopunds #########################
####################################################################################################
## find the optimal model and # of clusters
n.clusts.range <-  20:80
system.time({## 209 secs (1/23/22)
  bics <- mclustBIC(sele.gi,G=n.clusts.range,modelNames=c("EII","VII"),
                    prior=priorControl())
})
max.ind <- which(bics==max(bics),arr.ind=T)
model <- colnames(bics)[max.ind[2]] ## optimal model: VII
g <- as.numeric(rownames(bics)[max.ind[1]]) ## optimal # of clusters: 53 - not always the same?
## It also looks like mclust algorithm changed. g = 39 when run on 1/21/2022

if(0){
  ## Figure choose the optimal number of clusters
  pdf("choose_n_clusters_with_BIC.pdf")
  plot(n.clusts.range,bics[,model],type="b",
       xlab="# of clusters",ylab="Bayesian Information Criterion",
       main=paste("model:",model))
  abline(v=g,col=2)
  dev.off()
}

##
system.time(mc <- Mclust(sele.gi,G=g,modelNames=model,prior=priorControl())) ## 47 secs
save(mc,file="rda/mc.rda")

## Discover the ferroptosis cluster
mem <- mc$classification ## compounds' memberships
table(mem)
mem.fins <- mem[uniq.fin]
fin.clust <- names(sort(table(mem.fins),decreasing=T))[1] ## ferroptosis cluster is "1"-st cluster

## GI50 profile of each cluster
gi.clust <- sapply(seq(g),function(m)apply(sele.gi[mem==m,,drop=F],2,median)) ## gi50 profiles of the 53 clusters
colnames(gi.clust) <- as.character(seq(g))

if(0){
  ## Fig 4a heatmap
  setwd(fig.dir)
  pdf("figure_4a_heatmap.pdf",width=7,height=5) ##fig.4a
  colors <- c("white","red")[seq(as.numeric(g)) %in% fin.clust+1]
  xx <- myheatmap(gi.clust,cexRow=.5,cexCol=.5,keysize=1.5,RowSideColors=colors)
  dev.off()
  ## Fig 4a barplot (not used in the fig)
  pdf("figure_4a_barplot.pdf",width=2,height=5) ##fig.4a
  tab.fins <- table(mem,rep(0:1,c(2555,10)));rownames(tab.fins) <- as.character(1:53)
  sort.tab.fins <- tab.fins[xx$rowInd,]
  barplot(t(sort.tab.fins[,2:1]),beside=F,las=2,col=c(2,"lightblue"),horiz=T,cex.names=.3,border=NA)
  dev.off()
  ## heatmap, only for fin clusters (not used in the paper)
  pdf("figure_fin_heat.pdf",width=5,height=4) ##fig. fin_heat
  gi50.within.fi.clust <- t(sele.gi[mem==fin.clust,])
  myheatmap(gi50.within.fi.clust[xx$colInd,],cexRow=.7,cexCol=.5,keysize=2,
            Colv=FALSE,dendrogram="row")
  dev.off()
}

## more summary stats (not used)
if(0){
    ## iqrs (not used)
    pdf("figure_4d_boxplot_iqr.pdf",width=2,height=5) ##fig.4d
    sele.iqrs <- apply(sele.gi,1,IQR)
    clust.iqrs <- sapply(seq(g),function(m)sele.iqrs[mem==m]) ## gi50 profiles of the 53 clusters
    names(clust.iqrs) <- seq(g)
    boxplot(clust.iqrs[as.character(xx$rowInd)],las=2,col="lightblue",horizontal=T,cex=.1,pch=20,cex.axis=.3,lwd=.5)
    dev.off()

    ## iqrs - n.cmpds (not used)
    n.cmpds.clust <- sapply(clust.iqrs,length)
    clust.iqrs.med <- sapply(clust.iqrs,median)
    clust.iqrs.mean <- sapply(clust.iqrs,mean)
    clust.iqrs.iqr <- sapply(clust.iqrs,function(x)quantile(x,c(.25,.75)))
    plot(clust.iqrs.med,n.cmpds.clust,pch=20,col=rep(2:1,c(1,52)),type="n",xlim=range(clust.iqrs.iqr))
    x <- sapply(1:53,function(i)lines(clust.iqrs.iqr[,i],rep(n.cmpds.clust[i],2),col="grey80"))
    points(clust.iqrs.med,n.cmpds.clust,pch=20,col=rep(2:1,c(1,52)))
    ## average RMS - n.cmpds
}

##
if(0){
  setwd(dir)
  write.csv(gi,"logGI50_all.csv")
  write.csv(used.gi.fin,"logGI50_FIN.csv")
  write.csv(cbind(mem,sele.gi),"cell-line.selective.csv")
}

#######################################################################################
## Mechanisms of Action of 18 clusters
#######################################################################################
source("src/moa_cmpds.r") # Also, fig.4b

#######################################################################################
## Converging from 53 to 18 clusters
#######################################################################################
source("src/merge_to_18clusts.r") # Also, fig.4b

#######################################################################################
## Transcriptome analysis of the ferroptosis cluster ##################################
#######################################################################################
## transcriptome data are available at 
## https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data
## used GENELOGIC's U133 array data (RMA-normalized)
## Order of cell lines should be the sames between GI50 matrix and transcription matrix

## use subset of genes that are overlapped between MSigDB and the expression matrix
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
## 10850 genes x 53 compound clusters ## 48 secs with n.cores=8
system.time(corrs <- mclapply(rownames(t.gi50),function(cmpd){
  sapply(used.genes.ge10,function(gid){
    cor(txn[gid,],t.gi50[cmpd,],method="spearman",use="complete.obs")
  })
},mc.cores=(n.cores-1)))

corrs <- do.call(cbind,corrs)
colnames(corrs) <- rownames(t.gi50)
rownames(corrs) <- used.genes.ge10
save(corrs,file="rda/corrs.rda")

pos.ordered.genes <- names(sort(corrs[used.genes.ge10,fin.clust],decreasing=T))

if(0){
  ## not used
  ## scatterplots (txn.exppression vs GI50) of the largest and smallest correlation coefs 
  ## histogram of correlation coef with fin.clust
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
  # pos correlated gene
  par(mar=c(5,5,3,3)+.1)
  plot(txn[max.corr.txn,],t.gi50[fin.clust,],pch=20,
       xlab=paste(max.sym,"mRNA"),ylab="logGI50",main="",cex=.7)
  coef<- summary(lmrob(t.gi50[fin.clust,]~txn[max.corr.txn,]))$coef[1:2]
  abline(coef,col=2)
  plot(txn[min.corr.txn,],t.gi50[fin.clust,],pch=20,
       xlab=paste(min.sym,"mRNA"),ylab="logGI50",cex=.7)
  coef<- summary(lmrob(t.gi50[fin.clust,]~txn[min.corr.txn,]))$coef[1:2]
  abline(coef,col=2)
  dev.off()
}

#######################################################################################
## GSEA
#######################################################################################
## On 1/23/22:
## Note that this GSEA is a fairly primitive function.
## This code is maintained for an archival purpose.
## An improved implementation of GSEA, such as fgsea significantly improve the analysis.

source("src/gsea.r")  


####################################################################################################
## Analysis of Overlapped genes between GO Terms ########################################################
####################################################################################################
## Discover hierarchical structures among GO Terms
source("src/go_overlap.r")

## NADPH-related pathways in depth:
source("src/nadph.r")

##
sessionInfo() # originally in 2015
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
##
