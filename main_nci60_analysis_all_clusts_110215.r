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
## setwd(dir <- "/root/of/github/repo/")

## working directory
fig.dir <- paste(dir,"/figs",sep="") ## where you save figures under wd

## Load libraries -------------------------------------------------------------------------
library(gplots) ## heatmap.2 function
library(RColorBrewer) ## gradation of colors in heatmap
library(mclust) ## model-based clustering
library(robustbase) ## lmrob for robust linear regression
library(org.Hs.eg.db) ## Gene ID <=> Gene SYMBOL
library(parallel) ## mclapply to use multiple cores 

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
system.time({ ## 222 secs/492secs
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
system.time({## 244 secs
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

## more summary stats
if(0){
    ## iqrs
    pdf("figure_4d_boxplot_iqr.pdf",width=2,height=5) ##fig.4d
    sele.iqrs <- apply(sele.gi,1,IQR)
    clust.iqrs <- sapply(seq(g),function(m)sele.iqrs[mem==m]) ## gi50 profiles of the 53 clusters
    names(clust.iqrs) <- seq(g)
    boxplot(clust.iqrs[as.character(xx$rowInd)],las=2,col="lightblue",horizontal=T,cex=.1,pch=20,cex.axis=.3,lwd=.5)
    dev.off()

    ## iqrs - n.cmpds
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

## 53 -> 18 clusters
source("src/cutree_110215.r")

#######################################################################################
## Transcriptome analysis of the ferroptosis cluster ##################################
#######################################################################################
## transcriptome data are available at 
## https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data
## used GENELOGIC's U133 array data (RMA-normalized)
## Order of cell lines should be the sames between GI50 matrix and transcription matrix

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
save(corrs,file="rda/corrs.rda")

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

## for all classes
grps <- colnames(corrs)
system.time(es.grps <- sapply(grps,function(grp){
    pos <- names(sort(corrs[,grp],decreasing=T))
    ES.grp <- sapply(GOs.ge10,function(go){
        set.seed(12345)
        GSEA.EnrichmentScore3(pos,go)
    })
})) ## 16.8 secs

##
## 1/21/22 checked up to here.
##

sig.gos.all <- lapply(grps,function(this.clust){
    pos.ordered.genes <- names(sort(corrs[used.genes.ge10,this.clust],decreasing=T))
    ##[1st-criterion: significance] Computing the significance of the enrichment score with 1000 permutations
    es.pvals <- mclapply(GOs.ge10,function(go){ ## 418 secs
        cat("*")
        pvals.ES.GSEA(pos.ordered.genes,go,nperm=1000) ## enrichment score and P-values for each GO term
    },mc.cores=(n.cores-1))
    slogps <- sapply(es.pvals,function(x){ ## signed log P-values
        p <- x$p.val
        slogp <- -log10(p)*sign(x$ES)
        return(slogp)
    })

    ##[2nd-criterion: exclusiveness] Computing enrichment scores of each GO term against 53 clusters
    this.gos <- apply(es.grps,1,function(x){
        x <- unlist(x)
        s <- sort(x,decreasing=T)
        if (names(s)[1]==this.clust){
            return("pos")
        }else if(names(s)[length(s)]==this.clust){
            return("neg")
        }else{
            return("NA")
        }
    })
    
    ## summary: GOs satisfying both criteria.
    pos.corr.gos <- names(which(this.gos=="pos"&slogps==Inf)) ## 13 GOs (15??)
    neg.corr.gos <- names(which(this.gos=="neg"&slogps==-Inf)) ## 4 GOs

    cat("*")
    return(list(pos=pos.corr.gos,neg=neg.corr.gos))
})
save(sig.gos.all,file="rda/sig.gos.all-18.rda")

l2 <- sapply(sig.gos.all,function(x)sapply(x,length))
colnames(l2) <- grps

##
##source("~/Projects/src/R/hpc_temp/hpc_TEA/TEA.hpc_4.1.R")

##
if(0){
    used.db <- GOs.ge10
    if(0){
      jns1 <- sapply(1:18,function(i){
        sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
        ##
        jn <- TEA.hpc(sorted.e,used.db,nperm=10^4,wt="24:00:00",mem="2gb",
                      n.split=256,jobname=paste("clust18-",i,sep=""),verbose=T)
        return(jn)
      })
      save(jns1,file="jns_020115_5.rda")#10^5
    }else{
      load("jns_020115_5.rda");#jns1
    }
    
    ##
    slp18 <- sapply(jns1,function(jn){
        cat("*")
        outs <- call.TEA.out(jn)
        slps <- -log10(outs$p)* sign(outs$ES)
        return(slps[names(GOs.ge10)])
    })
    colnames(slp18) <- colnames(corrs)
    slp18[slp18 == Inf] <- 5
    slp18[slp18 == -Inf] <- -5
    cols <- colorRampPalette(c("cyan","black","yellow2"))(257)
    heatmap.2(slp18,trace="none",col=cols)
    
    ##
    ## 10^4 -> slps gt 3 or lt -3 -> 10^5
    totest5 <- lapply(1:18,function(i){
        x <- slp18[,i]
        which(x > 3 | x < -3)
    })
    jns5 <- sapply(1:18,function(i){
        sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
        used.db <- GOs.ge10[totest5[[i]]]
        ##
        n.split <- min(256,length(used.db))
        jn <- TEA.hpc(sorted.e,used.db,nperm=10^5,wt="24:00:00",mem="2gb",
                      n.split=n.split,jobname=paste("clust18-5-",i,sep=""),verbose=T)
        return(jn)
    })
    slp18.1 <- sapply(1:18,function(i){
        cat("*")
        jn <- jns5[i]
        outs <- call.TEA.out(jn)
        slps <- -log10(outs$p)* sign(outs$ES)
        this <- slp18[,i]
        this[totest5[[i]]] <- slps[names(GOs.ge10[totest5[[i]]])]
        return(this)
    })
    save(jns5,file="jns_020115_10.5.rda")#10^5
    ## 10^5 - slp gt 4 or lt -3 -> 10^6
    totest6 <- lapply(1:18,function(i){
        x <- slp18.1[,i]
        which(x > 4 | x < -4)
    })
    jns6 <- sapply(1:18,function(i){
        sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
        used.db <- GOs.ge10[totest6[[i]]]
        ##
        n.split <- min(256,length(used.db))
        jn <- TEA.hpc(sorted.e,used.db,nperm=10^6,wt="24:00:00",mem="2gb",
                      n.split=n.split,jobname=paste("clust18-6-",i,sep=""),verbose=T)
        return(jn)
    }) 
    save(jns6,file="jns_020115_10.6.rda")#10^6 
    
    ##
    for(i in 1:18){
        cat("*")
        jn <- jns6[i]
        outs <- call.TEA.out(jn)
    }

    slp18.2 <- sapply(1:18,function(i){
        cat("*")
        jn <- jns6[i]
        this <- slp18.1[,i]
        ##
        outs <- try(call.TEA.out(jn),silent=T)
        if(class(outs)!="try-error"){
            slps <- -log10(outs$p)* sign(outs$ES)
            this[totest6[[i]]] <- slps[names(GOs.ge10[totest6[[i]]])]
        }
        return(this)
    })
    
    ##
    ## 10^6 - slp gt 5 or lt -5 -> 10^7
    totest7 <- lapply(1:18,function(i){
        x <- slp18.2[,i]
        which(x > 5 | x < -5)
    })
    jns7 <- sapply(1:18,function(i){
        sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
    used.db <- GOs.ge10[totest7[[i]]]
        ##
        n.split <- min(256,length(used.db))
        jn <- TEA.hpc(sorted.e,used.db,nperm=10^7,wt="24:00:00",mem="2gb",
                      n.split=n.split,jobname=paste("clust18-7-",i,sep=""),verbose=T)
        return(jn)
    })
    save(jns7,file="jns_020115_10.7.rda")#10^6 ## waiting here (2/1/15)
    setwd(dir <- "~/Projects/TF_TEA/505_sd7/nci60_analysis/")
    if(0)save.image("020215-upto-jns7.rda")
}
    
setwd(dir <- "/Users/kenichi/Desktop/pharmacogenomics/nci60_analysis")
load("020215-upto-jns7.rda")

##################################################
if(0){
  slp18.3 <- sapply(1:18,function(i){
    cat("*")
    jn <- jns7[i]
    this <- slp18.2[,i]
    ##
    outs <- try(call.TEA.out(jn),silent=T)
    if(class(outs)!="try-error"){
      slps <- -log10(outs$p)* sign(outs$ES)
      this[totest7[[i]]] <- slps[names(GOs.ge10[totest7[[i]]])]
    }
    return(this)
  })
  save(slp18.3,file="slp18_3.rda")
}else{
  load("slp18_3.rda")
}

if(0){
  ## 10^7 - slp gt 6 or lt -6 -> 10^8
  totest8 <- lapply(1:18,function(i){
    x <- slp18.3[,i]
    which(x > 6 | x < -6)
  })
  jns8 <- sapply(1:18,function(i){
    sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
    used.db <- GOs.ge10[totest8[[i]]]
    ##
    n.split <- min(256,length(used.db))
    jn <- TEA.hpc(sorted.e,used.db,nperm=10^8,wt="24:00:00",mem="2gb",
                  n.split=n.split,jobname=paste("clust18-8-",i,sep=""),verbose=T)
    return(jn)
  })
  ################################################
  setwd(dir <- "~/Projects/TF_TEA/505_sd7/nci60_analysis/")
  save(jns8,file="jns_020315_10.7.rda")#10^6 ## waiting here (2/1/15)
  if(0)save.image("020315-upto-jns8.rda")
  if(0)load("020215-upto-jns7.rda")
  ################################################
  
  slp18.4 <- sapply(1:18,function(i){
    cat("*")
    jn <- jns8[i]
    outs <- call.TEA.out(jn)
    slps <- -log10(outs$p)* sign(outs$ES)
    this <- slp18.3[,i]
    this[totest8[[i]]] <- slps[names(GOs.ge10[totest8[[i]]])]
    return(this)
  })
}

## two criteria:
ns <- names(GOs.ge10)
own <- 3;others <- 3
sigs <- lapply(1:18,function(g){
  ## max
  m1 <- apply(slp18.3,1,function(x){
    i1 <- which.max(x);
    if(i1==g){
      tf <- x[i1]> own && max(x[-i1]) < others
      return(tf)
    }else{
      return(FALSE)
    }
  });maxi <- which(m1)
  m2 <- apply(slp18.3,1,function(x){
    i1 <- which.min(x);
    if(i1==g){
      tf <- x[i1]< -own && max(x[-i1]) > -others
      return(tf)
    }else{
      return(FALSE)
    }
  });mini <- which(m2)
  return(list(pos=ns[m1],neg=ns[m2]))
})
names(sigs) <- colnames(gi.clust18)
sig.names <- sub("^GO_","",unique(unlist(sigs)));names(sig.names) <- c()
setwd(dir <- "~/Projects/TF_TEA/505_sd7/nci60_analysis/")
save(sig.names,file="go.names.rda")

##
m <- sapply(sigs,function(x)length(unlist(x)))

if(0){
    sig.gos <- c(pos.corr.gos,neg.corr.gos)
    names(sig.gos) <- c("GO:0015674","GO:0006732","GO:0006807","GO:0030001","GO:0009308","GO:0006812",
                        "GO:0015672","GO:0016616","GO:0009055","GO:0016491","GO:0016651","GO:0005216",
                        "GO:0016614","GO:0007093","GO:0002376","GO:0019221","GO:0003729")
    names(pos.corr.gos) <- names(sig.gos[1:13])
    names(neg.corr.gos) <- names(sig.gos[14:17])
}



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

## NADPH
if(0){
  jns.nadph <- sapply(1:18,function(i){
    sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
    used.db <- list(redox=nadph.genes)
    ##
    n.split <- min(256,length(used.db))
    jn <- TEA.hpc(sorted.e,used.db,nperm=10^7,wt="24:00:00",mem="2gb",
                  n.split=n.split,jobname=paste("clust18-7-nadph-",i,sep=""),verbose=T)
    return(jn)
  })
  save(jns.nadph,file="jns_nadph.rda")#10^6 ## waiting here (2/1/15)
}else{
  load("jns_nadph.rda")
}
##
if(0){
  jns.nadph.1 <- sapply(1:18,function(i){
    sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
    used.db <- list(redox=nadph.genes)
    ##
    n.split <- min(256,length(used.db))
    jn <- TEA.hpc(sorted.e,used.db,nperm=10^7,wt="24:00:00",mem="2gb",
                  n.split=n.split,jobname=paste("clust18-7-nadph-",i,sep=""),verbose=T)
    return(jn)
  })
  save(jns.nadph,file="jns_nadph.rda")#10^6 ## waiting here (2/1/15)
}else{
  load("jns_nadph.rda")
}

slp.nadph <- sapply(1:18,function(i){
  cat("*")
  jn <- jns.nadph[i]
  ##
  outs <- try(call.TEA.out(jn),silent=T)
  if(class(outs)!="try-error"){
    slps <- -log10(outs$p)* sign(outs$ES)
  }
  return(slps)
})
## cmp.clust=1, 10^9 iterations
sorted.e <- names(sort(corrs[used.genes.ge10,1],decreasing=T))
used.db <- list(redox=nadph.genes)
##
n.split <- 1000##min(256,length(used.db))
jn <- TEA.hpc2(sorted.e,used.db,nperm=10^6,wt="24:00:00",mem="2gb",
               n.split=n.split,jobname=paste("nadph-6-",i,sep=""),verbose=T)
outs <- try(call.TEA.out(jn),silent=T)
  if(class(outs)!="try-error"){
    slps <- -log10(outs$p)* sign(outs$ES)
  }
  return(slps)
})
p1 <- outs$pval
slp.nadph[1] <- -log10(sum(head(outs)$pvals)/length(p1))*sign(outs$ES[1])
##
library(qvalue)
qvalue(0.1^abs(slp.nadph) ) : not working...
##
setwd(dir);setwd("figs")
pdf("slp.nadph.pdf",width=3,height=5)
set.seed(1234);xs <- jitter(rep(1,18),amount=.15)
ys <- slp.nadph
plot(xs,ys,pch=20,xlim=c(0,2),axes=F,type="n",ylim=c(-2,8),
     ylab="-log10 p-value")
abline(v=1,col="grey80")
abline(h=0)
points(xs,ys,pch=20)
axis(2,at=seq(-2,8,2))
dev.off()

## Running enrichment scores from GSEA

##
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
getwd()
##
