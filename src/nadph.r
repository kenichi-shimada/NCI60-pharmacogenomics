## NAD(P)(H)-dependent gene set (after Figure 4d)
nadph.genes <- unlist(GOs.ge10[ c("GO_OXIDOREDUCTASE_ACTIVITY_GO_0016616",
                                  "GO_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NADH_OR_NADPH")])

## NADPH
jns.nadph <- sapply(1:18,function(i){
	sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
	used.db <- list(redox=nadph.genes)
	##
	n.split <- min(256,length(used.db))
	jn <- TEA.hpc(sorted.e,used.db,nperm=10^7,wt="24:00:00",mem="2gb",
	              n.split=n.split,jobname=paste("clust18-7-nadph-",i,sep=""),verbose=T)
	return(jn)
})

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
slps <- -log10(outs$p)* sign(outs$ES)
p1 <- outs$pval

slp.nadph[1] <- -log10(sum(head(outs)$pvals)/length(p1))*sign(outs$ES[1])

## Fig. 6D
if(0){
	setwd(dir);setwd("figs")
	pdf("slp.nadph.pdf",width=3,height=5)
	set.seed(1234)
	xs <- jitter(rep(1,18),amount=.15)
	ys <- slp.nadph
	plot(xs,ys,pch=20,xlim=c(0,2),axes=F,type="n",ylim=c(-2,8),
	     ylab="-log10 p-value")
	abline(v=1,col="grey80")
	abline(h=0)
	points(xs,ys,pch=20)
	axis(2,at=seq(-2,8,2))
	dev.off()
}

## Running enrichment scores from GSEA

##
if(0){
  r1 <- GSEA.EnrichmentScore(pos.ordered.genes,nadph.genes,weighted.score.type=0,correl.vector=NULL)
  ## not used
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
