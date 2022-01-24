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
	## not used
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
	## not used
	setwd(fig.dir)
	ordered.sig.gos <- c(pos.corr.gos[rev(n1$rowInd)],neg.corr.gos[rev(n2$rowInd)])
	write.csv(cbind(Name=ordered.sig.gos,
	              Correlation=rep(c("positive","negative"),c(13,4)),
	              GO.ID=names(ordered.sig.gos),
	              N.Genes=sapply(GOs.ge10[ordered.sig.gos],length),
	              GO.Class=c(1,1,1,1,1,2,3,3,4,4,4,4,4,5:8)),## GO.Class was from overlap analysis
	        "EDFigure_9c.csv",row.names=FALSE)
}

