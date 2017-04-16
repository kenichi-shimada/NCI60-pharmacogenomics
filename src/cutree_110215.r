xx <- myheatmap(gi.clust,cexRow=.5,cexCol=.5,keysize=1.5)
xhc <- hclustfun(distfun(t(gi.clust)))

cut <- sapply(1:53,function(x){
    ct <- cutree(xhc,k=x)
    nf <- ct["2"]
    sum(ct==nf)
});names(cut) <- 1:53

## k=18 - let's go with it!
ct18 <- cutree(xhc,k=18)
m18 <- lapply(seq(18),function(x)names(ct18[ct18==x]))
names18 <- sapply(m18,function(x){
    if(length(x)==1)return(x)
    if(any(x=="2"))return("g1")
    if(any(x=="17"))return("g2")
    if(any(x=="8"))return("g3")
    if(any(x=="9"))return("g4")
    if(any(x=="5"))return("g5")
    if(any(x=="7"))return("g6")
})
                                 
## 
names(m18) <- names18##paste("n",1:18,sep="")
n18 <- rep(names(m18),sapply(m18,length));names(n18) <- unlist(m18)
new.col <- rev(unique(n18[as.character(xx$rowInd)]))
new.row <- rownames(gi.clust)[xx$colInd]
##
gi.clust18 <- sapply(m18,function(l){
    if(length(l)==1)return(gi.clust[,l])
    if(length(l)>1)return(apply(gi.clust[,l],1,mean))
})[new.row,]

##
n.clust18 <- sapply(m18,function(m)sum(as.character(mem) %in% m)) ## gi50 profiles of the 53 clusters

##
if(0){
    ## for making Fig. 4b
    cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(257)
    pdf("figure_2_heatmap18.pdf",width=7,height=5) ##heatmap
    yy <- heatmap.2(t(gi.clust18),col=cols,trace="none",
                    cexRow=.5,cexCol=.5,keysize=1.5,
                    Colv=FALSE, Rowv=TRUE,dendrogram="row")
    dev.off()
    ##
    pdf("figure_2_barplot18.pdf",width=2,height=5)
    barplot((n.clust18[yy$rowInd]),horiz=T,cex.names=.3,border=NA)##barplot
    dev.off()
}
