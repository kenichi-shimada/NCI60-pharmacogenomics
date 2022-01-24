xx <- myheatmap(used.gi.fin,cexRow=.5,cexCol=.5,keysize=1.5)
mat <- used.gi.fin[c("erastin","AE","PE","MEII","RSL3","1S3R","ML162","ML210","FIN56","CIL56"),]

if(0){
	## fig. 3a
	cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(513)
	setwd("/Users/ks2474/Projects/TF_TEA/505_sd7/nci60_analysis/figs")
	pdf("fin.gi.pdf")
	heatmap.2(t(mat),col=cols,trace="none",Colv=FALSE,dendrogram="row",density.info="none")
	dev.off()
}
##
if(0){
	## barplot of GI50 in 59 cell lines (not used in the paper)
	cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(513)
	fin56 <- sort(used.gi.fin["FIN56",])
	max.abs <- max(abs(range(used.gi.fin)))
	used.cols <- cols[match(round(fin56/max.abs*255),-255:255)]
	pdf("fin56.barplot.pdf")
	barplot(fin56,col=used.cols,border=NA,las=2,xlab="")
	abline(h=0)
	abline(h=quantile(fin56,c(.25,.75)),lty=2)
	dev.off()
}