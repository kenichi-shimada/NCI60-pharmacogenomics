fig_moa.dir <- paste(dir,"/figs_moa",sep="") ## where you save figures under wd

load("rda/GI50-75446.rda")

## trim data
uniq.gi50.info <- unique(gi50.info[c(1,4)])
moa.all <- c(uniq.gi50.info$Mech,rep("FIN",10))
names(moa.all) <- c(uniq.gi50.info$NSC,tail(names(mem),10))
moa.all[moa.all %in% c("-","","Unknown")] <- "-"
fin.iqr <- rep(c(FALSE,TRUE),c(1,9));names(fin.iqr) <- tail(names(mem),10)
large.iqr <- c(has.large.iqr,fin.iqr)
cmpds.all <- names(large.iqr)
moa.all <- (moa.all[names(moa.all) %in% cmpds.all])
uniq.moas <- unique(moa.all)
n.sele.moa <- rep(0,length(uniq.moas))
names(n.sele.moa) <- uniq.moas

## 
moa.all[moa.all=="Antifols"] <- "Df"

##
tbl <- table(large.iqr,moa.all)[c("TRUE","FALSE"),]
o <- order(colSums(tbl),decreasing=T)
tbl <- tbl[,o]
##
lbls <- names(which(apply(tbl,2,sum)>1))
tbl <- tbl[,lbls]
pvals <- sapply(colnames(tbl)[-1],function(x){
    i <- match(x,colnames(tbl))
    cont <- cbind(tbl[,i],apply(tbl[,-i],1,sum))
    p <- -log10(fisher.test(cont,alternative="greater")$p.val)
    names(p) <- c()
    return(p)
})
tbl <- tbl[,-1]

if(0){
    ## MOA-statistics - Fig. S2
    setwd(fig_moa.dir)
    pdf("moa_stats.pdf")
    par(mar=c(10,5,4,1)+.1)
    x <- barplot(tbl,beside=F,col=c("lightyellow","lightgreen"),las=2,ylim=c(0,60))
    for(i in seq(x)){
        if(pvals[i]>-log10(0.05))text(x[i],sum(tbl[,i]) + 2,"*")
    }
    legend("topright",c("selective","non-selective"),fill=c("lightyellow","lightgreen"))
    dev.off()

    ##
    pdf("moa_iqrs_plots.pdf")
    ##sele.gi <- rbind(used.gi[has.large.iqr,],used.gi.fin)
    iqrs.all <- c(used.iqrs,used.iqrs.fin)
    iqrs.w.moa <- lapply(colnames(tbl),function(moa){
        cmpds <- names(moa.all[moa.all==moa])
        this.iqrs <- iqrs.all[cmpds]
    })
    names(iqrs.w.moa) <- colnames(tbl)
    a <- boxplot(iqrs.w.moa,border="grey80",pch=NA,cex=.7,axes=F) ;box()
    axis(2)
    axis(1,las=2,at=seq(ncol(tbl)),labels=colnames(tbl))
    for(i in seq(ncol(tbl))){
        set.seed(12345)
        n.pts <- length(iqrs.w.moa[[i]])
        points(jitter(rep(i,n.pts),amount=.2),iqrs.w.moa[[i]],
               col=(names(iqrs.w.moa[[i]]) %in% rownames(sele.gi)[-2556])+1,pch=20,cex=.7)
    }
    abline(h=min(iqrs.w.moa[["FIN"]][-1]),lty=2)
    dev.off()
}

##
cmpds <- names(which(large.iqr))
large.iqr["CIL56"] <- TRUE
sele.moa <- c(moa.all[match(cmpds,names(moa.all))])
t.sele.moa <- table(sele.moa)
n.sele.moa[names(t.sele.moa)] <- t.sele.moa
if(length(has.large.iqr)==6249)has.large.iqr <- c(has.large.iqr,fin.iqr)
logp <- sapply(colnames(tbl),function(m){
    d2 <- table(moa.all==m,has.large.iqr)[c("TRUE","FALSE"),c("TRUE","FALSE")]
    ftest <- fisher.test(d2,alternative="greater")$p.value
    -log10(ftest)
})[colnames(tbl)]


##
if(0){
    ## Fig S4 A-B
    setwd(fig_moa.dir)
    pdf("1d_moa_barplot_1.pdf",width=2,height=5)
    for(moa in used.moa){
        i <- match(moa,used.moa)
        if(all(sele.moa!=moa))next;
        moa.table <- table(sele.moa==moa,mem)[c("TRUE","FALSE"),]
        sort.moa <- moa.table[1,xx$rowInd]
        barplot(sort.moa,beside=F,las=2,col=c(brewer.pal(12,"Paired")[-11][i]),
                horiz=T,cex.names=.3,border=NA,main=moa,xlim=c(0,25))
    }
    dev.off()

    ##
    pdf("1e_moa_all_barplot.pdf",width=2,height=5)
    used.moa <- names(sort(tbl[1,tbl[1,]>4],decreasing=T))
    sele.moa.1 <- sele.moa
    sele.moa.1[!sele.moa %in% used.moa] <- "-"
    sele.moa.tab <- table(sele.moa.1,mem)[c(used.moa,"-"),xx$rowInd]
    barplot(sele.moa.tab,beside=F,las=2,col=c(brewer.pal(9,"Paired"),"grey80"),
                horiz=T,cex.names=.3,border=NA,main="")
    dev.off()

    ##
    pdf("3_all_barplot.pdf",width=2,height=5)
    x <- table(mc$classification)[xx$rowInd]
    barplot(x,beside=F,las=2,col="grey80",horiz=T,cex.names=.3,border=NA,main="")
    dev.off()
}

moa.mem.pvals <- sapply(used.moa,function(x){
    sapply(as.character(seq(53)),function(y){
        ct <- table(sele.moa==x,mem==y)[c("TRUE","FALSE"),c("TRUE","FALSE")]
        fisher.test(ct,alternative="greater")$p.value
    })
})[rev(as.character(xx$rowInd)),]

if(0){
    ## Fig. S4C
    pdf("fdr_heatmap.pdf",width=4,height=5)
    moa.fdr <- qvalue(moa.mem.pvals)$qvalues
    cols <- rev(colorRampPalette(c("black","yellow"))(5))[c(1:5,rep(5,5))]
    image(1:10,1:53,t(moa.fdr)[ocol,rev(rownames(moa.fdr))],col=cols,
          asp=1,xlab="",ylab="",axes=F)
    par(tcl=0)
    axis(3,at=1:10,labels=colnames(moa.fdr)[ocol],las=2,cex.axis=.3)
    axis(2,at=1:53,labels=rev(rownames(moa.fdr)),cex.axis=.3,las=1)
    image(1,1:10,t(10:1),col=cols)
    dev.off()
}
