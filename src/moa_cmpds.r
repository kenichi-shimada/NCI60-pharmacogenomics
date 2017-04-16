##
## data imported
## mc

dir <- "~/Projects/TF_TEA/505_sd7/nci60_analysis/"

##
setwd("/Users/ks2474/Projects/CIL56/6_genomics/175_NCI60_analysis/pharmacology_data/old/original_data")
load("GI50-75446.rda")

setwd(dir)
load("nci60-analysis.rda")

## classification
cl <- mc$classification ##2565 -> 2555


## 
uniq.gi50.info <- unique(gi50.info[c(1,4)])
moa.all <- c(uniq.gi50.info$Mech,rep("FIN",10))
names(moa.all) <- c(uniq.gi50.info$NSC,tail(names(cl),10))
moa.all[moa.all %in% c("-","","Unknown")] <- "-"
fin.iqr <- rep(c(FALSE,TRUE),c(1,9));names(fin.iqr) <- tail(names(cl),10)
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

## MOA-statistics
setwd(dir);setwd("figs_moa")
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
par(mar=c(7,5,1,1))
plot(logp,axes=F,xlab="",ylab="-log(p-value)",type="n")
box()
axis(1,at=seq(logp),labels=names(logp),las=2,cex.axis=.7)
axis(2)
abline(h=seq(0,15,5),col="grey80")
abline(v=15,col=2)
points(logp,pch=20)
which(logp>2) ## "Ds","T1","STK"
d2 <- table(moa.all=="T1",has.large.iqr)[c("TRUE","FALSE"),c("TRUE","FALSE")]
d2 <- table(moa.all=="STK",has.large.iqr)[c("TRUE","FALSE"),c("TRUE","FALSE")]
d2 <- table(moa.all=="Ds",has.large.iqr)[c("TRUE","FALSE"),c("TRUE","FALSE")]
t1 <- names(which(moa.all=="T1"&has.large.iqr))
identical(names(sele.moa),names(cl))

##
setwd("~/Projects/TF_TEA/505_sd7/nci60_analysis/figs_moa")
##
##
## hist(moa.fdr)
ocol <- order(apply(moa.cl.pvals,2,which.min))
sort.moa.pval <- moa.cl.pvals[,ocol]
sort.moa <- moall.tab[used.moa,]
barplot(sort.moa,beside=F,las=2,col=c(brewer.pal(12,"Paired")[-11],"grey90"),
        horiz=T,cex.names=.3,border=NA,
        main="")
rowSums(moall.tab)
dev.off()
##
setwd("/Users/ks2474/Projects/TF_TEA/505_sd7/nci60_analysis/figs_moa")
pdf("1d_moa_barplot_1.pdf",width=2,height=5)
for(moa in used.moa){
    i <- match(moa,used.moa)
    if(all(sele.moa!=moa))next;
    moa.table <- table(sele.moa==moa,cl)[c("TRUE","FALSE"),]
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
sele.moa.tab <- table(sele.moa.1,cl)[c(used.moa,"-"),xx$rowInd]
barplot(sele.moa.tab,beside=F,las=2,col=c(brewer.pal(9,"Paired"),"grey80"),
            horiz=T,cex.names=.3,border=NA,main="")
dev.off()
##
pdf("3_all_barplot.pdf",width=2,height=5)
x <- table(mc$classification)[xx$rowInd]
barplot(x,beside=F,las=2,col="grey80",horiz=T,cex.names=.3,border=NA,main="")
dev.off()

##
hist(nlp,breaks=40,ylim=c(0,10))
##
mat <- moall.tab.ge3
nlp <- t(sapply(seq(nrow(mat)),function(x){
    sapply(seq(ncol(mat)),function(y){
        a <- mat[x,y]
        b <- sum(mat[-x,y])
        c <- sum(mat[x,-y])
        d <- sum(mat[-x,-y])
        ctg <- matrix(c(a,b,c,d),nrow=2)
        -log10(fisher.test(ctg,alternative="greater")$p.value)
    })
}))
rownames(nlp) <- rownames(mat)
colnames(nlp) <- colnames(mat)
cols <- colorRampPalette(c("black","yellow"))(10)
nlp[nlp>10] <- 10;nlp <- nlp[-12,]
pdf("1f_moa_matrix.pdf",width=4,height=5)
image(1:11,1:53,nlp[,as.character(xx$rowInd)],col=cols,asp=1,xlab="",ylab="",axes=F)
par(tcl=0)
axis(3,at=1:11,labels=rownames(nlp),las=2,cex.axis=.3)
axis(2,at=1:53,labels=xx$rowInd,cex.axis=.3,las=1)
dev.off()
##

moa.cl.pvals <- sapply(used.moa,function(x){
    sapply(as.character(seq(53)),function(y){
        ct <- table(sele.moa==x,cl==y)[c("TRUE","FALSE"),c("TRUE","FALSE")]
        fisher.test(ct,alternative="greater")$p.value
    })
})[rev(as.character(xx$rowInd)),]
library(qvalue)
pdf("fdr_heatmap.pdf",width=4,height=5)
moa.fdr <- qvalue(moa.cl.pvals)$qvalues
cols <- rev(colorRampPalette(c("black","yellow"))(5))[c(1:5,rep(5,5))]
image(1:10,1:53,t(moa.fdr)[ocol,rev(rownames(moa.fdr))],col=cols,
      asp=1,xlab="",ylab="",axes=F)
par(tcl=0)
axis(3,at=1:10,labels=colnames(moa.fdr)[ocol],las=2,cex.axis=.3)
axis(2,at=1:53,labels=rev(rownames(moa.fdr)),cex.axis=.3,las=1)
image(1,1:10,t(10:1),col=cols)
dev.off()


##
t1.clust <- c("37","40","42")
names(sig.gos.all) <- seq(53)
t1.pos <- unlist(sapply(sig.gos.all[t1.clust],function(x)x$pos))
t1.neg <- unlist(sapply(sig.gos.all[t1.clust],function(x)x$neg))
##
es.t1 <- t(es.grps[c(t1.pos,t1.neg),])
cols <- rep(1,53);cols[as.numeric(t1.clust)] <- 2
plot(c(1,18),range(es.t1),type="n")
abline(h=0)
for(i in 1:18){
    points(jitter(rep(i,50),amount=.2),es.t1[cols==1,i],pch=20,col=1,cex=.7)
    points(jitter(rep(i,3),amount=.2),es.t1[cols==2,i],pch=20,col=2,cex=.7)
}
##
head(gi50.info)
##
norm.t.stat <- apply(t.stat,2,function(x)x/sum(x))
barplot(t.stat,col=rainbow(53),las=2,ylim=c(0,30))
barplot(norm.t.stat,col=rainbow(53),las=2)

##
dtp2sid <- function(dtp){
    q <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/DTP.NCI/",
               dtp,"/sids/txt",sep="")
    z <- getURL(q)
    if(length(grep("^Status",z))){
        return(NA)
    }else{
        return(sub("\n$","",z))
    }
}
sele.cmpds <- names(which(has.large.iqr))
sids <- rep(NA,length(sele.cmpds))
names(sids) <- sele.cmpds
for(i in seq(sele.cmpds)){
    id <- sele.cmpds[i]
    sids[id] <- dtp2sid(id)
    if(!i %% 50){
        cat(paste("* ",i,"\n",sep=""))
    }else if(!i %% 10){
        cat("* ")
    }else{
        cat("*")
    }
}
##
sids[grep("^Status",sids)] <- NA
##
sid2cid <- function(sid){
    q <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/",sid,"/cids/txt",sep="")
    z <- getURL(q)
    if(length(grep("^Status",z))){
        return(NA)
    }else{
        return(sub("\n$","",z))
    }
}
cid2smiles <- function(cid){
    q <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
               cid,"/property/CanonicalSMILES/txt",sep="")
    z <- getURL(q)
    if(length(grep("^Status",z))){
        return(NA)
    }else{
        return(sub("\n$","",z))
    }
}
cid2tgts <- function(cid){
    require(XML)
    q <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
               cid,"/xrefs/RegistryID/XML",sep="")
    z <- getURL(q)
    ##if(length(grep("^Status",z))){return(NA)}else{return(sub("\n$","",z))}
    regid <- sub(" +<RegistryID>(.+)</Reg.+","\\1",
                 grep("RegistryID",strsplit(z,"\\n")[[1]],value=T))
    dbs <- grep("^DB",regid,value=T)
    if(!length(dbs)){
        return(NA)
    }else{
        tgts <- lapply(dbs,function(db){
            url <- paste("http://www.drugbank.ca/drugs/",db,sep="")
            html <- htmlTreeParse(getURL(url),useInternalNodes=T)
            tgts <- xpathSApply(html,"//div[@class='target well well-sm bonds']/div//tbody/tr/td",
                                xmlValue)
            tgts <- sub(" +$","",tgts)
            return(matrix(tgts,byrow=T,ncol=3)[,1:2])
        })
        if(length(dbs)>0)tgts <- unique(do.call(rbind,tgts))
        if(nrow(tgts)==0){
            return(NA)
        }else{
            return(tgts)
        }
    }
}

##
cids <- sids
smiles <- cids;smiles[] <- NA
smiles[] <- cids[] <- NA
for(i in seq(sids)){
    sid <- sids[i]
    if(is.na(sid))next;
    cids[i] <- sid2cid(sid)
    if(!i %% 50){
        cat(paste("* ",i,"\n",sep=""))
    }else if(!i %% 10){
        cat("* ")
    }else{
        cat("*")
    }
}
##
cids[grep("^Status",cids)] <- NA
for(i in seq(cids)){
    if(!i %% 50){
        cat(paste("* ",i,"\n",sep=""))
    }else if(!i %% 10){
        cat("* ")
    }else{
        cat("*")
    }
    cid <- cids[i]
    if(is.na(cid))next;
    smiles[i] <- cid2smiles(cid)
}

pchem <- cbind(sele.cmpds,cids,sids,cl)
save(pchem,file="pubchem_data.rda")
write.csv(pchem,"pubchem_data.csv")
##
tgts <- smiles;tgts[] <- NA; tgts <- as.list(tgts)
for(i in seq(tgts)){
    if(!i %% 50){
        cat(paste("* ",i,"\n",sep=""))
    }else if(!i %% 10){
        cat("* ")
    }else{
        cat("*")
    }
    cid <- cids[i]
    tgts[[i]] <- try(cid2tgts(cid),TRUE)
}
nr <- sapply(tgts,nrow)
hist(unlist(nr[!sapply(nr,is.null)]))
##
unique(unlist(sapply(tgts,function(x)if(is.na(x)){NA}else{x[,1]})))
