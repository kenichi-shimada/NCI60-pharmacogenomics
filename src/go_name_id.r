setwd(dir <- "~/Projects/TF_TEA/505_sd7/nci60_analysis/")



load("go.names.rda") ## where is it???






##
library(RCurl)
gos <- rep(NA,length(sig.names));names(gos) <- sig.names
for(sn in sig.names){
    i <- match(sn,sig.names)
    if(!i %% 50){
        cat("*\n")
    }else if(!i %% 10){
        cat("* ")
    }else{
        cat("*")
    }
    url <- paste("www.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=",sn,sep="")
    html <- htmlTreeParse(getURL(url),useInternalNodes=T)
    tgts <- xpathSApply(html,"//table[@class='lists4']",xmlValue)
    gos[sn] <- sub(".+(GO:[0-9]+)\\..+","\\1",tgts)
}
names(gos) <- sig.names

## all GOID
gos.ge10 <- sub("GO_","",names(GOs.ge10))
gos.all <- rep(NA,length(gos.ge10));names(gos.all) <- gos.ge10
for(sn in gos.ge10){
    i <- match(sn,gos.ge10)
    if(!i %% 50){
        cat("*\n")
    }else if(!i %% 10){
        cat("* ")
    }else{
        cat("*")
    }
    url <- paste("www.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=",sn,sep="")
    html <- htmlTreeParse(getURL(url),useInternalNodes=T)
    tgts <- xpathSApply(html,"//table[@class='lists4']",xmlValue)
    gos.all[sn] <- sub(".+(GO:[0-9]+)\\..+","\\1",tgts)
}
names(gos.all) <- c()

##
library(GO.db)

goterms <- mget(gos,GOTERM,ifnotfound=NA)
goterms <- goterms[!sapply(goterms,is.na)]
mget(names(goterms[sapply(goterms,is.na)]),GOOBSOLETE,ifnotfound=NA)
notfound <- sig.names[which(sapply(goterms,is.na))]
apply(slp18.3[match(notfound,sub("^GO_","",names(GOs.ge10))),],1,range)
head(gt)

##
onts <- sapply(goterms,Ontology)


##
sig.got <- lapply(sigs,function(x){
    go.pos <- gos[sub("^GO_","",x$pos)]
    go.neg <- gos[sub("^GO_","",x$neg)]
    return(list(pos=go.pos,neg=go.neg))
})

sig.cats <- lapply(sigs,function(x){
    list(pos=mget(gos[sub("^GO_","",x$pos)],
             neg=gos[sub("^GO_","",x$neg)]))
})

## BP
bp.pos <- lapply(sig.got,function(x){
    y <- x$pos
    ont <- Ontology(y)
    y[ont=="BP" & !is.na(ont)]
})
bp.neg <- lapply(sig.got,function(x){
    y <- x$neg
    ont <- Ontology(y)
    y[ont=="BP" & !is.na(ont)]
})

## MF
mf.pos <- lapply(sig.got,function(x){
    y <- x$pos
    ont <- Ontology(y)
    y[ont=="MF" & !is.na(ont)]
})
mf.neg <- lapply(sig.got,function(x){
    y <- x$neg
    ont <- Ontology(y)
    y[ont=="MF" & !is.na(ont)]
})

##CC
cc.pos <- lapply(sig.got,function(x){
    y <- x$pos
    ont <- Ontology(y)
    y[ont=="CC" & !is.na(ont)]
})
cc.neg <- lapply(sig.got,function(x){
    y <- x$neg
    ont <- Ontology(y)
    y[ont=="CC" & !is.na(ont)]
})

##
bp.pos.off <- lapply(bp.pos,function(x){
    pa <- mget(x,GOBPANCESTOR)
    pas <- lapply(pa,function(y)y[y%in% x]) ## remove the nodes
    return(x[!x %in% unlist(pas)])
})
bp.neg.off <- lapply(bp.neg,function(x){
    pa <- mget(x,GOBPANCESTOR)
    pas <- lapply(pa,function(y)y[y%in% x]) ## remove the nodes
    return(x[!x %in% unlist(pas)])
})
mf.pos.off <- lapply(mf.pos,function(x){
    pa <- mget(x,GOMFANCESTOR)
    pas <- lapply(pa,function(y)y[y%in% x]) ## remove the nodes
    return(x[!x %in% unlist(pas)])
})
mf.neg.off <- lapply(mf.neg,function(x){
    pa <- mget(x,GOMFANCESTOR)
    pas <- lapply(pa,function(y)y[y%in% x]) ## remove the nodes
    return(x[!x %in% unlist(pas)])
})
cc.pos.off <- lapply(cc.pos,function(x){
    pa <- mget(x,GOCCANCESTOR)
    pas <- lapply(pa,function(y)y[y%in% x]) ## remove the nodes
    return(x[!x %in% unlist(pas)])
})
cc.neg.off <- lapply(cc.neg,function(x){
    pa <- mget(x,GOCCANCESTOR)
    pas <- lapply(pa,function(y)y[y%in% x]) ## remove the nodes
    return(x[!x %in% unlist(pas)])
})
##
library("Rgraphviz")
rEG <- new("graphNEL", nodes=c("A", "B"), edgemode="directed")
rEG <- addEdge("A", "B", rEG, 1)
plot(rEG)

##
setwd("/Users/ks2474/Dropbox/NCI60 analysis/gos")
save.image("02-04-15-GOs.rda")

gocc <- as.list(GOCCOFFSPRING)
length(gocc)
