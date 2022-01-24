## for all classes
grps <- colnames(corrs)
system.time(es.grps <- sapply(grps,function(grp){
    pos <- names(sort(corrs[,grp],decreasing=T))
    ES.grp <- sapply(GOs.ge10,function(go){
        set.seed(12345)
        GSEA.EnrichmentScore3(pos,go)
    })
})) ## 28 sec

## GSEA using pvals.ES.GSEA(), nperm=1000 (up to p=1e-3)
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
used.db <- GOs.ge10

jns1 <- sapply(1:18,function(i){
	sorted.e <- names(sort(corrs[used.genes.ge10,i],decreasing=T))
	##
	jn <- TEA.hpc(sorted.e,used.db,nperm=10^4,wt="24:00:00",mem="2gb",
	              n.split=256,jobname=paste("clust18-",i,sep=""),verbose=T)
	return(jn)
})

# save(jns1,file="rda/jns_020115_5.rda")#10^5
    
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

# save(jns5,file="rda/jns_020115_10.5.rda")#10^5

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

# save(jns6,file="rda/jns_020115_10.6.rda")#10^6 
    
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

save(jns7,file="rda/jns_020115_10.7.rda")#10^6 ## waiting here (2/1/15)

##
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
# save(slp18.3,file="rda/slp18_3.rda")


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

save(jns8,file="rda/jns_020315_10.7.rda")#10^6 ## waiting here (2/1/15)
  
slp18.4 <- sapply(1:18,function(i){
	cat("*")
	jn <- jns8[i]
	outs <- call.TEA.out(jn)
	slps <- -log10(outs$p)* sign(outs$ES)
	this <- slp18.3[,i]
	this[totest8[[i]]] <- slps[names(GOs.ge10[totest8[[i]]])]
	return(this)
})


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
    });
    maxi <- which(m1)

    m2 <- apply(slp18.3,1,function(x){
        i1 <- which.min(x);
        if(i1==g){
          tf <- x[i1]< -own && max(x[-i1]) > -others
          return(tf)
        }else{
          return(FALSE)
        }
    });
    mini <- which(m2)
    return(list(pos=ns[m1],neg=ns[m2]))
})

names(sigs) <- colnames(gi.clust18)

sig.names <- sub("^GO_","",unique(unlist(sigs)));names(sig.names) <- c()
save(sig.names,file="rda/go.names.rda")

m <- sapply(sigs,function(x)length(unlist(x)))

##
sig.gos <- c(pos.corr.gos,neg.corr.gos)
names(sig.gos) <- c("GO:0015674","GO:0006732","GO:0006807","GO:0030001","GO:0009308","GO:0006812",
                    "GO:0015672","GO:0016616","GO:0009055","GO:0016491","GO:0016651","GO:0005216",
                    "GO:0016614","GO:0007093","GO:0002376","GO:0019221","GO:0003729")
names(pos.corr.gos) <- names(sig.gos[1:13])
names(neg.corr.gos) <- names(sig.gos[14:17])

