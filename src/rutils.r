## This file is called from 'main_nci60_analysis.r'. 

## hierarchical clustering and distance computing functions
hclustfun <- function(x)hclust(x,method="average")
distfun <- function(x)as.dist(1-cor(t(x),method="pearson",use="everything"))

## a wrapper of heatmap.2 function in 'gplots' package
myheatmap <- function(mat,...){
  require(gplots)
  rn <- range(mat,na.rm=T)
  rng <- round(rn/max(abs(rn))*256)
  cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(513)##[match(seq(rng[1],rng[2]),-256:256)]
  heatmap.2(t(mat),trace="none",col=cols,na.color="grey20",
            hclustfun=hclustfun,distfun=distfun,
            margins=c(5,7),density.info="none",...) 
}

## The following functions (GSEA.EnrichmentScore and GSEA.EnrichmentScore2) were cloned
## from "GSEA-P-R" package at the Broad Institute's website.
GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  

   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 
   if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
   }
   alpha <- weighted.score.type
   correl.vector <- abs(correl.vector**alpha)
   sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
   norm.tag    <- 1.0/sum.correl.tag
   norm.no.tag <- 1.0/Nm
   RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
   max.ES <- max(RES)
   min.ES <- min(RES)
   if (max.ES > - min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
   } else {
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
   }
   return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

GSEA.EnrichmentScore2 <- function (gene.list, gene.set, weighted.score.type = 0, correl.vector = NULL)
{
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  loc.vector <- vector(length = N, mode = "numeric")
  peak.res.vector <- vector(length = Nh, mode = "numeric")
  valley.res.vector <- vector(length = Nh, mode = "numeric")
  tag.correl.vector <- vector(length = Nh, mode = "numeric")
  tag.diff.vector <- vector(length = Nh, mode = "numeric")
  tag.loc.vector <- vector(length = Nh, mode = "numeric")
  loc.vector[gene.list] <- seq(1, N) 
  tag.loc.vector <- loc.vector[gene.set] 
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector] * 
      correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]^weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  norm.tag <- 1/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
  return(list(ES = ES))
}

## The following function is a modification of GSEA.EnrichmentScore2 for faster computation
GSEA.EnrichmentScore3 <- function (gene.list, gene.set)  ##weighted.score.type=0
{
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    tag.loc.vector <- sort(match(gene.set,gene.list),decreasing=F)
    tag.correl.vector <- rep(1/Nh,Nh)
    tag.diff.vector <- (diff(c(0,tag.loc.vector))-1)/Nm ## didn't count the genes in gene.set themselves
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector ## count one value before peaks
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
    return(ES)
}

##
pvals.ES.GSEA <- function(gene.list,gene.set,nperm=1000){
  ## set seed
  original.ES <- GSEA.EnrichmentScore3(gene.list,gene.set)
  set.seed(12345)
  perm.ES <- unlist(lapply(seq(nperm),function(i)
                           GSEA.EnrichmentScore3(sample(gene.list),gene.set)))
  p.value <- min(sum(perm.ES > original.ES),sum(perm.ES < original.ES))/nperm
  return(list(ES=original.ES,p.value=p.value))
}

##
TEA.hpc2 <- function (sorted.genes, target.sets, nperm = 10^6, all = F,
                      weighted.score.type=0, ## 0:NULL or 1:actvity.scores to use as correl.vector
                      ## parameters related to TEA up to here
                      n.split=1, ## n.split is explicitly specified in TEA
                      src.dir="/Users/kenichi/Desktop/pharmacogenomics/nci60_analysis/hpc_temp/hpc_TEA", 
                      r.file="temp.TEA2.R",
                      sh.file="temp.sh1.txt",
                      jobname=NULL,nodes=1,ppn=1,wt="24:00:00",mem="2gb",
                      username="ks2474",
                      server="hpcsubmit.cc.columbia.edu",
                      remote.dir="~/hpc_wd/TEA/",
                      verbose=FALSE,
                      scp=TRUE,...) {
    ## store working directory, wd
    wd <- getwd()
    
    ## create 'jobname' directory as root for the project
    if(verbose)cat("Creating directories...")
    tstamp <- sub(" ([0-9]{6})","-\\1",gsub(":","",Sys.time()))
    if(is.null(jobname)){
        jobname <- tstamp
    }else{
        jobname <- paste(tstamp,"-",jobname,sep="")
    }

    ## 
    setwd(tempdir())
    if(!jobname %in% dir())dir.create(jobname)
    setwd(tmpdir <- paste(tempdir(),"/",jobname,sep=""))
    dir.create("log")
    dir.create("src")
    dir.create("input_data")
    dir.create("output_data")
    if(verbose)cat("Done.\n")
    
    ## read script templates
    if(verbose)cat("Loading tempalte files...")
    setwd(src.dir)
    r.src <- readLines(r.file)
    sh <- readLines(sh.file)
    if(verbose)cat("Done.\n")
    
    ## R script - no modification (level 2)
    if(verbose)cat("Saving src files...")
    setwd(tmpdir);setwd("src")
    r.src <- sub("JOBNAME",jobname,r.src)
    r.src <- sub("NPERM",nperm,r.src)
    writeLines(r.src,r.file) ## done
    
    ## Set parameters for Torque
    n.proc <- nodes * ppn
    ##mc <- n.proc - 1
    vars <- c(jobname,nodes,ppn,wt,mem,n.split)
    names(vars) <- c("JOBNAME","NODES","PPN","WT","MEM","NSPLIT")
    for(i in seq(vars)){
        a <- names(vars)[i]
        b <- vars[i]
        sh <- gsub(a,b,sh)
    }
    
    ## Writing master shell script (level 1)
    setwd(tmpdir) ;setwd("src")
    writeLines(sh,"master-TEA.sh") ## save shell script w mail
    if(verbose)cat("Done.\n")
    
    ## split list of target.sets into n.split vars
    ## save rda files
    if(verbose)cat("Saving \"input.rda\"...")
    setwd(tmpdir);setwd("input_data")
    save(sorted.genes,file="input.rda")
    if(verbose)cat("Done.\n")
    
    if(verbose)cat("Splitting and Saving \"target sets\"...")
    tgt.sets <- target.sets
    save(tgt.sets,file="tgts.rda")
    
    if(verbose)cat("Done.\n")
    
    ## Configuration log
    if(verbose)cat("Writing \"init.conf\" under \"log\" dir...")
    setwd(tmpdir);setwd("log")
    writeLines(paste("n.split=",n.split,
                     "\nweighted.score.type=",weighted.score.type,"\n",sep=""),"init.conf")
    if(verbose)cat("Done.\n")
    
    ## compress jobname directory and send it to the server (if scp is TRUE)
    setwd(tempdir())
    if(scp){
        if(verbose)cat("SCPing to the server...\n")
        remote <- paste(username,"@",server,":",remote.dir,sep="")
        system(paste("scp -rC ",jobname," ",remote,sep=""),ignore.stdout=(!verbose))
        
        ## Running master.sh
        if(verbose)cat("Start \"master-TEA.sh\" execution...\n")
        system(paste("ssh ",username,"@",server," \"qsub ",remote.dir,"/",jobname,"/src/master-TEA.sh\"",sep=""),
               ignore.stdout=(!verbose))
        ##system(paste("ssh ",username,"@",server," id",sep=""));return();
        if(verbose)cat("Done.\n")
    }
    
    ## Return absolute path to output directory
    setwd(wd)
    if(verbose)cat("**********************************************\n")
    if(verbose)cat("\"qsub master-TEA.sh\" executed.\n")
    if(verbose)cat("**********************************************\n")
    return(jobname)
}
