rl <- readLines("GSE5846_series_matrix.txt")
ni <- sapply(rl,function(x)length(strsplit(x,"\\t")[[1]]))
##22317 genes (not the same verison with 17,186 genes)
