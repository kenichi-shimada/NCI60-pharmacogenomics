# NCI-60-pharmacogenomics
This repo contains the codes and data to generate figures for the following paper:
- Shimada *et al.*, Cell-Line Selectivity Improves the Predictive Power of Pharmacogenomic Analyses and Helps Identify NADPH as Biomarker for Ferroptosis Sensitivity. __*Cell Chemical Biology*__, 2016.

### Overview
Please see [the html](rmd/main_nci60_analysis_all_clusts_012017.html) (generated from [this rmd file](rmd/main_nci60_analysis_all_clusts_012017.rmd)) for the analysis overview. 

### Full analysis
See [main_nci60_analysis.r](main_nci60_analysis.r) for the full analysis.

### R `sessionInfo()`

Note that some packages (e.g., *mclust* package) requires to use the specific version. The latest version couldn't recapitulate some of my analysis. For switching older versions of R and using a specific library, refer to [this](https://github.com/kenichi-shimada/Simple_R_envs).

	## R version 3.3.2 (2016-10-31)
	## Platform: x86_64-apple-darwin13.4.0 (64-bit)
	## Running under: OS X El Capitan 10.11.6
	## 
	## locale:
	## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
	## 
	## attached base packages:
	## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
	## [8] methods   base     
	## 
	## other attached packages:
	##  [1] org.Hs.eg.db_3.4.0   AnnotationDbi_1.36.1 IRanges_2.8.1       
	##  [4] S4Vectors_0.12.1     Biobase_2.34.0       BiocGenerics_0.20.0 
	##  [7] robustbase_0.92-7    mclust_5.2           RColorBrewer_1.1-2  
	## [10] gplots_3.0.1         rmarkdown_1.3        knitr_1.15.1        
	## 
	## loaded via a namespace (and not attached):
	##  [1] Rcpp_0.12.8        magrittr_1.5       stringr_1.1.0     
	##  [4] caTools_1.17.1     tools_3.3.2        KernSmooth_2.23-15
	##  [7] DBI_0.5-1          htmltools_0.3.5    gtools_3.5.0      
	## [10] yaml_2.1.14        rprojroot_1.1      digest_0.6.10     
	## [13] bitops_1.0-6       memoise_1.0.0      RSQLite_1.1-2     
	## [16] evaluate_0.10      gdata_2.17.0       stringi_1.1.2     
	## [19] DEoptimR_1.0-8     backports_1.0.4