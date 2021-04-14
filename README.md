# BAAE_51350
Code to reproduce analysis for Callaghan et al. 2021. Urbanization negatively impacts frog diversity at continental, regional, and local scales. Basic and Applied Ecology. https://doi.org/10.1016/j.baae.2021.04.003

**NOTE** This repository will not fully reproduce the figures found in Callaghan et al. 2021. This is because all FrogID data cannot be made Open Access due to data sensitivity/privacy of the underlying recordings and localities of threatened or otherwise sensitive species (see here for details: https://zookeys.pensoft.net/article/38253/). But in order to make our workflow reproducible, we provide data and code to reproduce our analyses using the publicly available data that does not have data generalizations (i.e., removing those species that are sensitive or threatened). These data can be downloaded from here: https://www.frogid.net.au/explore. Further data can also be requested from the Australian Museum (see here for details: https://www.frogid.net.au/science).

The Australian frog list was up to date as of April 2020.

The following R information was used to run the code in this repository:

R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows Server x64 (build 17763)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252 
[2] LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] broom_0.7.3       lmerTest_3.1-3    lme4_1.1-26       Matrix_1.2-18    
 [5] purrr_0.3.4       lubridate_1.7.9.2 readr_1.4.0       patchwork_1.1.1  
 [9] phylobase_0.8.10  picante_1.8.2     nlme_3.1-149      ape_5.4-1        
[13] scales_1.1.1      tibble_3.0.5      vegan_2.5-7       lattice_0.20-41  
[17] permute_0.9-5     ggplot2_3.3.3     sf_0.9-7          tidyr_1.1.2      
[21] dplyr_1.0.3      

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6          prettyunits_1.1.1   class_7.3-17       
 [4] assertthat_0.2.1    R6_2.5.0            plyr_1.8.6         
 [7] backports_1.2.1     e1071_1.7-4         httr_1.4.2         
[10] pillar_1.4.7        rlang_0.4.10        progress_1.2.2     
[13] lazyeval_0.2.2      uuid_0.1-4          rstudioapi_0.13    
[16] minqa_1.2.4         nloptr_1.2.2.2      RNeXML_2.4.5       
[19] splines_4.0.3       statmod_1.4.35      stringr_1.4.0      
[22] munsell_0.5.0       tinytex_0.29        compiler_4.0.3     
[25] numDeriv_2016.8-1.1 xfun_0.20           pkgconfig_2.0.3    
[28] mgcv_1.8-33         tidyselect_1.1.0    XML_3.99-0.5       
[31] crayon_1.3.4        withr_2.4.1         MASS_7.3-53        
[34] grid_4.0.3          gtable_0.3.0        lifecycle_0.2.0    
[37] DBI_1.1.1           magrittr_2.0.1      units_0.6-7        
[40] KernSmooth_2.23-17  stringi_1.5.3       reshape2_1.4.4     
[43] xml2_1.3.2          ellipsis_0.3.1      generics_0.1.0     
[46] vctrs_0.3.6         boot_1.3-25         tools_4.0.3        
[49] ade4_1.7-16         rncl_0.8.4          glue_1.4.2         
[52] hms_1.0.0           parallel_4.0.3      colorspace_2.0-0   
[55] cluster_2.1.0       classInt_0.4-3     
