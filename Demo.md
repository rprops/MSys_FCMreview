---
title: <font size = "6">Demo of microbial fcm analysis pipeline</font>
output:
  html_document:
    code_folding: show
    highlight: haddock
    keep_md: yes
    theme: flatly
    toc: yes
    number_sections: true
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 4
editor_options:
  chunk_output_type: console
---

<style type="text/css">
h1.title {
  font-size: 38px;
  color: black;
  text-align: center;
}

</style>



# Libraries


```r
library("Phenoflow")
library("plyr")
library("dplyr")
library("ggplot2")
library("flowAI")
library("scales")
library("plyr")
library("cowplot")
library("ggcyto")
library("RColorBrewer")
library("grid")
library("tidyr")
# library("Rphenograph")
library("FlowSOM")

# Set seed for reproducible analysis
set.seed(777)

# For avoiding verbose console output
run_quiet <- function (x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
```

# Preprocessing

## Import data 


```r
# Import data
flowData <- read.flowSet(path = "data/zishu_2018/")
param <- c("FS Log", "SS Log","FL 4 Log")
sampleNames(flowData) <- gsub(".5", "_5", sampleNames(flowData), fixed = TRUE)
```

## Transform data 


```r
# Transform data
flowData_transformed_log <- transform(flowData,`FS Log`=log10(`FS Log`),
                                   `SS Log`=log10(`SS Log`),
                                   `FL 4 Log`=log10(`FL 4 Log`))
flowData_transformed_asinh <- transform(flowData,`FS Log`=asinh(`FS Log`),
                                   `SS Log`=asinh(`SS Log`),
                                   `FL 4 Log`=asinh(`FL 4 Log`))
```

## Denoise data - step 1 {.tabset}

### Untransformed data 


```r
r_sam <- sample(1:length(flowData), 6)

#  scatter plot
p_scatter1 <- ggcyto::ggcyto(flowData[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter1)
```

<img src="./Figures/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

### Log-transformed data 


```r
#  scatter plot
p_scatter2 <- ggcyto::ggcyto(flowData_transformed_log[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter2)
```

<img src="./Figures/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

### asinh-transformed data 


```r
#  scatter plot
p_scatter3 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter3)
```

<img src="./Figures/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

## Denoise data - step 2 {.tabset}

### Without gate


```r
#  scatter plot
p_scatter4 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))

print(p_scatter4)
```

<img src="./Figures/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

### With gate


```r
# Create a PolygonGate for denoising the dataset
polycut <- matrix(c(1.5, 1.5, 4, 4, 9, 9, 6.75, 6.75, 6.25, 6.25,
  2.25, 6, 6, 8.75, 8.75, 2.25, 2.25, 2.75, 2.75, 2.25),
     ncol = 2,
     nrow = 10)
colnames(polycut) <- c("FS Log", "FL 4 Log")
polyGate <- polygonGate(.gate = polycut, filterId = "Cell population")

#  scatter plot
p_scatter5 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `FS Log`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Forward scatter (a.u.)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))+
  # coord_cartesian(xlim = c(0,6), ylim = c(0,6))+
  geom_gate(polyGate, col = "#CB0001", fill = "#ffa401", alpha = 0.8, size = 1)
  
print(p_scatter5)
```

<img src="./Figures/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

```r
# Retain only events in gate
flowData_transformed_asinh <- Subset(flowData_transformed_asinh, polyGate)
```

## Denoise data - step 3

Doublet/clump removal not possible for this dataset due to lack of combined
-A/-H parameters.

## Denoise data - step 4 {.tabset}

### Before flowAI


```r
#  scatter plot
p_scatter6 <- ggcyto::ggcyto(flowData_transformed_asinh[r_sam], aes(x = `TIME`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Time (ms)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))

print(p_scatter6)
```

<img src="./Figures/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />



### After flowAI 


```r
#  scatter plot
p_scatter7 <- ggcyto::ggcyto(flowData_transformed_asinh_dn[r_sam], aes(x = `TIME`, y = `FL 4 Log`)) + 
  geom_hex(bins = 300) +
  theme_bw()+
  labs(x = "Time (ms)", y = "DAPI (a.u.)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "white"))

print(p_scatter7)
```

<img src="./Figures/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

### QC stats


```r
QC_stats <- read.delim("./data/QC_flowAI/QCmini.txt")

p_qc1 <- QC_stats %>%
  arrange(-X..anomalies) %>% 
  mutate(Name.file = factor(Name.file, levels = unique(Name.file))) %>% 
  ggplot(., aes(y = X..anomalies, x = Name.file))+
  geom_point(shape = 21, size = 3, fill = "#ffa401", alpha = 0.5)+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
        axis.title = element_text(size = 12))+
  labs(y = "% anomalies",  x = "Ranked samples")+
  geom_hline(yintercept = 5, linetype = 2)

p_qc2 <- QC_stats %>% 
  ggplot(., aes(x = X..anomalies))+
  geom_histogram(fill = "#ffa401", alpha = 0.5, col = "black")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(x = "% anomalies",  y = "Number of samples")+
  geom_vline(xintercept = median(QC_stats$X..anomalies), linetype = 2)+
  scale_x_continuous(breaks = seq(from = 0, to = 100, by = 10), limits = c(-5,100))

cowplot::plot_grid(p_qc1, p_qc2, align = "hv", ncol = 1)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 2 rows containing missing values (geom_bar).
```

<img src="./Figures/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

# Cell density measurement

Cell concentration cannot be estimated for this dataset due to lack of 
volumetric measurement.



# Cell population identification {.tabset}

## FlowSOM - fig1

```r
fSOM <- FlowSOM(
  input = flowData_transformed_asinh[1:3],
  # Input options:
  compensate = FALSE,
  transform = FALSE,
  scale = FALSE,
  # SOM options:
  colsToUse = param,
  xdim = 7,
  ydim = 7,
  nClus = 10
  )

PlotStars(fSOM[[1]], backgroundValues = as.factor(fSOM[[2]]))
```

<img src="./Figures/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

## FlowSOM - fig2

```r
PlotMarker(fSOM[[1]], param[3])
```

<img src="./Figures/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

# Cytometric fingerprinting



# Conclusions



# Session info


```r
sessionInfo()
```

```
## R version 4.0.1 (2020-06-06)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17134)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Dutch_Belgium.1252  LC_CTYPE=Dutch_Belgium.1252   
## [3] LC_MONETARY=Dutch_Belgium.1252 LC_NUMERIC=C                  
## [5] LC_TIME=Dutch_Belgium.1252    
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] FlowSOM_1.20.0            igraph_1.2.5             
##  [3] tidyr_1.1.0               RColorBrewer_1.1-2       
##  [5] ggcyto_1.16.0             flowWorkspace_4.0.6      
##  [7] ncdfFlow_2.34.0           BH_1.72.0-3              
##  [9] RcppArmadillo_0.9.900.2.0 cowplot_1.0.0            
## [11] scales_1.1.1              ggplot2_3.3.2            
## [13] dplyr_1.0.0               plyr_1.8.6               
## [15] Phenoflow_1.1.2           foreach_1.5.0            
## [17] flowAI_1.18.5             flowFDA_0.99             
## [19] mclust_5.4.6              multcomp_1.4-13          
## [21] TH.data_1.0-10            MASS_7.3-51.6            
## [23] survival_3.2-3            mvtnorm_1.1-1            
## [25] flowFP_1.46.0             flowViz_1.52.0           
## [27] lattice_0.20-41           flowClean_1.26.0         
## [29] flowCore_2.0.1           
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.4-1            ellipsis_0.3.1             
##  [3] class_7.3-17                cytolib_2.0.3              
##  [5] XVector_0.28.0              base64enc_0.1-3            
##  [7] farver_2.0.3                hexbin_1.28.1              
##  [9] IDPmisc_1.1.20              CytoML_2.0.5               
## [11] prodlim_2019.11.13          lubridate_1.7.9            
## [13] xml2_1.3.2                  codetools_0.2-16           
## [15] splines_4.0.1               knitr_1.29                 
## [17] ade4_1.7-15                 jsonlite_1.7.0             
## [19] phyloseq_1.32.0             pROC_1.16.2                
## [21] caret_6.0-86                cluster_2.1.0              
## [23] png_0.1-7                   sfsmisc_1.1-7              
## [25] graph_1.66.0                compiler_4.0.1             
## [27] Matrix_1.2-18               htmltools_0.5.0            
## [29] tools_4.0.1                 gtable_0.3.0               
## [31] glue_1.4.1                  reshape2_1.4.4             
## [33] Rcpp_1.0.5                  Biobase_2.48.0             
## [35] vctrs_0.3.2                 Biostrings_2.56.0          
## [37] multtest_2.44.0             ape_5.4                    
## [39] nlme_3.1-148                iterators_1.0.12           
## [41] changepoint_2.2.2           timeDate_3043.102          
## [43] gower_0.2.2                 xfun_0.16                  
## [45] stringr_1.4.0               lifecycle_0.2.0            
## [47] XML_3.99-0.5                zlibbioc_1.34.0            
## [49] zoo_1.8-8                   ipred_0.9-9                
## [51] RProtoBufLib_2.0.0          RBGL_1.64.0                
## [53] parallel_4.0.1              biomformat_1.16.0          
## [55] sandwich_2.5-1              rhdf5_2.32.2               
## [57] yaml_2.2.1                  gridExtra_2.3              
## [59] rpart_4.1-15                latticeExtra_0.6-29        
## [61] stringi_1.4.6               S4Vectors_0.26.1           
## [63] permute_0.9-5               BiocGenerics_0.34.0        
## [65] boot_1.3-25                 lava_1.6.7                 
## [67] rlang_0.4.7                 pkgconfig_2.0.3            
## [69] matrixStats_0.56.0          evaluate_0.14              
## [71] purrr_0.3.4                 Rhdf5lib_1.10.1            
## [73] labeling_0.3                recipes_0.1.13             
## [75] bit_1.1-15.2                tidyselect_1.1.0           
## [77] magrittr_1.5                R6_2.4.1                   
## [79] IRanges_2.22.2              generics_0.0.2             
## [81] pillar_1.4.6                withr_2.2.0                
## [83] mgcv_1.8-31                 nnet_7.3-14                
## [85] tsne_0.1-3                  tibble_3.0.3               
## [87] crayon_1.3.4                KernSmooth_2.23-17         
## [89] rmarkdown_2.3               jpeg_0.1-8.1               
## [91] data.table_1.13.0           vegan_2.5-6                
## [93] Rgraphviz_2.32.0            ConsensusClusterPlus_1.52.0
## [95] ModelMetrics_1.2.2.2        digest_0.6.25              
## [97] RcppParallel_5.0.2          stats4_4.0.1               
## [99] munsell_0.5.0
```
