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



# Preprocessing

## Import, normalize & denoise data


```r
# Create a PolygonGate for denoising the dataset
sqrcutWater<- matrix(c(7.25, 7.25, 14, 14,
     2.5, 5.5, 13, 2.5),
     ncol = 2,
     nrow = 4)
colnames(sqrcutWater) <- c("BL1-H", "BL3-H")
polyGateWater <- polygonGate(.gate = sqrcutWater, filterId = "Total Cells")

# Import & preprocess
```

## Extract metadata


```r
# Create metadata from sample names
```

# Cell density assessment



# Cell population identification



# Cytometric fingerprinting



# Conclusions


