# prime-seq 2

This repository contains scripts used for the analysis performed in our manuscript

**Improving RNA-seq protocols**
*Felix Pförtner, Eva Briem, Wolfgang Enard, Daniel Richter*

For the full prime-seq protocol please visit [protocols.io](https://www.protocols.io/view/prime-seq-2-dsyx6fxn).

---

We introduce an approach to optimize bulk RNA-seq protocols by systematically minimizing read loss across processing steps while keeping sensitivity and complexity uncompromised. We applied this **"funnel strategy"** to prime-seq, addressing critical stages of the protocol—including DNA digestion, reverse transcription, adapter ligation, and amplification—to achieve notable efficiency enhancements. 

The optimized protocol increases final unique molecular identifiers (UMIs) by 60~\% at equal sequencing costs or, conversely, reduces sequencing costs by 38~\% at equal counts. This further improves one of the most cost-efficient bulk RNA-seq protocols available, creating **prime-seq 2**, and outlines strategies and mechanisms of potential relevance for other protocols amplifying complex nucleic acid samples for sequencing.

![Figure A: Overview of the funnel strategy applied to prime-seq protocol optimization.](figures/combined_figures/FigA.png)

---
The main code creating the figures: [Paper_Data_&_Plotting.Rmd](scripts/Paper_Data_&_Plotting.Rmd)

---
## Session Info

```r
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Devuan GNU/Linux 4 (chimaera)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.13.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Berlin
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        purrr_1.0.4        car_3.1-2          generics_0.1.4     ggpubr_0.6.0       glue_1.8.0         backports_1.5.0   
[10] rprojroot_2.0.4    scales_1.4.0       grid_4.4.1         abind_1.4-5        carData_3.0-5      tibble_3.2.1       rstatix_0.7.2      lifecycle_1.0.4    ggsignif_0.6.4    
[19] compiler_4.4.1     dplyr_1.1.4        RColorBrewer_1.1-3 pkgconfig_2.0.3    tidyr_1.3.1        here_1.0.1         rstudioapi_0.16.0  farver_2.1.2       R6_2.6.1          
[28] dichromat_2.0-0.1  tidyselect_1.2.1   pillar_1.10.2      magrittr_2.0.3     tools_4.4.1        gtable_0.3.6       broom_1.0.6        ggplot2_3.5.2   
```
