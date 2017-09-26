# Genomic Feature & Coordinate Utilities 

This R library includes various utilities to handle genomic feature
coordinates, analyze overlaps between genomic segments, or overlaps
between categorical annotations (GO etc) of genomic features.

The library provides coordinate indexing and sequencing/feature
handling utilities for the
[vsRGB genome browser](https://gitlab.com/raim/genomeBrowser).

The repository currently also holds the command-line scripts
(directory `scripts`) that were used for running and analyses of results
from [Karl, the segmenTier](https://github.com/raim/segmenTier), a
(genomic) segmentation algorithm working with abstract similarities,
e.g., derived from RNA-seq time series
([Machne, Murray & Stadler 2017](www.nature.com/articless41598-017-12401-8)).


## Quick Guide

```R
library(devtools)
install_github("raim/segmenTools")
```
