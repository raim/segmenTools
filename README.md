# Coordinate Indexing and Overlap Analysis of Genomic Features 

This R library includes various utilities to handle genomic feature
coordinates, analyze overlaps between genomic segments, or overlaps
between categorical annotations (GO etc) of genomic features.

It was used for analysis of results from
[Karl, the segmenTier](https://github.com/raim/segmenTier), a
(genomic) segmentation algorithm working with abstract similarities,
e.g., derived from RNA-seq time series
([Machne, Murray & Stadler 2017](www.nature.com/articless41598-017-12401-8)).

And it provides coordinate indexing and sequencing/feature handling
utilities for the
[vsRGB genome browser](https://gitlab.com/raim/genomeBrowser).

## Quick Guide

```R
library(devtools)
install_github("raim/segmenTools")
```
