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
([Machne, Murray & Stadler 2017](http://www.nature.com/articles/s41598-017-12401-8)).


## Installation

```R
library(devtools)
install_github("raim/segmenTools")
```

... or conventionally via the source files, cloned from github.

## Capabilities

The package provides various tools for cross-correlating genomic
segmentations and annotations. 

`segmenTools' were developed specifically for analysis of genome-wide
time-series data, more specifically time series with periodic properties
such as circadian data sets. But many functionalities are broadly applicable.

### Gene-Based

#### Time-Series Analysis for Features

Fourier-based clustering of periodic time-series, after [Machne &
Murray 2012](https://doi.org/10.1371/journal.pone.0037906) and as used
in [Machne, Murray & Stadler
2017](http://www.nature.com/articles/s41598-017-12401-8).

#### Categorical Analysis

Comparing different gene categories (clusters) by cumulative
hypergeometric distribution tests, and plotting overlap enrichments.

### Coordinate-Based

#### Segment Overlap Analysis

Jaccard index statistics and relative positioning of distinct genome
segmentations (interval definitions and annotations); used in [Machne,
Murray & Stadler
2017](http://www.nature.com/articles/s41598-017-12401-8) for analysis
of segmentations by [Karl](https://github.com/raim/segmenTier).

#### Genomic Coordinate Indexing

Accessing genomic coordinates efficiently by indexing, used by
[vsRGB genome browser](https://gitlab.com/raim/genomeBrowser) and
[Karl](https://github.com/raim/segmenTier).

#### Sequence Spectra 

Analyzing periodic enrichment of oligomers and DNA structural
parameters, after [Lehmann, Machne & Herzel
2014](https://doi.org/10.1093/nar/gku641).