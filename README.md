# Genomic Feature & Coordinate Utilities 

The package provides various tools for cross-correlating genome
segmentations and annotations. 


`segmenTools' were developed specifically for analysis of genome-wide
time-series data, more specifically time series with periodic properties
such as circadian data sets. But many functionalities are broadly applicable.

Its coordinate indexing and feature annotation
utilities used in various publications (see *Capabilities*), and by
[u'r gene bro](https://gitlab.com/raim/genomeBrowser).

The git repository also holds the command-line scripts (directory
`scripts`) that were used for running and analyses of results from
[Karl, the segmenTier](https://github.com/raim/segmenTier), a
(genomic) segmentation algorithm working with abstract similarities,
e.g., derived from RNA-seq time series ([Machne, Murray & Stadler
2017](http://www.nature.com/articles/s41598-017-12401-8)).

The drawing is the most unconstrained method of modeling in biology,
therefore many functionalities in `segmenTools' provide exploratory as well
as publication-quality plotting utilities.

## Installation

```R
library(devtools)
install_github("raim/segmenTools")
```

... or conventionally via the source files, cloned from github.

## Capabilities

### Gene-Based

#### Time-Series Analysis 

Via [Karl](https://github.com/raim/segmenTier): Fourier-based
clustering of periodic time-series, after [Machne & Murray
2012](https://doi.org/10.1371/journal.pone.0037906) and as extended in
[Machne, Murray & Stadler
2017](http://www.nature.com/articles/s41598-017-12401-8) for
similarity-based segmentation of coordinate-based time-series
(RNA-seq).

TODO: Cluster-wise oscillation parameters

```R
library(segmenTier) # for clustering 
library(segmenTools) # for plots

## download & parse data
rawdata.url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5612/matrix/GSE5612_series_matrix.txt.gz"
rawdata <- gsub( ".*/","",rawdata.url)
if ( !file.exists(rawdata) )
  utils::download.file(url=rawdata.url, dest=rawdata)
dat <- read.delim(gzfile(rawdata),comment.char="!",row.names=1)

## process time-series (Discrete Fourier Transform)
tset <- processTimeseries(dat, use.fft=TRUE, dc.trafo="ash",use.snr=TRUE)
## cluster (by kmeans)
cset <- clusterTimeseries(tset,K=7) # CLUSTERING! takes a while

## and inspect clustered time-series via the versatile
## cluster time series plotter
pdf("edwards06.pdf")
plotClusters(tset, cset, norm="lg2r", each=TRUE, type="all", ylim="all")
plotClusters(tset, cset, norm="lg2r", each=TRUE, q=0.8)
## selected clusters in all-in-one plo
plotClusters(tset, cset, norm="lg2r", each=FALSE, type="rng", cls.srt=c(3,5,7))
dev.off()
```

#### Categorical Analysis 

Comparing different gene categories (clusters) by cumulative
hypergeometric distribution tests, and plotting overlap enrichments
after [Machne & Murray
2012](https://doi.org/10.1371/journal.pone.0037906).

### Coordinate-Based

#### Segment Overlap Analysis

Jaccard index statistics and relative positioning of distinct genome
segmentations (interval definitions and annotations); used in [Machne,
Murray & Stadler
2017](http://www.nature.com/articles/s41598-017-12401-8) for analysis
of segmentations by [Karl](https://github.com/raim/segmenTier).

#### Genomic Coordinate Indexing

Accessing genomic coordinates efficiently by indexing, used by
[u'r gene bro](https://gitlab.com/raim/genomeBrowser) and
[Karl](https://github.com/raim/segmenTier).

#### Positional Alignment and DNA Structural Patterns

Align genomic intervals around specific genomic sites, such as transcription
start sites, and calculate position-specific statistics. E.g. to generate
sequence or DNA motif enrichment, or average DNA binding data profiles.

#### Sequence Spectra 

... coming soon

Analyzing periodic enrichment of oligomers and DNA structural
parameters, after [Lehmann, Machne & Herzel
2014](https://doi.org/10.1093/nar/gku641).


### General Utilities

#### Parsers

* `parseGEOSoft` parses GEO Soft family files of microarray data sets
into data matrices, and accompanying probe-ID mapping, and sample/data
annotation
* `summarizeGEOSoft` offers a light-weight summarization function, to
average probe data for features with multiple probes 
* `gff2tab` parses a GFF file into tabular format, including
collection of attributes into data columns

## TODO

* Vignettes:
    - clusterTools v segmenTools,
    - time series clustering v genomic intervals,
    - used in: genomeBrowser; input from: segmenTier, dpseg,
    - class `clusterOverlaps`: sort and plot overlap enrichment profiles,
      produced by clusterCluster, clusterAnnotation, clusterProfile, 
      segmentOverlaps.

* clusterCluster: add fields for statistical corrections,
* clusterTimeseries: re-cluster cset by kmeans with centers initialized by flowclust cluster centers,
* segmenTier: change indexing to allow segments of length 1, see dpseg.
