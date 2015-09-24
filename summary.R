#'---
#'title: "Summarize new features of genomation"
#'author: "Katarzyna Wreczycka"
#'output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    fig_width: 7
#'    fig_height: 7
#'    keep_md: true
#'    cache: true
#'---


#' [Genomation](https://github.com/BIMSBbioinfo/genomation) is an R package to summarize, annotate
#' and visualize genomic intervals. It contains a collection of tools for visualizing and analyzing genome-wide data sets,
#' i.e. RNA-seq, reduced representation bisulfite sequencing (RRBS) or chromatin-immunoprecipitation followed by sequencing 
#' (Chip-seq) data.
#'
#' We recently added new fetures to genomation and here are they presented on example of 
#' binding profiles of 6 transcription factors around the Ctcf binding sites.
#'
#' All new functionality are available on the latest version of genomation available on github.
#'
# install the package from github
# library(devtools)
# install_github("BIMSBbioinfo/genomation",build_vignettes=FALSE)
# library(genomation)
library(devtools)
load_all("../genomation")
library(GenomicRanges)


#' # Extending genomation to work with paired-end BAM files
#'
#' Genomation can work with paired-end BAM files. Mates from a pair
#' are treated as fragments (are stitched together).

genomationDataPath = system.file('extdata',package='genomationData')
bam.files = list.files(genomationDataPath, full.names=TRUE, pattern='bam$')
bam.files = bam.files[!grepl('Cage', bam.files)]

#' # Accelerate of function responsible for reading files
#' This is achived by using readr::read_delim function to read genomic files
#' instead of read.table.
#' Additionally if skip="auto" in readGeneric or track.line="auto" other functions
#' like readBroadPeak
#' then they detect UCSC header (and first track). #TODO

ctcf.peaks = readBroadPeak(file.path(genomationDataPath, 
                                     'wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz'))

ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == 'chr21']
ctcf.peaks = ctcf.peaks[order(-ctcf.peaks$signalValue)]
ctcf.peaks = resize(ctcf.peaks, width=1000, fix='center')

#' # Parallelizing data processing
#' We use ScoreMatrixList function to extract coverage values of all transcription factors 
#' around chipseq peaks. ScoreMatrixList was improved by adding new argument cores
#' that indicated number of cores to be used at the same time by using parallel:mclapply.

sml = ScoreMatrixList(bam.files, ctcf.peaks, bin.num=50, type='bam', cores=2)

# Names of .. we stored in the SamplesInfo.txt file.
sampleInfo = read.table(system.file('extdata/SamplesInfo.txt',
                                    package='genomationData'),header=TRUE, sep='\t')
names(sml) = sampleInfo$sampleName[match(names(sml),sampleInfo$fileName)]


#' # Arithmetic, indicator and logic operations as well as subsetting work on score matrices
#' Arithmetic, indicator and logic operations work on ScoreMatrix, ScoreMatrixBin and ScoreMatrixList<br />
#' objects, e.i.:<br />
#' Arith: "+", "-", "*", "^", "%%", "%/%", "/" <br />
#' Compare: "==", ">", "<", "!=", "<=", ">=" <br />
#' Logic: "&", "|"   <br />
sml1 = sml * 100
sml1

#' Subsetting:
sml[[6]] = sml[[1]]
sml 
sml[[6]] <- NULL


#' # New arguments in visualizing functions

#' Because of large signal scale the rows of each element in the ScoreMatrixList.
sml.scaled = scaleScoreMatrixList(sml)


#' Heatmap profile of scaled coverage shows a colocalization of Ctcf, Rad21 and Znf143. 
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500))


#' ## clustfun in multiHeatMatrix
#' clustfun allow to add more clustering functions and integrate them with heatmap function
#' multiHeatMatrix. clustfun argument should be a function that returns a vector
#' of integers indicating the cluster to which each point is allocated.
#' It's an extention of previous version that could cluster rows of heatmaps using k-means algorithm.

# k-means algorithm, 2 clusters
cl1 <- function(x) kmeans(x, centers=2)$cluster
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1)

# hierarchical clustering with Ward's method for agglomeration, 2 clusters
cl2 <- function(x) cutree(hclust(dist(x), method="ward"), k=2)
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl2)

#' ## clust.matrix in multiHeatMatrix
#' clust.matrix argument indicates which matrices are used for clustering.
#' It can be a numerical vector of indexes of matrices or a character vector of
#' names of the ‘ScoreMatrix’ objects in 'sml'. By default all matrices are clustered.

multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1, clust.matrix = 1)

#' ## centralTend in plotMeta

#' extending visualization capabilities for meta-plots - improving the plotMeta func-
#' tion to plot not only mean, but also median as a central tendency, adding possibil-
#' ity to plot dispersion bands around the central tendency and to smoothing central
#' tendency and dispersion bands,

plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),
	 xcoords=c(-500, 500),
	 winsorize=c(0,99),
	 centralTend="mean")

#' ## smoothfun in plotMeta
	 	 
plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),
	 xcoords=c(-500, 500),
	 winsorize=c(0,99),
	 centralTend="mean",  
	 smoothfun=function(x) stats::smooth.spline(x, spar=0.5))

#' ## dispersion in plotMeta	 
	 
plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),
	 xcoords=c(-500, 500),
	 winsorize=c(0,99),
	 centralTend="mean",  
	 smoothfun=function(x) stats::smooth.spline(x, spar=0.5),
	 dispersion="se", lwd=4)
	    
           
#' # Integration with Travis CI for auto-testing
#' Recently we integrated genomation with [Travis CI](travis-ci.org) which allows users to see current status
#' of the package which is updated during every change of the package. Travis
#' automatically runs R CMD CHECK and reports it. Such shields are visible on the genomation github site:<br />
#' [https://github.com/BIMSBbioinfo/genomation](https://github.com/BIMSBbioinfo/genomation)
#' <br />
#' Status [![Build Status](https://api.travis-ci.org/BIMSBbioinfo/genomation.svg)](https://travis-ci.org/BIMSBbioinfo/genomation)   [![codecov.io](https://codecov.io/github/BIMSBbioinfo/genomation/coverage.svg)](https://codecov.io/github/BIMSBbioinfo/genomation?branch=master)     [![BioC_years](http://www.bioconductor.org/shields/years-in-bioc/genomation.svg)](http://www.bioconductor.org/packages/release/bioc/html/genomation.html)     [![BioC_availability](http://www.bioconductor.org/shields/availability/release/genomation.svg)](http://www.bioconductor.org/packages/release/bioc/html/genomation.html)
#' <br />
#' <br />
# <br />
sessionInfo()
 