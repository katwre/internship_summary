#'---
#'title: "Summary of new features of genomation"
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


#' [Genomation](https://bioconductor.org/packages/devel/bioc/html/genomation.html) is an R package to summarize, annotate
#' and visualize genomic intervals. It contains a collection of tools for visualizing and analyzing genome-wide data sets,
#' i.e. RNA-seq, reduced representation bisulfite sequencing (RRBS) or chromatin-immunoprecipitation followed by sequencing 
#' (Chip-seq) data.
#'
#' Recently were added new fetures to genomation and here we present them on example of 
#' binding profiles of 6 transcription factors around the Ctcf binding sites derived from Chip-seq.
#' All new functionalities are available in the latest version of genomation that can be found on 
#' [it's github website](https://github.com/BIMSBbioinfo/genomation).
#'
# install the package from github
library(devtools)
install_github("BIMSBbioinfo/genomation",build_vignettes=FALSE)
library(genomation)

#' # Extending genomation to work with paired-end BAM files
#'
#' Genomation can work with paired-end BAM files. Mates from reads
#' are treated as fragments (they are stitched together).

genomationDataPath = system.file('extdata',package='genomationData')
bam.files = list.files(genomationDataPath, full.names=TRUE, pattern='bam$')
bam.files = bam.files[!grepl('Cage', bam.files)]

#' # Accelerate functions responsible for reading genomic files
#' This is achived by using _readr::read_delim_ function to read genomic files
#' instead of _read.table_.
#' Additionally if skip="auto" argument is provided in _readGeneric_or track.line="auto" in other functions 
#' that read genomic files, e.g. _readBroadPeak_ then
#' UCSC header is detected (and first track).
library(GenomicRanges)

ctcf.peaks = readBroadPeak(file.path(genomationDataPath, 
                                     'wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz'))
ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == 'chr21']
ctcf.peaks = ctcf.peaks[order(-ctcf.peaks$signalValue)]
ctcf.peaks = resize(ctcf.peaks, width=1000, fix='center')

#' # Parallelizing data processing in ScoreMatrixList
#' We use _ScoreMatrixList_ function to extract coverage values of all transcription factors 
#' around ChIP-seq peaks. _ScoreMatrixList_ was improved by adding new argument _cores_
#' that indicates number of cores to be used at the same time (by using _parallel:mclapply_).

sml = ScoreMatrixList(bam.files, ctcf.peaks, bin.num=50, type='bam', cores=2)

# file descriptions of transcription factors
sampleInfo = read.table(system.file('extdata/SamplesInfo.txt',
                                    package='genomationData'),header=TRUE, sep='\t')
names(sml) = sampleInfo$sampleName[match(names(sml),sampleInfo$fileName)]

#' # Arithmetic, indicator and logic operations as well as subsetting work on score matrices
#' Arithmetic, indicator and logic operations work on _ScoreMatrix_, _ScoreMatrixBin_ and _ScoreMatrixList_<br />
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


#' # New arguments and improvements in visualization functions

#' Due to large signal scale of rows of each element in the _ScoreMatrixList_ 
#' we scale them.
sml.scaled = scaleScoreMatrixList(sml)

#' # Faster heatmaps
#'
#'  _HeatMatrix_ and _multiHeatMatrix_ function works faster by faster assigning colors.

#' Heatmap profile of scaled coverage shows a colocalization of Ctcf, Rad21 and Znf143. 
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500))

#' ## New clustering possibilities in heatmaps: "clustfun" argument in multiHeatMatrix
#' _clustfun_ allow to add more clustering functions and integrate them with the heatmap function
#' multiHeatMatrix. It has to be a function that returns a vector
#' of integers indicating the cluster to which each point is allocated.
#' Previous version of _multiHeatMatrix_ 
#' could cluster rows of heatmaps using only k-means algorithm.

# k-means algorithm with 2 clusters
cl1 <- function(x) kmeans(x, centers=2)$cluster
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1)

# hierarchical clustering with Ward's method for agglomeration into 2 clusters
cl2 <- function(x) cutree(hclust(dist(x), method="ward"), k=2)
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl2)

#' ## Defining which matrices are used for clustering: "clust.matrix" in multiHeatMatrix
#' clust.matrix argument indicates which matrices are used for clustering.
#' It can be a numerical vector of indexes of matrices or a character vector of
#' names of the ‘ScoreMatrix’ objects in 'ScoreMatrixList'. 
#' Matrices that are not in clust.matrix are ordered 
#' according to the result of the clustering algorithm.
#' By default all matrices are clustered.

multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1, clust.matrix = 1)

#' ## Central tendencies in line plots: centralTend in plotMeta
#' We extended visualization capabilities for meta-plots.
#' _plotMeta_ function can plot not only mean, but also median as central tendency
#' and it can be set up using _centralTend_ argument.
#' Previously user could plot only mean.

plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),
	 xcoords=c(-500, 500),
	 winsorize=c(0,99),
	 centralTend="mean")

#' ## Smoothing central tendency: smoothfun in plotMeta
#' We added _smoothfun_ argument to smooth central tendency as well as dispersion bands around
#' it which is shown in the next figure. Smoothfun has to be a function that returns a 
#' list that contains a vector of y coordinates (list with numeric vector named '$y').

plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),
	 xcoords=c(-500, 500),
	 winsorize=c(0,99),
	 centralTend="mean",  
	 smoothfun=function(x) stats::smooth.spline(x, spar=0.5))

#' ## Plotting dispersion around central lines in line plots: dispersion in plotMeta	 
#' We added new argument _dispersion_ to plotMeta that shows dispersion bands around _centralTend_.
#' It can take one of the arguments:
#'          
#'* "se"  shows standard error of the mean and 95 percent
#'              confidence interval for the mean
#'* "sd"  shows standard deviation and 2*(standard deviation)
#'* "IQR" shows 1st and 3rd quartile and confidence interval
#'              around the median based on the median +/- 1.57 *
#'              IQR/sqrt(n) (notches)

plotMeta(mat=sml, profile.names=names(sml),
	 xcoords=c(-500, 500),
	 winsorize=c(0,99),
	 centralTend="mean",  
	 smoothfun=function(x) stats::smooth.spline(x, spar=0.5),
	 dispersion="se", lwd=4)
  
            	    
#' # Calculating scores that correspond to k-mer or PWM matrix occurence: patternMatrix object 
#' We added new function _patternMatrix_ that calculates
#' k-mer and PWM occurrences over predefined equal width windows.
#' If one pattern (character of length 1 or PWM matrix) is given then it returns ScoreMatrix, 
#' if more than one then ScoreMatrixList.
#' It finds either positions of pattern hits above a specified threshold and creates score matrix filled with
#' 1 (presence of pattern) and 0 (its absence) or
#' matrix with score themselves.
#' Windows can be a DNAStringList object or GRanges object (but then genome argument has to be provided,
#' a BSgenome object).
#' It is still under development, but it can be installed using:

install_github("katwre/genomation",ref="patternMatrix",build_vignettes=FALSE)	    

#ctcf motif from the JASPAR database
ctcf.pfm = matrix( as.integer(c(87, 167, 281,  56,   8, 744,  40, 107 ,851  , 5 ,333 , 54 , 12,  56, 104, 372 , 82, 117 ,402, 
		     291, 145 , 49, 800 ,903,  13, 528, 433 , 11 ,  0 ,  3 , 12,   0 ,  8, 733 , 13, 482 ,322, 181, 
		     76 ,414 ,449  ,21 ,  0 , 65 ,334 , 48 , 32, 903, 566, 504 ,890 ,775  , 5 ,507 ,307 , 73, 266, 
		     459 ,187, 134  ,36,   2 , 91,11, 324 , 18,   3 ,  9 ,341 ,  8 , 71 , 67 , 17 , 37, 396,  59 )), 
		     ncol=19,byrow=TRUE)
rownames(ctcf.pfm) <- c("A","C","G","T")

library(Biostrings)
prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25)
priorProbs = prior.params/sum(prior.params)
postProbs = t( t(ctcf.pfm + prior.params)/(colSums(x)+sum(prior.params)) )
ctcf.pwm = unitScale(log2(postProbs/priorProbs))

library(BSgenome.Hsapiens.UCSC.hg19)
hg19 = BSgenome.Hsapiens.UCSC.hg19

p = patternMatrix(pattern=ctcf.pwm, windows=ctcf.peaks, genome=hg19, min.score=0.8)

#' Visualization of the patternMatrix
#' _patternMatrix_ (here as ScoreMatrix object) can be visualized using i.e. heatMatrix, heatMeta or plotMeta functions.
heatMatrix(p, xcoords=c(-500, 500), main="CTCF motif")

plotMeta(mat=p, xcoords=c(-500, 500), smoothfun=function(x) stats::lowess(x, f = 1/10), 
         line.col="red", main="ctcf motif")

#' # Integration with Travis CI for auto-testing
#' Recently we integrated genomation with [Travis CI](travis-ci.org). It allows users to see current status
#' of the package which is updated during every change of the package. Travis
#' automatically runs R CMD CHECK and reports it. Shields shown below are on the genomation github site:<br />
#' [https://github.com/BIMSBbioinfo/genomation](https://github.com/BIMSBbioinfo/genomation)
#' <br />
#' Status [![Build Status](https://api.travis-ci.org/BIMSBbioinfo/genomation.svg)](https://travis-ci.org/BIMSBbioinfo/genomation)   [![codecov.io](https://codecov.io/github/BIMSBbioinfo/genomation/coverage.svg)](https://codecov.io/github/BIMSBbioinfo/genomation?branch=master)     [![BioC_years](http://www.bioconductor.org/shields/years-in-bioc/genomation.svg)](http://www.bioconductor.org/packages/release/bioc/html/genomation.html)     [![BioC_availability](http://www.bioconductor.org/shields/availability/release/genomation.svg)](http://www.bioconductor.org/packages/release/bioc/html/genomation.html)
#' <br />
#' <br />
# <br />
sessionInfo()
 
