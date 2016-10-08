#!/usr/bin/env Rscript
## Kiley Graim
## For generating HOCUS patient-patient networks

## Load libraries
library(getopt)#, quietly=T)

## Calculate similarity between 2 patient vectors- Hamming
get.sim <- function( a, b ) {
  if( length(a) != length(b) ) {
    cat("ERROR: Vectors not of same length")
    return(-1)
  }
  return( sum(a==b)/length(a) )
}

# Calculate similarity between 2 patient vectors - Jaccard
#get.sim <- function( a, b ) {
#  if( length(a) != length(b) ) {
#    cat("ERROR: Vectors not of same length")
#    return(-1)
#  }
#  ## Get number of disagreeing voxels
#  ##  ( and B) / (A or B)
#  val <- sum((a==b)&(a==1)) / (sum(a==1)+sum(b==1))
#
#  return( val )
#}

# Calculate similarity between 2 patient vectors - Pearson 
#get.sim <- function( a, b ) {
#  if( length(a) != length(b) ) {
#    cat("ERROR: Vectors not of same length")
#    return(-1)
#  }
#  return( cor(a,b,method='pearson') )
#} 

# Calculate similarity between 2 patient vectors - TFIDF 
#get.sim <- function( a, b ) {
#  if( length(a) != length(b) ) {
#    cat("ERROR: Vectors not of same length")
#    return(-1)
#  }
#  
#  return( sum(a[a==b])/length(a) )
#} 


## Main
main <- function(argv) {
  if (missing(argv)) argv <- commandArgs(trailingOnly=T)

  ## Parse options and set defaults
  spec <- matrix(c(    'help',     'h', 0, "logical", "Usage information",
    'input', 'i', 1, "character", "An input file containing the feature data, where rows are samples.",
    'delim', 'd', 1, "character", "Delimiter, default is tab.",
    'header', 'H', 0, "logical", "Indicates the file has a header. Default False.",
    'order', 'r', 0, "integer", "Indicates the file has a header. Default False.",
    'prefix', 'o', 1, 'character', 'Output file prefix'
  ),ncol=5,byrow=TRUE);
  o = getopt(spec)

  ## Print help and exit nicely, if asked
  if ( !is.null(o$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  }

  ## Set defaults, exit if input file is missing.
  if( is.null(o$prefix) ) { o$prefix = 'hocus_network' }
  if( is.null(o$H) ) { o$H=F }
  if( is.null(o$d) ) { o$d='\t' }
  if( is.null(o$order) ) { o$order=2 }
  if( is.null(o$input) ) { 
    cat(getopt(spec, usage=T))
    stop("Input file required")
  }

  ## Make the output file name
  out.prefix <- paste( o$prefix, 'tab', sep='.' )

  ## Load data
  cat( paste("Opening",o$input,"for reading...\n",sep=' ') );flush.console()
  if(!is.null(o$H)) { 
    dat.voxels <- data.matrix( read.table(o$input, sep=o$d, header=o$H, row.names=1, check.names=F) )
    dat.voxels <- dat.voxels[-1,] 
  } else {
    dat.voxels <- data.matrix( read.table(o$input, sep=o$d, header=F, row.names=1, check.names=F) )
  }
  if (!is.numeric(dat.voxels)) {
    stop("Input file isn't numeric. Check for strings after the first column.")
  }

  ## Initialize the patient-patient network
  hocus.net <- matrix(NA, nrow=nrow(dat.voxels), ncol=nrow(dat.voxels))
  rownames(hocus.net) <- rownames(dat.voxels)
  colnames(hocus.net) <- rownames(dat.voxels)

  ## Calculate pairwise patient similarities
  cat( "Calculating patient-patient similarities for all patients" ); flush.console()
  for(i in 1:nrow(hocus.net)) {
    hocus.net[i,] <- apply(dat.voxels, 1, function(x) { get.sim(dat.voxels[i,],x) } )
    cat('.')
  }

  ## Raise to the X order
  if(!(o$order <3)){
    for(i in 3:o$order) {
      hocus.net <- hocus.net %*% hocus.net
    }
    hocus.net <- hocus.net/max(hocus.net)
  }

  ## Write patient-patient network to specified output file
  cat( paste("\nWriting patient-patient network to", out.prefix, "\n") );flush.console()
  write.table(hocus.net, file=out.prefix, sep='\t', row.names=T, col.names=T, quote=F)
}

if (Sys.getenv("RGETOPT_DEBUG") != "") {debug(main)}
main()
