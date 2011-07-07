#####################################################################
#
# happyrqtl.R
#
# copyright (c) 2010, Danny Arends
# last modified sep, 2010
# first written Sep, 2010
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the iqtl package
# Contains: 
#   -  loadhappyobject
#         Wrapper for loading the happydatabase object it, returns object of type condensed.happy
#   -  probstomostlikely
#         Convert a genotype probability list for a marker to genotype with the highest probability
#   -  convertprobtogenotype
#         Convert a genotype probability matrix for a single individual to the individuals genotype with the highest probability
#   -  getchromosomesplits
#         Gets ids of the markers at which a new chromosome starts
#   -  convertallprobabilities
#         Converts all pobability matrices into a genotypematrix
#   -  converthappytocross
#         Convert a happy object to a cross object
#   -  doconversion
#         Wrapper for the entire conversion from happy to cross (sourcing external scripts, loading libs)
#   -  dotest
#         Tests some basic functionality of R/qtl on the new cross object
#
######################################################################

loadhappyobject <- function(datalocation="d:/happydata/condensed.x"){
  a <- load.condensed.database(datalocation)
  a
}

probstomostlikely <- function(g){
  as.numeric(which(g == max(g)))
}

convertprobtogenotype <- function(genotypeprobs){
  a <- apply(genotypeprobs,1,probstomostlikely)
  r <- NULL
  for(x in 1:length(a)){
    r <- c(r,a[[x]][1])
  }
  r
}

getchromosomesplits <- function(locations){
  l <- 0
  r <- NULL
  for(x in 1:length(locations)){
    if(!(l <= locations[x])){
      r <- c(r,x)
    }
    l <- locations[x]
  }
  r
}

convertallprobabilities <- function(happyobject){
  #Converting happy coding to R/qtl
  fullgenotypes <- NULL
  nmarkers <- length(happyobject$full[[2]])
  cat("Start converting ",nmarkers," marker codings\n")
  for(x in 1:nmarkers){
    if(x%%1000==0)cat("Conversion done for",x,"/",nmarkers,"\n")
    fullgenotypes <- cbind(fullgenotypes,convertprobtogenotype(happyobject$full[[2]][[x]]$mat))
  }
  cat("Conversion complete\n")
  fullgenotypes
}

converthappytocross <- function(happyobject){
  
  fullgenotypes <- convertallprobabilities(happyobject)

  #chreate new objects
  nchr <- length(happyobject$chrs)
  ntraits <- 3
  nindividuals <- nrow(happyobject$full[[2]][[1]]$mat)
  locations <- happyobject$full$map
  genotypes <- vector("list", nchr)

  chr <- 1
  prev <- 1
  cat("Creating genotype part of Cross object\n")
  for(x in c(getchromosomesplits(locations),length(locations))){
    cat("Chromosome:",chr,"from",prev,"to",x,"\n")
    genotypes[[chr]]$data <- fullgenotypes[,prev:(x-1)]
    colnames(genotypes[[chr]]$data) <- happyobject$full$markers[prev:(x-1)]
    
    genotypes[[chr]]$map <- locations[prev:(x-1)]
    names(genotypes[[chr]]$map) <- happyobject$full$markers[prev:(x-1)]
    
    if(happyobject$chrs[chr]=="chrX"){
      class(genotypes[[chr]]) <- "X"
    }else{
      class(genotypes[[chr]]) <- "A"
    }
    
    chr <- chr+1
    prev <- x
  }
  names(genotypes) <- happyobject$chrs
  cat("Creating Cross object\n")
  cross <- NULL
  cross$geno <- genotypes
  cross$pheno <- as.data.frame(matrix(runif(nindividuals*ntraits),nindividuals,ntraits))
  #Hide 2 QTL in phenotype 1
  cross$pheno[,1] <- genotypes[[5]]$data[,5] + genotypes[[1]]$data[,200] + cross$pheno[,1]
  class(cross) <- c("f2","cross")

  cat("Cleaning up Genotype intermediate\n")
  rm(happyobject)
  rm(fullgenotypes)
  r <- gc()
  cat("Done\n")
  cross
}

doconversion <- function(sourcedir="d:/happy",happydata="d:/happydata/condensed.x",output="cross.rdata"){
  library(hash)
  library(qtl)
  setwd(sourcedir)
  source("hackinghappy.R")
  happyobject <- loadhappyobject(happydata)
  cross <- converthappytocross(happyobject)
  #Save the cross file
  save(cross,file=output)
  cat("Please close R to cleanup the objects created\n")
  #close R to cleanup the objects created ( the rm(object) + gc() function isn't propperly woking under windows2008 RC2 )
  cross
}

dotest <- function(sourcedir="d:/happy",input="cross.rdata"){
  #Restart R
  library(qtl)
  setwd(sourcedir)
  load(input)
  #Test some basic functions (These need to be updated to support HS and a coding number class(cross) <- c("hs",8,"cross"))
  summary(cross)
  nf <- layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), TRUE)
  layout.show(nf)
  plot.map(cross)
  plot.info(cross)
  #Try to find our QTL (Not working ofcourse because most genotypes are treated as missing )
  r <- scanone(cross)
  plot(r)
  cross
}
