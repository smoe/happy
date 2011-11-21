#####################################################################
#
# happyrqtl.R
#
# copyright (c) 2002-2010,  Richard Mott
# Modified by:              Danny Arends
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
# Part of the happy package and iqtl package
# Contains:
#   -  load.condensed.database (version september 2010)
#           loads a pre-created database of condensed HAPPY design matrices into memory
#           returns an object of class "condensed.happy" with the following fields
#           subjects
#           strains
#             chr$additive$$chr$chr1 ....chr$additive$chr19 ....
#             additive$condensed.summary
#             chr$full$chr$chr1 ....chr$full$chr19 ....
######################################################################

load.condensed.database <- function( dir="./CONDENSED.20/" ) {

  load( paste(dir, "/db.RData",sep=""))

  load( paste(dir, "/strains.RData", sep=""))
  db$strains = strains;
  
  load( paste(dir, "/subjects.RData", sep=""))
  db$subjects = subjects;

  cat( "loading condensed db", dir, "\n")

  models = db$models 
  for( model in models ) {
    n = 0;
    db[[model]] = list()
    db[[model]]$chr = list()
    db[[model]]$matrices = list()
    h = hash()
    files = paste(dir, "/chr/", model, "/", db$chrs, ".RData", sep="")
    for ( i in 1:length(files)) {
      chr.data = files[i]
      load( chr.data )
      db[[model]]$chr[[db$chrs[i]]]=condensed
      db[[model]]$matrices = c( db[[model]]$matrices, condensed) 
      n = n + length(db[[model]]$chr[[db$chrs[i]]])
      lapply( condensed, function( item ) { h[[item$from.marker]] <<- item$mat } )
    }
    db[[model]]$hash = h
    cat( "model", model, "read ", n, " matrices\n");
    cat( "loading genome summary\n")
    db[[model]]$summary = NULL
    files = paste(dir, "/chr/", model, "/", db$chrs, ".summary.RData", sep="")
    for (chr.summary in files ) {
      load( chr.summary )
      db[[model]]$summary = rbind( db[[model]]$summary, genome )
    }

    cat ( "loading condensed summary\n")
    db[[model]]$condensed.summary = NULL
    files = paste(dir, "/chr/", model, "/", db$chrs, ".condensed.summary.RData", sep="")
    for (condensed.chr.summary in files ) {
      load( condensed.chr.summary )
      db[[model]]$condensed.summary = rbind( db[[model]]$condensed.summary, condensed.summary )
    }

# create some data objects so that the resulting object looks like a genome.cache...

    idx = match( db[[model]]$condensed.summary$from.marker, db[[model]]$summary$marker)
    db[[model]]$subjects = subjects
    db[[model]]$subjects = strains
    db[[model]]$genome = db[[model]]$summary[idx,]
    db[[model]]$genome$map = db[[model]]$genome$cm
    db[[model]]$genome$chr = gsub( "chr", "", db[[model]]$genome$chromosome )
    db[[model]]$markers = db[[model]]$genome$marker
    db[[model]]$chromosome = db[[model]]$genome$chromosome
    db[[model]]$map = db[[model]]$genome$map
  }
  load( paste(dir, "/strains.RData", sep=""))
  db$strains = strains;
  
  load( paste(dir, "/subjects.RData", sep=""))
  db$subjects = subjects;
  
  class(db) = "condensed.happy";
  return(db)
}		
	