#' Create table with ligand-receptor pairs between cells.
#' This function combines gene expression output data from SEURAT with
#' a database of Ligand-Receptor pairs to identify potential
#' cell-cell interactions in a given dataset.
#' @import dplyr circlize reshape2
#' @param seuratDEGS Exact output from Seurat::FindAllMarkers function. Must contain a 'gene' column
#' @param LRdatabase Table of ligand-receptor pairs
#' @param subsetgenes Vector containing a subset of genes
#' @param from vector of cell type names to subset outgoing interactions. Default is all cells.
#' @param to vector of cell type names to subset incoming interactions. Default is all cells.
#' @export
LigandReceptorPairsTable <- function(seuratDEGS, LRdatabase, subsetgenes = seuratDEGS$gene, from = celltypelabels, to = celltypelabels){
  
  ncells <- as.numeric(length(unique(seuratDEGS$cluster))) #number of distinct cell types in SEURAT object
  celltypelabels <- as.character(unique(seuratDEGS$cluster)) #save names as vector of strings
  #LRdatabase:Must have 5 colums in this order: 'Pair', "Ligand", "Ligand.name", "Receptor" and "Receptor.name"
  #subsetgenes: subset of genes to create ligand-receptor pairs matix with. Must be a vector.
  
  # create a matrix with the fold change values from SEURAT output using reshape library
  foldchanges <- reshape2::dcast(seuratDEGS,formula = gene~cluster,fun.aggregate = sum,value.var = "avg_logFC") #this function creates the matrix
  FC.receptors <- merge(foldchanges, LRdatabase, by.x = "gene", by.y = "Receptor") #combine with ligand receptor database to create matrix of receptors
  
  potential.Pairs <- merge(FC.receptors,foldchanges, by.x = "Ligand", by.y = "gene", no.dups = F) # add ligand information to matrix to identify pairs
  names(potential.Pairs)[2] <- "Receptor" # correct column name
  
  
  ## create subset of potential pairs with desired genes
  potential.Pairs <- potential.Pairs[potential.Pairs$Receptor %in% subsetgenes | potential.Pairs$Ligand %in% subsetgenes, ]
  
  counter = replicate(ncells,0) #initialize a vector of zeroes with length=ncells
  number.of.pairs <- matrix(ncol = ncells, # initialize matrix of dimensions ncells x ncells to store number of potential pairs
                            nrow = ncells,
                            dimnames = list(colnames(potential.Pairs[3:(2+ncells)]),
                                            colnames(potential.Pairs[3:(2+ncells)])))
  
  
  # create matrix with number of potential pairs between cell types
  for (j in 1:ncells){
    counter = replicate(ncells, 0) #initialize a vector of zeroes with length=ncells
    for (i in 1:NROW(potential.Pairs)){
      if(potential.Pairs[i,j+2]>0){ #the first two columns should contain Ligand and Receptor genes
        for(k in 1:ncells){
          if(potential.Pairs[i,k+(5+ncells)]>0){ #number 5 is from the number of columns in LR Database which is maintained in the pair matrix
            counter[k] <- counter[k] + 1
          }
        }
      }
    }
    number.of.pairs[j,] <- counter[1:ncells]
  }
  
  number.of.pairs <- t(number.of.pairs)
  # Transposes matrix with 'Number of Pairs' is organized with recipient (receptor-expressing) cells in columns and provider (ligand-expressing) cells in rows
  colnames(number.of.pairs) <-  celltypelabels # fix colnames to remove ".x /.y"
  rownames(number.of.pairs) <- celltypelabels # fix colnames to remove ".x /.y"
  
  melted<- reshape2::melt(data = number.of.pairs, varnames = c("from", "to")) #use melt function from reshape2 library to create a table from matrix
  
  ### Add ligand-receptor pair labels to melted table:
  melted$pairs <- NA #add a "pairs" column to add information to
  ncols = NCOL(number.of.pairs) + NCOL(LRdatabase) # number of columns in first half of potential pair table
  
  for (i in 1:ncells){
    
    cell1 <- celltypelabels[i]
    
    temp <- potential.Pairs[potential.Pairs[,i+ncols] >0,c(1:ncols, i+ncols)]  #Create a subset from the matrix table that contains all the rows that have
    # ligand values >0 in column i + ncols (tested ligand), and all columns with gene labels and receptor expression values (columns 1:ncol)
    
    colnames(temp)[3:(2+ncells)] <- celltypelabels
    
    temp$Pair <- as.character(temp$Pair) # make sure column with pair information is of character type
    
    
    if(NROW(temp)>0){ # Look for values > 0 in the expression matrix
      
      for (k in 1:ncells) {
        pair <- c()
        cell2 <- colnames(temp)[2+k] # cell type potentially expressing receptors for the tested ligand
        for (j in 1:NROW(temp)) {
          if(temp[j,k+2]>0 ){
            pair <- c(temp$Pair[j], pair) # save pair value in vector
          }
        }
        pair <- paste(pair, collapse=", ") # convert vector to a string of characters
        melted[melted$from %in% cell1 & melted$to %in% cell2, ]$pairs <- pair # save pair information to melted table
      }
    }
  }
  melted <- melted[melted$from %in% from & melted$to %in% to, ]
  return(melted)
  
}




#' Chord Plot of ligand-receptor interactions
#' This function combines gene expression output data from SEURAT with
#' a database of Ligand-Receptor pairs to identify potential
#' cell-cell interactions in a given dataset and uses the circlize package
#' to produce a chord plot representative of those interactions
#'
#' @import dplyr circlize reshape2
#' @param cellcolors vector of color names or codes for each cell type
#' @param seuratDEGS Exact output from Seurat::FindAllMarkers function. Must contain a 'gene' column
#' @param LRdatabase Table of ligand-receptor pairs
#' @param subsetgenes Vector containing a subset of genes
#' @param from vector of cell type names to subset outgoing interactions. Default is all cells.
#' @param to vector of cell type names to subset incoming interactions. Default is all cells.
#' @export
PairsPlot <- function(seuratDEGS, LRdatabase, cellcolors, subsetgenes=seuratDEGS$gene, from = celltypelabels, to = celltypelabels){
  #ncells: number of distinct cell types in SEURAT object
  #seuratDEGs: direct output from the SEURAT "FindAllMarkers()" function
  #LRdatabase: Table of ligands and receptor pairs. Must have 5 colums in this order: 'Pair', "Ligand", "Ligand.name", "Receptor" and "Receptor.name"
  #subsetgenes: subset of genes to create ligand-receptor pairs matix with
  # create a matrix with the fold change values from SEURAT output using reshape library
  ncells <- as.numeric(length(unique(seuratDEGS$cluster))) #number of distinct cell types in SEURAT object
  celltypelabels <- as.character(unique(seuratDEGS$cluster))
  
  if(length(cellcolors) == ncells){
    
    foldchanges <- reshape2::dcast(seuratDEGS,formula = gene~cluster,fun.aggregate = sum,value.var = "avg_logFC") #this function creates the matrix
    FC.receptors <- merge(foldchanges, LRdatabase, by.x = "gene", by.y = "Receptor") #combine with ligand receptor database to create matrix of receptors
    
    potential.Pairs <- merge(FC.receptors,foldchanges, by.x = "Ligand", by.y = "gene", no.dups = F) # add ligand information to matrix to identify pairs
    names(potential.Pairs)[2] <- "Receptor" # correct column name
    
    ## create subset of potential pairs with desired genes
    potential.Pairs <- potential.Pairs[potential.Pairs$Receptor %in% subsetgenes | potential.Pairs$Ligand %in% subsetgenes, ]
    
    counter = replicate(ncells,0) #initialize a vector of zeroes with length=ncells
    number.of.pairs <- matrix(ncol = ncells, # initialize matrix of dimensions ncells x ncells to store number of potential pairs
                              nrow = ncells,
                              dimnames = list(colnames(potential.Pairs[3:(2+ncells)]),
                                              colnames(potential.Pairs[3:(2+ncells)])))
    
    # create matrix with number of potential pairs between cell types
    for (j in 1:ncells){
      counter = replicate(ncells, 0) #initialize a vector of zeroes with length=ncells
      for (i in 1:NROW(potential.Pairs)){
        if(potential.Pairs[i,j+2]> 0){ #the first two columns should contain Ligand and Receptor genes
          for(k in 1:ncells){
            if(potential.Pairs[i,k+(5+ncells)]>0){ #number 5 is from the number of columns in LR Database which is maintained in the pair matrix
              counter[k] <- counter[k] + 1
            }
          }
        }
      }
      number.of.pairs[j,] <- counter[1:ncells]
    }
    
    number.of.pairs <- t(number.of.pairs)
    # Transposes matrix with 'Number of Pairs' is organized with recipient (receptor-expressing) cells in columns and provider (ligand-expressing) cells in rows
    colnames(number.of.pairs) <-  celltypelabels # fix colnames to remove ".x /.y"
    rownames(number.of.pairs) <- celltypelabels # fix colnames to remove ".x /.y"
    
    
    melted<- reshape2::melt(data = number.of.pairs, varnames = c("from", "to")) #use melt function from reshape2 library to create a table from matrix
    
    circlize::circos.par(gap.degree=5, gap.after=5)
    circlize::chordDiagram(melted[melted$from %in% from & melted$to %in% to, ], grid.col = cellcolors, link.lwd = 1, link.lty = 1, link.border = "black",
                           symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1,
                           link.arr.length = 0.1, link.arr.type = "big.arrow",
                           annotationTrack = c("grid"), link.largest.ontop = T, link.arr.lty = 1,grid.border = 1)
    circlize::circos.clear()
    
  }
  else {
    print("Number of colors does not mach number of cell types")
  }
}


