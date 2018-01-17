library(matrixStats)
library(parallel)
library(entropy)
library(Rcpp)

renameMat <- function(mat, ref.mat){
  rownames(mat) <- rownames(ref.mat)
  colnames(mat) <- colnames(ref.mat)
  mat
}

which.in <- function(gene, listP){
    inn <- sapply(listP, function(xx){ gene %in% xx })
    return(which(inn == TRUE))
}

assignGenes <- function(genes, lisst){
  assigns <- sapply(1:length(genes), function(ii){
    which.in(genes[ii], lisst)
  })
  names(assigns) <- genes
  assigns
}

pmarkToPat <- function(pmarkList, genes, ncores=3){
  do.call(rbind, mclapply(pmarkList, function(xx){
    assignGenes(genes, xx$PatternMarkers)
  }, mc.cores=ncores))
}

pMarkArray <- function(a.snap, p.snap, a.mean, p.mean, scaledP=FALSE, ncores=3){
  pmarkList <- mclapply(1:dim(a.snap)[3], function(ii){
    atmp <- renameMat(a.snap[,,ii], a.mean)
    ptmp <- renameMat(p.snap[,,ii], p.mean)
    patternMarkersC(A=atmp, P=ptmp)
  }, mc.cores=ncores)
  return(simplify2array(pmarkList))
}

## should take pattern array (genes x patterns x iters) and patternMarkers on mean mat
## return vector (#genes)
pVals <- function(patt.mat, ref.pats){
  sapply(1:nrow(patt.mat), function(ii){ mean(patt.mat[ii,] == ref.pats[ii]) })
}

scoreVars <- function(scores.arr){
  sapply(1:dim(scores.arr)[1], function(ii){ rowVars(scores.arr[ii,,]) })
}

pat.entropy <- function(pat.mat){ apply(pat.mat, 2, function(xx){ entropy(table(xx)) }) }

heatPlot <- function(mat, name){
  Heatmap(mat, cluster_rows=FALSE, cluster_columns=FALSE, name=name,
          show_row_names=FALSE, show_column_names=FALSE,
          row_names_side="left")
}

PUMP <- function(asnaps, psnaps, amean, pmean, scaledP=FALSE){
  mean.pats <- patternMarkers(A=amean, P=pmean, scaledPmatrix=scaledP)
  snap.pat <- pMarkArray(asnaps, psnaps, amean, pmean, scaledP)
  ps <- pVals(snap.pat, mean.pats)
  return(cbind(mean.pats+1, ps))
}

