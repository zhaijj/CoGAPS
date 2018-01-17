
#' patternMarkers
#'
#' @param Amatrix A matrix of genes by weights resulting from CoGAPS or other NMF decomposition
#' @param scaledPmatrix logical indicating whether the corresponding pattern matrix was fixed to have max 1 during decomposition
#' @param Pmatrix the corresponding Pmatrix (patterns X samples) for the provided Amatrix (genes x patterns). This must be supplied if scaledPmatrix is FALSE.
#' @param threshold # the type of threshold to be used. The default "all" will distribute genes into pattern with the lowest ranking. The "cut" thresholding by the first gene to have a lower ranking, i.e. better fit to, a pattern.
#' @param lp a vector of weights for each pattern to be used for finding markers. If NA markers for each pattern of the A matrix will be used.
#' @param full logical indicating whether to return the ranks of each gene for each pattern
#'
#' @return By default a non-overlapping list of genes associated with each \code{lp}. If \code{full=TRUE} a data.frame of
#' genes rankings with a column for each \code{lp} will also be returned.
#' @export
#'
#' @examples \dontrun{
#' patternMarkers(Amatrix=AP$Amean,scaledPmatrix=FALSE,Pmatrix=NA,threshold="All",full=TRUE)
#' }
#'
patternMarkers <- function(
    Amatrix=NA, #A matrix of genes by weights resulting from CoGAPS or other NMF decomposition
    scaledPmatrix=FALSE, # logical indicating whether the corresponding pattern matrix was fixed to have max 1 during decomposition
    Pmatrix=NA, #the corresponding Pmatrix (patterns X samples) for the provided Amatrix (genes x patterns). This must be supplied if scaledPmatrix is FALSE.
    threshold="all", # the type of threshold to be used. The default "All" will distribute genes into pattern with the highest ranking. 
    lp=NA, # a vector of weights for each pattern to be used for finding markers. If NA markers for each pattern of the A matrix will be used.
    full=FALSE #logical indicating whether to return the ranks of each gene for each pattern.
){

  dnames <- dimnames(Amatrix)
    if(scaledPmatrix==FALSE){
        if(!is.na(Pmatrix)){
            pscale <- apply(Pmatrix,1,max)   # rescale p's to have max 1
            Amatrix <- sweep(Amatrix, 2, pscale, FUN="*")   # rescale A in accordance with p's having max 1
        }
        else(warning("P values must be provided if not already scaled"))
    }
    # find the A with the highest magnitude
  dimnames(Amatrix) <- dnames
  Arowmax <- t(apply(Amatrix, 1, function(x){
    ## this prevents NA values from creeping due to rows with
    ## zero means
    if (mean(x) == 0){ return(rep(0, times=length(x))) }
    else { return(x/max(x)) }
  }))
    pmax<-apply(Amatrix, 1, max)
    # determine which genes are most associated with each pattern
    sstat<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=dimnames(Amatrix))#list()
    ssranks<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=dimnames(Amatrix))#lst
    ssgenes<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=NULL)
    nP=dim(Amatrix)[2]
    if(!is.na(lp)){
        if(length(lp)!=dim(Amatrix)[2]){
            warning("lp length must equal the number of columns of the Amatrix")
        }
        ## this I is missing a loop
        for (i in 1:nP){
          sstat[,i] <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
          ssranks[order(sstat[,i]),i] <- 1:length(sstat[,i])
          ssgenes[,i]<-names(sort(sstat[,i],decreasing=FALSE))
        }
    } else {
        for(i in 1:nP){
          lp <- rep(0,dim(Amatrix)[2])
          lp[i] <- 1
          sstat[,i] <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
          ssranks[order(sstat[,i]),i] <- 1:length(sstat[,i])
          ssgenes[,i]<-names(sort(sstat[,i],decreasing=FALSE))
        }
    }
    if(threshold=="cut"){
        geneThresh <- sapply(1:nP,function(x) min(which(ssranks[ssgenes[,x],x] >
                                                        apply(ssranks[ssgenes[,x],],1,min))))
        ssgenes.th <- sapply(1:nP,function(x) ssgenes[1:geneThresh[x],x])
        ## geneThresh <- apply(sweep(ssranks,1,t(apply(ssranks, 1, min)),"-"),2,
        ## function(x) which(x==0))
        #ssgenes.th <- lapply(geneThresh,names)
    }
    if(threshold=="unique"){
        pIndx<-apply(sstat,1,which.min) 
        gBYp <- lapply(sort(unique(pIndx)),function(x) names(pIndx[pIndx==x]))
        ssgenes.th <- lapply(1:nP, function(x){ ## 
            if (x > length(gBYp)){ return(NA) }
            ssgenes[which(ssgenes[,x] %in% gBYp[[x]]),x]
        })
    }
    if(threshold=="all"){
        ssgenes.th <- ssgenes
    }
    if(full){
        return(list("PatternMarkers"=ssgenes.th,"PatternRanks"=ssranks,"PatternScores"=sstat))
    } else{return("PatternMarkers"=ssgenes.th)}
}

