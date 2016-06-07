#' Returns a vector of seeds based on the number of chains and initial seed
#'
#' @param chains number of parallel runs
#' @param seed initial seed. Negative values use time in seconds as seed
#' @export
generateSeeds <- function(chains=2, seed=-1) {
    if (chains < 2 || (as.integer(chains) != chains)) {
        stop("chains must be >= 2 and an integer")
    }

    if (seed < 1) {
        secs <- as.numeric(difftime(Sys.time(),
                                    paste(Sys.Date(), "00:00"), 
                                    units="secs"))
        secs <- round(secs)
        seeds <- seq_len(chains) * secs

        return(seeds)
    } else {
        seeds <- seq_len(chains) * seed

        return(seeds)
    }
}
