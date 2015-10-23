#' Basic Metropolis algorithm
#'
#' This is a basic Metropolis algorithm
#'
#' @export
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a dataset to compute the log-density on (likelihood, posterior, etc.)
#' @param f a function computing log-densities
#' @param param a list of named parameters passed on to 'f'
#' @param n a number of iterations for the MCMC
#'
#' @importFrom coda mcmc
metro <- function(f, param=list(x=0), n=1e4,
                  sd=1){
    ## FIND WHICH ARGUMENTS NEED TO MOVE ##
    toMove <- which(sapply(param, is.numeric))
    nParam <- length(toMove)

    ## DEFINE OVERAL DENSITY COMPUTATION ##
    logdens <- function(param){
        return(sum(do.call(f, param)))
    }

    ## pre-generate unif random variates
    RAND <- log(runif(length(toMove)*n))
    COUNTER <- 1

    ## build output
    out.log <- double(n)
    out.param <- matrix(double(n*nParam),ncol=nParam)
    colnames(out.param) <- names(toMove)

    ## init param
    out.param[1,] <- rnorm(nParam, mean=0, sd=20)

    ## for each iteration
    for(i in 2:n){
        ## move each param
        param <- out.param[i-1,]
        for(j in 1:nParam){
            new <- param
            new[j] <- rnorm(1, mean=new[j],sd=sd)
            logratio <- logdens(new) - logdens(param)

            ## accept/reject
            if(logratio >= RAND[COUNTER]){param <- new}
            COUNTER <- COUNTER+1
        }
        ## store param
        out.param[i,] <- param
    }

    ## SHAPE OUTPUT ##
    out <- cbind(out.log, out.param)
    names(out)[1] <- "logdens"
    out <- mcmc(out, start=1, end=n, thin=1)
    return(out)
} # end metro
