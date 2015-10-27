#' Basic Metropolis algorithm
#'
#' This is a basic Metropolis algorithm
#'
#' @export
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a dataset to compute the log-density on (likelihood, posterior, etc.)
#' @param f a function computing log-densities
#' @param param a list of named parameters passed on to 'f'; parameters must be given their initial values
#' @param to.move a logical or an integer vector indicating which parts of 'param' should be moved in the MCMC
#' @param n a number of iterations for the MCMC
#'
#' @importFrom coda mcmc
metro <- function(f, param=list(x=1), to.move=TRUE, n=1e4,
                  sd=1){
    ## FIND WHICH ARGUMENTS NEED TO MOVE ##
    nParam <- length(param)
    if(nParam<1) stop("At least one parameter is needed to define the domain of f")
    if(!is.null(param)) to.move <- rep(to.move, length=nParam)
    if(is.logical(to.move)) to.move <- which(to.move)
    param.names <- names(param)
    
    ## DEFINE OVERAL DENSITY COMPUTATION ##
    logdens <- function(param){
        return(sum(do.call(f, param)))
    }

    ## pre-generate unif random variates
    RAND <- log(runif(length(to.move)*n))
    COUNTER <- 1

    ## build output
    out.log <- double(n)
    out.param <- list()
    out.param[1:n] <- param

    ## for each iteration
    for(i in 2:n){
        ## move each param
        param <- out.param[i-1]
        for(j in to.move){
            new <- param
            new[[j]] <- rnorm(1, mean=new[[j]],sd=sd)
            logratio <- logdens(new) - logdens(param)

            ## accept/reject
            if(logratio >= RAND[COUNTER]){param <- new}
            COUNTER <- COUNTER+1
        }
        ## store param
        out.param[i] <- param

        ## store log dens
        out.log[i] <- logdens(param)
    }

    ## SHAPE OUTPUT ##
    out.param <- matrix(unlist(out.param), byrow=TRUE, ncol=nParam)
    colnames(out.param) <- param.names
    out <- cbind(out.log, out.param)
    out <- mcmc(out, start=1, end=n, thin=1)
    return(out)
} # end metro
