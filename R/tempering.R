#' Variation on the tempering algorithm
#'
#' This is a variation on the tempering algorithm
#'
#' @export
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param f a function computing log-densities
#' @param param a list of named parameters passed on to 'f'; parameters must be given their initial values
#' @param to.move a logical or an integer vector indicating which parts of 'param' should be moved in the MCMC
#' @param n a number of iterations for the MCMC
#' @param sd the standard deviation of the proposal normal distribution
#' @param omega the probability of sampling from the cool distribution
#'
#' @examples
#'
#' ## try with a basic normal density
#' f1 <- function(x){dnorm(x, log=TRUE)}
#' mc1 <- metro(f1)
#' head(mc1)
#' tail(mc1)
#' plot(mc1)
#'
#' ## change initial point
#' mc2 <- metro(f1, list(x=-10))
#' plot(mc2)
#'
#' ## try with a mixture of densities
#' fmix <- function(x){log(dnorm(x, mean=-3)+dnorm(x, mean=1, sd=.5))}
#' mc3 <- metro(fmix, list(x=0))
#' plot(mc3)
#'
#' ## try harder example
#' fmix2 <- function(x){log(dnorm(x, mean=-5)+dnorm(x, mean=5, sd=.5))}
#' mc4 <- metro(fmix2, list(x=0), sd=3, n=1e4)
#' plot(mc4)
#'
#' @importFrom coda mcmc
#' @importFrom stats rnorm runif
#'
metro <- function(f, param=list(x=0), to.move=TRUE, n=1e3,
                  sd=1, omega=0.8){
    ## FIND WHICH ARGUMENTS NEED TO MOVE ##
    nParam <- length(param)
    if(nParam<1) stop("At least one parameter is needed to define the domain of f")
    if(!is.null(param)) to.move <- rep(to.move, length=nParam)
    if(is.logical(to.move)) to.move <- which(to.move)
    param.names <- names(param)

    ## CHECK OMEGA ##
    if(omega<0 || omega>1) stop("omega must be between 0 and 1")
    ## DEFINE OVERAL DENSITY COMPUTATION ##
    logdens <- function(param){
        return(sum(do.call(f, param)))
    }

    ## pre-generate unif random variates
    RAND <- log(runif(length(to.move)*n))
    IS.COOL <- sample(c(TRUE,FALSE), prob=c(omega, 1-omega)) # TRUE if cool temperature
    COUNTER <- 1

    ## build output
    out.log <- double(n)
    out.param <- list()
    out.param[1:n] <- param

    ## for each iteration
    for(i in 2:n){
        ## CASE 1: WE SAMPLE FROM THE COOL DISTRIBUTION ##
        ## (classical Metropolis algorithm)
        if(IS.COOL[i]){
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
    } else{
        ## CASE 2: WE SAMPLE FROM HEATED DISTRIBUTION ##

    }

    ## SHAPE OUTPUT ##
    out.param <- matrix(unlist(out.param), byrow=TRUE, ncol=nParam)
    colnames(out.param) <- param.names
    out <- cbind(out.log, out.param)
    out <- mcmc(out, start=1, end=n, thin=1)
    return(out)
} # end metro
