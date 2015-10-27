#' Eruptive MCMC
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
#' @param max.temp the maximum temperature to be used
#'
#' @examples
#'
#' ## try with a basic normal density
#' f1 <- function(x){dnorm(x, log=TRUE)}
#'
#' ## basic MCMC
#' plot(mc1)
#' mc1 <- metro(f1)
#'
#' ## eMCMC
#' emc1 <- emcmc(f1)
#' plot(emc1[[1]]) # all chains
#' plot(emc1[[2]]) # thinned
#'
#' ## try with a mixture of densities
#' fmix <- function(x){log(dnorm(x, mean=-5) +
#'                         dnorm(x, mean=5, sd=.5))}
#'
#' plot(function(x) exp(fmix(x)), xlim=c(-10,10))
#' title("Target density")
#'
#' ## basic MCMC
#' mc2 <- metro(fmix, list(x=0), n=6e4)
#' plot(mc2)
#'
#' ## eMCMC
#' emc2 <- emcmc(fmix)
#' plot(emc2[[1]]) # all chains
#' plot(emc2[[2]]) # thinned
#'
#' @importFrom coda mcmc
#' @importFrom stats rnorm runif
#'
emcmc <- function(f, param=list(x=0), to.move=TRUE, n=200,
                  cold.block.size=200, hot.block.size=100,
                  sd=1, min.temp=10, max.temp=100){
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

    ## compute full output size
    N <- n * (cold.block.size + hot.block.size)

    ## pre-generate unif random variates for Metropolis
    RAND <- log(runif(length(to.move)*N))
    COUNTER <- 1

    ## pre-generate hot temperatures
    alpha <- runif(n, min=min.temp, max=max.temp)

    ## define temperature scheme
    out.alpha <- unlist(lapply(1:n, function(i) rep(c(1,alpha[i]), c(cold.block.size, hot.block.size))))

    ## build output
    out.log <- double(N)
    out.param <- list()
    out.param[1:N] <- param

    ## for each iteration
    for(i in 2:N){
        ## move parameters
        param <- out.param[i-1]
        for(j in to.move){
            new <- param
            new[[j]] <- rnorm(1, mean=new[[j]], sd=sd)
            logratio <- (logdens(new) - logdens(param))/out.alpha[i]

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
    ## with all chains
    out.param <- matrix(unlist(out.param), byrow=TRUE, ncol=nParam)
    colnames(out.param) <- param.names
    out.all <- cbind(out.log, out.alpha, out.param)
    colnames(out.all)[1:2] <- c("logdens", "alpha")

    ## thined version (keeping last of cold blocks)
    toKeep <- seq(from=cold.block.size, by=cold.block.size+hot.block.size, length=n)
    out.thin <- out.all[toKeep,]

    ## convert to mcmc
    out.all <- mcmc(out.all, start=1, end=N, thin=1)
    out.thin <- mcmc(out.thin, start=1, end=n, thin=1)

    ## return
    return(list(all=out.all, thin=out.thin))
} # end emcmc
