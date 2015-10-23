#' Simulate multimodal univariate distribution
#'
#' Simulate multimodal univariate distribution using a mixture of 'k' normal distributions.
#'
#' @export
#'
#' @rdname multimod
#'
#' @aliases rmultimod dmultimod
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param n a number of values to simulate
#' @param k a number of modes for the distribution
#' @param mean a vector of means (optional)
#' @param sd a vector of standard deviations (optional)
#' @param mean.range the range of the randomly generated means
#' @param sd.range the range of the randomly generated standard deviations
#'
#' @details This is essentially a wrapper for rnorm with different default values
#'
#' @examples
#' set.seed(1)
#' x <- rmultimod(100)
#' hist(x, col="grey", nclass=50, border="white")
#'
#' @importFrom stats rnorm runif
#'
rmultimod <- function(n, k=2, mean=NULL, sd=NULL,
                      mean.range=c(0,100), sd.range=c(1,1)){
    ## GET VECTOR OF MEANS ##
    if(is.null(mean)){
        mean <- runif(k, min=min(mean.range), max=(max(mean.range)))
    }

    ## GET VECTOR OF SDS ##
    if(is.null(sd)){
        sd <- runif(k, min=min(sd.range), max=(max(sd.range)))
    }

    ## GENERATE VARIABLES ##
    out <- rnorm(n, mean=mean, sd=sd)

    ## RETURN ##
    return(out)
} # end rmultimod




#' @rdname multimod
#' @export
#'
dmultimod <- function(x, k=2, mean=NULL, sd=NULL, log=TRUE,
                      mean.range=c(0,100), sd.range=c(1,1)){
    ## GET VECTOR OF MEANS ##
    if(is.null(mean)){
        mean <- runif(k, min=min(mean.range), max=(max(mean.range)))
    }

    ## GET VECTOR OF SDS ##
    if(is.null(sd)){
        sd <- runif(k, min=min(sd.range), max=(max(sd.range)))
    }

    ## COMPUTE DENSITIES ##
    out <- double(length(x))
    for(i in 1:k){
        out <- out + dnorm(x, mean=mean[i], sd=sd[i] , log=log)
    }

    ## RETURN ##
    return(out)
} # end rmultimod

