#' @name mixmodel
#' @title mixmodel
#' @description Specify the type of model for use in \code{\link{mixstream}}.
#' @param x a list to be coerced or object to be verified.
#' @param family a character string representing the exponential family of the
#' mixture densities. The current supported families are \code{"normal"},
#' \code{"poisson"}, \code{"binomial"} and \code{"mvnormal"}.
#' @param type a character string representing the type of mixture. The current
#' supported types are \code{"mixef"} for mixtures of exponential families,
#' \code{"mixglm"} for mixtures of generalised linear models.
#' @param invlink a function representing the inverse link function in the
#' generalised linear model, i.e. the function that maps the linear predictor to
#' the mean. This is NULL by default, which chooses the canonical inverse link
#' of the exponential family provided by the \code{family} argument.
#' @param dlink a function representing the derivative of the link function in
#' the generalised linear model, i.e. the derivative of the function that maps
#' the mean to the linear predictor. This is NULL by default, which chooses the
#' derivative of the canonical link of the exponential family provided by the
#' \code{family} argument.
#' @return Returns an object of class \code{mixmodel}, which is essentially a
#' list with the following components:
#' \item{family}{a character string equal to the value provided to the 
#' \code{family} argument.}
#' \item{type}{a character string equal to the value provided to the \code{type}
#' argument.}
#' \item{samp}{a function that samples from the exponential family used for this
#' mixture model.}
#' \item{sampmap}{a named vector of character strings that represent the names
#' of the arguments required by the sampling function. The names correspond to
#' the names used for these parameters when specifying parameters (e.g. in
#' \code{\link{mixstream}} or \code{\link{rmixture}}).}
#' \item{dens}{a function that computes the product of a mixture probability and
#' the density of the exponential family used for the mixture model. Note that
#' this density only needs to be correct up to a constant of proportionality.}
#' \item{suff}{a function that returns the sufficient statistic of the supplied
#' data given the exponential family of the mixture model.}
#' \item{upd}{a function that updates the parameter estimates give the
#' sufficient statistics}
#' \item{suffdims}{a numeric vector representing the }
#' @section Specifying custom mixture models
#' 
#' @seealso \code{\link{mixstream}}
mixmodel <- function(family, type, invlink = NULL, dlink = NULL) {
  
  if (missing(family) | missing(type)) {
    stop("both family and type must be specified in mixmodel().")
  }
  
  family <- match.arg(family, c("normal", "poisson", "binomial", "mvnormal"))
  type <- match.arg(type, c("mixef", "mixglm"))
  out <- structure(list(), class = "mixmodel")
  out$family <- family
  out$type <- type
  out$sgd <- FALSE
  
  # check mvtnorm is present for mvnormal family 
  if (family == "mvnormal" & !requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("the package 'mvtnorm' is required for the multivariate normal ",
         "family, please ensure it is installed.", call. = FALSE)
  }
  
  # sampling functions
  samplist <- list(
    normal = rnorm,
    poisson = rpois,
    binomial = rbinom,
    mvnormal = mvtnorm::rmvnorm)
  
  # map the parameter names to the arguments of the sampling function
  sampmap <- list(
    normal = c(props = "n", mu = "mean", sigma = "sd"),
    poisson = c(props = "n", lambda = "lambda"),
    binomial = c(props = "n", size = "size", p = "prob"),
    mvnormal = c(props = "n", mu = "mean", sigma = "sigma")
    )
  
  out$samp <- samplist[[family]]
  out$sampmap <- sampmap[[family]]
  
  if (type == "mixef") {
    
    # trimmed density functions
    denslist <- list(
      normal = function(y, params) {
        sigma <- params$sigma
        params$props * exp(-0.5 * (y - params$mu)^2 / sigma^2) / sigma
      },
      poisson = function(y, params) {
        params$props * exp(-params$lambda) * params$lambda^y
      },
      binomial = function(y, params) {
        n <- params$size
        params$props * params$p^(n * y) * (1 - params$p)^(n - n * y)
      },
      mvnormal = function(y, params) {
        params$props * mapply(function(mu, sigma) {
          tmp <- y - mu
          det(sigma)^(-0.5) * exp(-0.5 * (tmp %*% solve(sigma) %*% tmp))
        }, params$mu, params$sigma)
      })
    
    # sufficient statistic functions
    sufflist <- list(
      normal = function(y, p) {
        s <- c(1, y, y^2)
        lapply(s, `*`, p)
      },
      poisson = function(y, p) {
        s <- c(1, y)
        lapply(s, `*`, p)
      },
      binomial = function(y, p) {
        s <- c(1, y)
        lapply(s, `*`, p)
      },
      mvnormal = function(y, p) {
        list(p, lapply(p, `*`, y), lapply(p, `*`, y %*% t(y)))
      })
    
    # parameter update functions
    updlist <- list(
      normal = function(s) {
        mu <- s[[2]] / s[[1]]
        list(props = s[[1]], mu = mu, sigma = sqrt(s[[3]] / s[[1]] - mu^2))
      },
      poisson = function(s) {
        list(props = s[[1]], lambda = s[[2]] / s[[1]])
      },
      binomial = function(s) {
        list(props = s[[1]], p = s[[2]] / s[[1]])
      },
      mvnormal = function(s) {
        mu <- mapply(`/`, s[[2]], s[[1]], SIMPLIFY = FALSE)
        sigma <- mapply(function(s1, s3, mu) {
          s3 / s1 - mu %*% t(mu)
        }, s[[1]], s[[3]], mu, SIMPLIFY = FALSE)
        list(props = s[[1]], mu = mu, sigma = sigma)
      })
    
    # dimensions of sufficient statistic
    suffdimslist <- list(
      normal = c(1, 1, 1),
      poisson = c(1, 1),
      binomial = c(1, 1),
      mvnormal = c(1, expression(d), expression(c(d, d))))
    
    # names of parameters to be estimated
    paramnameslist <- list(
      normal = c("props", "mu", "sigma"),
      poisson = c("props", "lambda"),
      binomial = c("props", "p"),
      mvnormal = c("props", "mu", "sigma"))
    
    outnames <- c("dens", "suff", "upd", "suffdims", "paramnames")
    outlistnames <- paste(outnames, "list", sep = "")
    out[outnames] <- lapply(outlistnames, function(x) get(x)[[family]])
    out$suffmv <- if(family == "mvnormal") TRUE else FALSE
    
  }
  
  if (type == "mixglm") {
    
    if (family == "mvnormal") {
      stop("multivariate normal regression not yet supported.", .call = FALSE)
    }
    
    # default invlink
    if (is.null(invlink) | !is.function(invlink)) {
      invlink <- switch(family, 
                        normal = identity,
                        poisson = exp,
                        binomial = function(x) 1 / (1 + exp(-x)))
    }
    
    # default dlink
    if (is.null(dlink) | !is.function(dlink)) {
      dlink <- switch(family,
                      normal = function(x) 1,
                      poisson = function(x) 1 / x,
                      binomial = function(x) 1 / (x * (1 - x)))
    }
    
    # trimmed density functions calculated via invlink
    
    denslist <- list(
      normal = function(y, params) {
        x <- y[-1]
        eta <- sapply(params$beta, `%*%`, x)
        mu <- invlink(eta)
        sigma <- params$sigma
        params$props * exp(-0.5 * (y[1] - mu)^2 / sigma^2) / sigma
      },
      poisson = function(y, params) {
        x <- y[-1]
        eta <- sapply(params$beta, `%*%`, x)
        lambda <- invlink(eta)
        tmp <- exp(-lambda + y[1] * log(lambda))
        if (any(is.infinite(tmp))) {
          # Stirling's approximation to prevent overflow
          tmp <- exp(y[1] - lambda) * (lambda / y[1])^y[1] / sqrt(2*pi*y[1])
        }
        params$props * tmp
      },
      binomial = function(y, params) {
        x <- y[-1]
        eta <- sapply(params$beta, `%*%`, x)
        n <- params$size
        p <- invlink(eta)
        params$props * p^(n * y[1]) * (1 - p)^(n - n * y[1])
      })
    
    out$dens <- denslist[[family]]
    formals(out$dens)$invlink <- invlink
    out$invlink = invlink
    
    # linear regression can exploit closed form of sufficient statistics
    if (family == "normal" & isTRUE(all.equal(invlink, identity))) {
      
      out$suff <- function(y, p) {
        x <- y[-1]
        xtx <- x %*% t(x)
        list(p, lapply(p, `*`, x * y[1]), lapply(p, `*`, xtx), p * y[1]^2)
      }
      out$upd <- function(s) {
        beta <- mapply(function(s3, s2) as.vector(solve(s3) %*% s2), s[[3]], 
                       s[[2]], SIMPLIFY = FALSE)
        sigma <- mapply(function(s4, beta, s2, s1) {
          sqrt((s4 - sum(beta * s2)) / s1)
        }, s[[4]], beta, s[[2]], s[[1]])
        list(props = s[[1]],
             beta = beta,
             sigma = sigma)
      }
      out$suffdims <- expression(1, d - 1, c(d - 1, d - 1), 1)
      out$paramnames <- c("props", "beta", "sigma")
      out$suffmv <- TRUE
      0
    } else {
      
      out$suff <- function(y, p) list(p)
      out$suffdims <- 1
      out$suffmv <- FALSE
      out$paramnames <- "props"
      out$upd <- function(s) { list(s[[1]]) }
      
      scorelist <- list(
        normal = function(y, params) {
          # not yet supported (have to deal with nuisance sigma^2)   
          NULL
        },
        poisson = function(y, params) {
          x <- y[-1]
          eta <- sapply(params$beta, `%*%`, x)
          lambda <- sapply(eta, invlink)
          tmp <- sapply(lambda, function(l) exp((y[1] - 1) * log(l) - l))
          if (any(is.infinite(tmp))) {
            # Stirling's approximation
            tmp <- exp(y[1] - lambda) * (lambda / y[1])^y[1] / lambda / 
              sqrt(2*pi * y[1])
          } 
          lapply(1:length(lambda), function(j) {
            params$props[j] * x * (y[1] - lambda[j]) * tmp[j] / dlink(lambda[j])
          })
        },
        binomial = function(y, params) {
          x <- y[-1]
          eta <- lapply(params$beta, `%*%`, x)
          p <- lapply(eta, invlink)
          n <- params$size
          lapply(1:length(p), function(j) {
            params$props[j] * n * x * p[[j]]^(n * y[1] - 1) *
              (1 - p[[j]])^(n - n * y[1] -1) * (y[1] - p[[j]]) / dlink(p[[j]])
          })
        })
      
      out$score <- scorelist[[family]]
      formals(out$score)$invlink <- invlink
      formals(out$score)$dlink <- dlink
      out$sgd <- TRUE
    }
    
  }
  
  return(out)
  
}


#' @rdname mixmodel
is.mixmodel <- function(x) {
  if (!inherits(x, "mixmodel") | !is.list(x) | length(x) == 0) return(FALSE)
  req_names <- list(
    c("dens", "paramnames", "suff", "suffdims", "upd", "suffmv"),
    c("samp"))
  if (all(!sapply(req_names, function(a) all(a %in% names(x))))) return(FALSE)
  TRUE
}

#' @rdname mixmodel
as.mixmodel <- function(x) {
  if (is.mixmodel(x)) return(x)
  if (is.character(x)) x <- as.list(x)
  if (is.list(x)) {
    class(x) <- "mixmodel"
    if (is.mixmodel(x)) return(x)
    return(do.call(mixmodel, x))
  }
  stop("Unable to coerce x to mixmodel object.", call. = FALSE)
}