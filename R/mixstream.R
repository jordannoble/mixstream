#' @title Online Mixture Modelling
#' @name mixstream
#' @param data a matrix or data frame of the variables in the model, where the
#' first column should be the response variable.
#' @param m an integer representing the number of finite mixtures.
#' @param model a list or \code{\link{mixmodel}} object specifying the type of
#' mixture model.
#' @param params a list of parameters for \code{mixstream} to start at that
#' should correspond to \code{model}. See 'Details' for more.
#' @param settings a list or \code{\link{mixsettings}} object specifying the settings
#' for mixstream.
#' 
#' @seealso \code{\link{mixmodel}} for specifying mixture models, 
#' \code{\link{mixsettings}} for configuring mixstream and 
#' \code{\link{rmixture}} for simulating from mixture models.
mixstream <- function(data, m, model, params, settings) {
  UseMethod('mixstream')
}

mixstream.default <- function(data, m, model, params, settings) {
  
  data <- as.data.frame(data)
  n <- nrow(data)
  d <- ncol(data)
  
  # default mixstream settings
  settings <- as.mixsettings(settings, tours = n, toursize = 1)
  
  # inherit appropriate mixture model properties
  model <- as.mixmodel(model)    
  
  # save model, initial values and settings for output
  init <- list(model = model, params = params, m = m, settings = settings)
  
  # WITH (keep settings and model elements as local vars)
  with(data = c(settings, model), {
  
  # initialise and allocate (constant) storage
  p <- rep(0, m)
  s <- lapply(suffdims, function(x) replicate(m, array(0, eval(x)),
                                           simplify = identical(x, 1)))
  prevs <- s
  names(s) <- paramnames
  r <- length(suffdims)
  
  history <- list(params)
  prevparams <- params
  loglik <- 0
  its <- tours * toursize
  prevsign <- 0
  zeta <- 0
  adaG <- 0
  terminate <- FALSE
  
  ### new touring system
  
  for(tour in 1:tours) {
    
    if (toursize == 1) {
      batch <- as.matrix(data[tour, ])
    } else {
      firstit <- (tour - 1) * toursize + 1
      lastit <- tour * toursize
      batch <- as.matrix(data[firstit:lastit, ])
    }
    
    batchsbar <- 0
    batchscore <- 0
    
    for(i in 1:toursize) {
      
      wg <- dens(batch[i, ], params)
      if(all(wg == 0)) { 
        warning("Iter: ", i, " in tour ", tour,
                " has 0 probability for all components.")
        terminate <- TRUE
        break
      }
      
      lik <- sum(wg)
      loglik <- loglik + log(lik)
      p <- wg / lik
      sbar <- suff(batch[i, ], p)
      batchsbar <- mixstream_agg(sbar, batchsbar)
      
      if (sgd) {
        # zero probability mixture components don't contribute to beta updates
        A <- which(wg != 0)
        wA <- wg[A]
        tmp <- rep(0, m)
        tmp[A] <- ((1 + log(wA)) * lik - sum(wA * log(wA))) / lik / lik
        grad <- mapply(`*`, score(unname(batch[i, ]), params), tmp,
                             SIMPLIFY = FALSE)
        batchscore <- mixstream_agg(grad, batchscore)
        ## add termination check for infinite batch score?
        
        if (any(is.na(unlist(batchscore)))) {
          warning("beta diverged to infinity.")
          break
          terminate <- TRUE
        }
      }
      
    }
    
    # check for any convergence issues
    if (terminate) break
    
    gamma <- if (delay < tour * toursize) (tour)^(-rate) else 1 / tour
    
    # Stochastic E-step
    s <- mixstream_agg(s, batchsbar, 1 - gamma, gamma / toursize)
    
    # M-step (update params)
    if(delay <= tour * toursize) params[paramnames] <- upd(s)
    
    # SGD update
    if (sgd) {
    
      # decreasing learning rate eta
      if(all(!is.null(c(eta1, eta2)))) eta <- eta1 / (1 + eta1 * eta2 * tour)
      
      #AdaGrad
      if (adagrad) {
        adaG <- mapply(function(u, v) {
          sqrt(u^2 + (v / toursize)^2)
        }, adaG, batchscore, SIMPLIFY = FALSE)
        batchscore <- mapply(`/`, batchscore, adaG, SIMPLIFY = FALSE)
      }      
      
      # RPROP
      if (rprop > 0) {
        newsign <- lapply(batchscore, sign)
        zeta <- mapply(function(i, j, k) 1 + rprop * sign(i*j) * k,
                       newsign, prevsign, zeta, SIMPLIFY = FALSE)
        batchscore <- mapply(`*`, batchscore, zeta, SIMPLIFY = FALSE)
        prevsign <- newsign 
      }

      params$beta <- mixstream_agg(params$beta, batchscore, 1, eta / toursize)
      
    }
    
    # Polyak-Ruppert averaging
    if(PRavg > 0) {
      if (tour >= PRavg) {
        params <- rolling_avg(params, prevparams, tour, PRavg)
      }
      prevparams <- params
    }
    
    # store parameters in history
    if (keep_history) history[[tour + 1]] <- params
    
    # debug
    if(tour == debugtour) browser()
    
  }
  
  output <- list(history = NULL, init = init, loglik = loglik, params = params)
  class(output) <- "mixstream"
  
  # reformat history
  if (keep_history) {
    history <- lapply(names(params), function(pn) lapply(history, `[[`, pn))
    names(history) <- names(params)
    output$history <- history
  }
  
  output
  
  }) # end with
  
}

mixstream_agg <- function(x, y, a = 1, b = 1) {
  if (identical(y, 0)) return(x)
  mapply(function(xi, yi) {
    if (!is.list(xi)) {
      a * xi + b * yi
    } else {
      mapply(function(xij, yij) {
        a * xij + b * yij 
      }, xi, yi, SIMPLIFY = FALSE) 
    }
  }, x, y, SIMPLIFY = FALSE) 
}

rolling_avg <- function(newparams, prevparams, it, n0) {
  
  # multivariate params?
  if (any(sapply(c(newparams, prevparams), is.list))) {
    out <- mapply(function(newpn, prevpn) { 
      mapply(function(newparam, prevparam) {
        (newparam + (it - n0) * prevparam) / (it - n0 + 1)
      }, newpn, prevpn, SIMPLIFY = !is.list(newpn))
    }, newparams, prevparams, SIMPLIFY = FALSE)
  } else {
    out <-mapply(function(newpn, prevpn) { 
      (newpn + (it - n0) * prevpn) / (it - n0 + 1)
    }, newparams, prevparams, SIMPLIFY = FALSE)
  }
  out
}

plot.mixstream <- function(x, paramname = "props", trueparams = NULL, ...) {
  
  if (is.null(x$history)) stop("mixstream object has no history.")
  
  ph <- x$history[[paramname]]
  if (is.null(ph)) stop("no parameter named ", paramname, ".")
  if (!is.null(dim(ph))) stop("matrix parameters not supported.")
  
  m <- x$init$m
  pl <- length(ph[[1]][[1]])
  
  trueline <- !is.null(trueparams)
  if (!is.null(trueparams[[paramname]])) trueparams <- trueparams[[paramname]]
  
  # prettier plotting
  par(mar = c(5.1, 6.15, 4.1, 2.1), cex.lab = 1.5)
  if(paramname == 'props') paramname <- 'omega'

  if(pl == 1) {
    par(mfrow=c(m, 1))
    for(i in 1:m) {
      plot(sapply(1:length(ph), function(j) ph[[j]][i]), type='l',
           ylab = bquote(hat(.(as.symbol(paramname)))[.(i)]),
           xlab = "Tour",
           ...)
      if (trueline) abline(h = trueparams[i], lty = 2, col = 'red')
    }
  } else {
    par(mfrow=c(pl*m, 1))
    for(i in 1:m) {
      for(k in 1:pl) {
        plot(sapply(1:length(ph), function(j) ph[[j]][[i]][k]), type='l',
             ylab = bquote(hat(.(as.symbol(paramname)))[list(.(i),.(k))]),
             xlab = "Tour",
             ...) 
        if (trueline) {
          abline(h = trueparams[[i]][k], lty = 2, col = 'red')
        }
      }
    }
  }

}

mixstream_tune <- function(data, m, model, params,
                           settings = NULL, n = NULL, etas = NULL) {
  
  data <- as.data.frame(data)
  if (is.null(n)) n <- nrow(data) else data <- data[1:n, ]
  
  if (is.null(etas)) etas <- 1*10^-(6:0)
  logliks <- rep(-Inf, length(etas))
  newsettings <- settings
  for(i in 1:length(etas)) {
    newsettings$eta <- etas[i]
    ms <- mixstream(data, m, model, params, newsettings)
    logliks[i] <- ms$loglik
  }
  
  bestit <- which.max(logliks)
  list(loglik = logliks[bestit], eta = etas[bestit])
  
}

test_fit <- function(x, newdata) {
  
  newdata <- as.data.frame(newdata)
  model <- x$init$model
  sum(log(colSums(apply(newdata, 1, function(r) model$dens(r, x$params)))))
  
}

mixstream_init <- function(data, m, model, params = NULL,
                           settings = NULL, n = NULL, its = 60,
                           noise = 0) {
  
  data <- as.data.frame(data)
  if (is.null(n)) n <- nrow(data) else data <- data[1:n, ]
  
  if(is.null(params$props)) {
    propslst <- lapply(1:its, function(x) {
      tmp <- runif(m)
      tmp / sum(tmp)
    })
  } else {
    propslst <- list(params$props)
  }
  
  if(model$type == "mixglm") {
    
    if(is.null(params$beta)) {
      
      data_s <- data[order(data[, 1]), ]    
      betalst <- lapply(propslst, function(props) {
        bins <- round(n * cumsum(c(0, props)))
        coeffs <- lapply(1:m, function(j) {
          binsq <- (bins[[j]] + 1):bins[[j + 1]]
          tmpdf <- data.frame(y = data_s[binsq, 1], data_s[binsq, 2:ncol(data_s)])
          unname(glm(y ~ . + 0, family = model$family,
                     data = tmpdf)$coefficients) * (1 + rnorm(1) * noise)
        })
        
      })
      
    } else {
      betalst <- params$beta
    }
  } else if(model$type == "mixef") {
    
    # NEED TO DO K-MEANS HERE TO ESTIMTE MIXEF PARAMS
#     if(is.null(params$)) {
#       
#       data_s <- data[order(data[, 1]), ]    
#       betalst <- lapply(propslst, function(props) {
#         bins <- round(n * cumsum(c(0, props)))
#         coeffs <- lapply(1:m, function(j) {
#           binsq <- (bins[[j]] + 1):bins[[j + 1]]
#           tmpdf <- data.frame(y = data_s[binsq, 1], data_s[binsq, 2:ncol(data_s)])
#           unname(glm(y ~ . + 0, family = model$family,
#                      data = tmpdf)$coefficients) * (1 + rnorm(1) * noise)
#         })
#         
#       })
#       
#     }
#     
  }  
  
  logliks <- rep(-Inf, its)
  newparams <- params
  for(i in 1:length(propslst)) {
    newparams$props <- propslst[[i]]
    newparams$beta <- betalst[[i]]
    if (!all(is.na(unlist(newparams$beta)))) {
      ms <- mixstream(data, m, model, newparams, settings)
      logliks[i] <- ms$loglik 
    }
  }
  
  bestit <- which.max(logliks)
  out <- params
  out$props <- propslst[[bestit]]
  out$beta <- betalst[[bestit]]
  
  out
  
}