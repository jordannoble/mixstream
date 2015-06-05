#' @name mixsettings
#' @title mixsettings
#' @description Settings to configure \code{\link{mixstream}}.
#' @param x a list to be coerced or object to be verified.
#' @param tours an integer specifying the number of tours to run.
#' @param toursize an integer specifying the number of data points in each tour.
#' @param rate a double in the open interval (0.5, 1). This represents the
#' exponent to be used for the learning rate in the online EM algorithm where
#' the learning rate is of the form gamma = i^(-rate) where i is the current
#' iteration. Larger values lead to quicker convergence but less stability and
#' vice-versa for smaller values. The default is 0.6, which is often a suitable
#' compromise.
#' @param delay an integer representing the length of the inhibition phase, i.e.
#' the number of iterations of the online EM algorithm to run before the
#' parameters of the mixture model are updated. Having a small inhibition phase
#' provides a level of robustness for the early parameter estimates. The default
#' is 20.
#' @param PRavg an integer representing the tour to start Polyak-Ruppert
#' averaging on, which computes a rolling average of the the parameter estimates
#' from this value onwards. This is useful for reducing oscillatory behaviour
#' from the parameter estimates. The default is 0, which results in no
#' Polyak-Ruppert averaging being applied.
#' @param keep_history a logical. If \code{TRUE}, the parameter estimates at
#' each tour are stored and returned as the \code{history} component in the 
#' output of the \code{\link{mixstream}} function. Warning: this setting is not
#' recommended for a large number of data points due to the large memory
#' requirement. The default is \code{FALSE}.
#' @param debugtour an integer representing the tour at which to pause the
#' algorithm and open up the R browser. Useful for debugging. The default is 0,
#' i.e the algorithm is not interrupted by default.
#' @param eta a double representing a constant learning rate when using SGD to
#' estimate the regression coefficients. Defaults to 0.01.
#' @param eta1 a double representing the first parameter in a decreasing
#' learning rate of the form eta = eta1 / (1 + eta1 * eta2 * it) when using SGD
#' to estimate the regression coefficients. Default is \code{NULL}, i.e. not to
#' use the decreasing learning rate form.
#' @param eta2 a double representing the second parameter in a decreasing
#' learning rate of the form eta = eta1 / (1 + eta1 * eta2 * it) when using SGD
#' to estimate the regression coefficients. The default is \code{NULL}, i.e. not
#' to use the decreasing learning rate form.
#' @param rprop a double representing the proportion to increase or decrease the
#' learning rate when using SGD. Based on RPROP, the learning rate will be
#' increased if the sign of the previous score agrees with the sign of the
#' current score, otherwise the learning rate will be decreased. This ought to
#' provide faster convergence. Here the score refers to the derivative of the 
#' expectation of the complete log-likelihood with respect to the parameters.
#' The default is 0, i.e. by default the learning rate won't be increased or
#' decreased.
#' @param adagrad logical. If \code{TRUE}, the AdaGrad algorithm is used to 
#' continuously adapt the learning rate to the data.
#' 
#' @details Note that the \code{\link{mixstream}} function will call
#' \code{as.mixsettings} on the value of the settings argument by default and
#' so it is fine to provide a list instead of a \code{mixsettings} object.
#' Note also that \code{tours} and \code{toursize} need to be specified when
#' calling the \code{mixsettings} and \code{as.mixsettings} functions, and by
#' default \code{\link{mixstream}} sets \code{tours} to \code{n} (when \code{n} is
#' known) and \code{toursize} to 1.
#' @return Both \code{mixsettings} and \code{as.mixsettings} functions return
#' an object of class \code{mixsettings}, which is a list of all the settings
#' that can be used by the \code{\link{mixstream}} function.
#' 
#' \code{is.mixsettings} returns a logical, which is \code{TRUE} if \code{x}
#' is a valid \code{mixsettings} object.
#' @seealso \code{\link{mixstream}}
mixsettings <- function(tours, toursize, rate = 0.6, delay = 20, PRavg = 0,
                         keep_history = FALSE, debugtour = 0, eta = 0.01,
                         eta1 = NULL, eta2 = NULL, rprop = 0, adagrad = FALSE) {
  structure(mget(ls()), class='mixsettings')
}

#' @rdname mixsettings
is.mixsettings <- function(x) {
  if (!inherits(x, "mixsettings") | !is.list(x) | length(x) == 0) return(FALSE)
  req_names <- c('tours', 'toursize', 'rate', 'delay', 'PRavg', 'keep_history',
                 'debugtour', 'eta', 'eta1', 'eta2', 'rprop', 'adagrad')
  if (any(!sapply(req_names, function(a) a %in% names(x)))) return(FALSE)
  TRUE
}

#' @rdname mixsettings
as.mixsettings <- function(x, tours, toursize) {
  # default settings
  out <- mixsettings(tours, toursize)
  # replace default settings with user specified settings
  out[names(x)] <- x
  out
}