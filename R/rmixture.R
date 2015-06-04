#' @title Simulating from a mixture model.
#' @description \code{rmixture} is used to simulate from a mixture model by
#' specifying the type of model (e.g. a mixture of exponential families or
#' generalised linear models) and the necessary parameters of the model type.
#' @param n an integer representing the number of data points to simulate.
#' @param m an integer representing the number of mixture components.
#' @param model a \code{\link{mixmodel}} object representing the type of mixture
#' and exponential family to be used.
#' @param params a list of parameters that should correspond to \code{model} and
#' may include the mixing proportions. See 'Details' for more.
#' @param data a matrix or data frame representing the design matrix of the
#' covariates. Only required when simulating mixtures of generalised linear
#' models.
#' @param shuffle logical. If \code{TRUE} then shuffle the output.
#' 
#' @details The value provided to the \code{params} argument should be a list
#' with components representing the parameters of the mixture model given by
#' \code{model} (see \code{\link{mixmodel}} for more information). The
#' proportions of each mixture can be specified in the \code{"props"} component
#' of \code{params}. This should be a numeric vector, where each entry is
#' between 0 and 1 such that the sum of all entries is 1. If the \code{"props"} 
#' component is missing then it assumed that all mixing propotions are equal.
#' Similarly, if any proportions are not specified, the remaining proportion
#' will be evenly distributed amongst the unspecified proportions.
#' 
#' @return Returns a list containing the following components:
#'  \item{data}{a data frame, which consists of the simulated mixture model data
#'  . When the model is of type \code{"mixef"} this will be the simulated random
#'  variables. For models of type \code{"mixglm"} the data frame will contain
#'  the response variables in the first column and the covariates in the
#'  remaining columns.}
#'  \item{component}{a vector of integers from 1 to \code{m} representing the
#'  mixture component of the associated row in \code{data}.}
#' @seealso \code{\link{mixmodel}}
rmixture <- function(n, m, model, params, data = NULL, shuffle = TRUE) {
  
  rfn <- model$samp
  props <- params$props
  if (is.null(props)) props <- c()
  
  messages <- c()
  nparams <- unname(length(formals(rfn)) - 1)
  
  # ensure model type is valid
  type = model$type
  stopifnot(type %in% c('mixef', 'mixglm'))
  
  # ensure mixglm critera satisfied
  if (type == "mixglm") {
    if (is.null(data)) stop("no covariate data provided.")
    if (is.null(params$beta)) stop("no regression coefficients provided.")
  }
  
  # ensure proportions are probabilities
  if (any(props > 1) | any(props < 0)) {
    stop("Proportions should be between 0 and 1.")
  }
  
  # distribute remaining proportion evenly between unprovided props
  props_length <- length(props)
  if (props_length != m) {
    props_missing <- m - props_length
    props_estimate <- (1 - sum(props)) / props_missing
    props <- c(props, rep(props_estimate, props_missing))
    messages <- c(messages,
                  paste("Replaced",
                        props_missing,
                        "missing proportion(s) with the value",
                        props_estimate))
  }
  
  # ensure proportions sum to 1
  if (!isTRUE(all.equal(sum(props), 1))) {
    stop("Proportions should sum to 1.")
  }
  
  # number of points to sample per mixture (rounded)
  mixture_split <- round(n * props)

  # adjust to ensure n points are sampled in total
  # (add/remove randomly to mixture_split according to proportions)
  points_diff <- n - sum(mixture_split)
  to_adjust <- 1 + findInterval(runif(abs(points_diff)), cumsum(props))
  mixture_split <- mixture_split + sign(points_diff) * tabulate(to_adjust, m)
  component <- rep(1:m, mixture_split)
  
  # sample from each mixture accordingly
  pos <- c(0, cumsum(mixture_split))
  resp <- data.frame()
  for (j in 1:m) {
    mix_params <- lapply(params[which(names(params) != "props")], `[[`, j)
    mix_params <- c(props = mixture_split[j], mix_params)
    if (type == "mixglm") {
      mix_params$beta <- NULL
      # mean = invlink(linear predictor)
      mu <- model$invlink(rowSums(sweep(data[(pos[j] + 1):pos[j + 1], ],
                                        2, params$beta[[j]], `*`)))
      missing_name <- setdiff(names(model$sampmap), names(params))
      mix_params[[missing_name]] <- mu
    }
    names(mix_params) <- model$sampmap[names(mix_params)]
    resp <- rbind(resp, data.frame(y = do.call(rfn, mix_params)))
  }
  
  # data frame formatting
  if (ncol(resp) > 1) colnames(resp) <- paste(rep("y"), 1:ncol(resp), sep = "")
  if (type == "mixglm") resp <- cbind(resp, data)
  data <- resp
  
  if (shuffle) {
    tmp <- sample(1:n)
    data <- data[tmp, ]
    component <- component[tmp]
  }
  
  # print messages
  for(msg in messages) {
    message(msg, ".")
  }
  
  list(data = data, component = component)
  
}