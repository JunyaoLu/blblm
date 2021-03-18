#' @aliases blblm-package
#' @import purrr
#' @import stats
#' @import furrr
#' @import future
#' @import utils
#' @import RcppArmadillo
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))



#' Fitting Linear Models with Bag of Little Bootstrap
#'
#' `blblm` is used to fit linear models with bag of little bootstrap.
#' It can be used to carry out regression, analysis of variance, coefficients, confident interval,
#' and to compute prediction based on the regression results.
#'
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted (for example, y ~ x).
#' @param data a data frame, list or environment containing variables in the model.
#' @param m an integer specifying the number of subsamples to be used for the bag of little bootstrap process. By default, `m = 10`.
#' @param B an integer specifying the number of bootstraps for each subsample. By default, `B = 5000`.
#' @param nthreads an integer which indicates the number of CPUs to be used for parallelization. By default, `nthreds = 1`.
#' @param fast_option logical. `TRUE` if using Rcpp code in `lm1()`. By default, `fast_option = FALSE`.
#'
#' @return `blblm` returns an object of class "blblm".
#' @export
#' @examples
#' blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' blblm(mpg ~ wt * hp, data = mtcars, m = 10, B = 1000, nthreads = 4) # using parallelization
#' @export
blblm <- function(formula, data, m = 10, B = 5000, nthreads = 1, fast_option = FALSE) {
  data_list <- split_data(data, m)
  if (nthreads > 1) {
    suppressWarnings(plan(multiprocess, workers = nthreads))
    options(future.rng.onMisuse = "ignore")
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, fast_option)
    )
  }
  else {
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, fast_option)
    )
  }
  if (!inherits(plan(), "sequential")) plan(sequential)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' Splitting Data
#'
#' Split data into m parts of approximated equal sizes.
#'
#'
#' @param data a data frame, list or environment containing variables in the model.
#' @param m an integer specifying the number of approximated equal size parts to be used for splitting data.
#'
#' @return a list of approximated equal size data.
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' Computing Estimates
#'
#' Compute the estimates for each subsample.
#'
#'
#' @param formula a symbolic description of the model to be fitted.
#' @param data a set of data.
#' @param n number of rows of the data.
#' @param B number of bootstraps.
#' @param fast_option logical. `TRUE` if using Rcpp code in `lm1()`. By default, `fast_option = FALSE`.
#'
#' @return a list of estimates for the subsample.
lm_each_subsample <- function(formula, data, n, B, fast_option) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n, fast_option), simplify = FALSE)
}


#' Computing Regression Estimates
#'
#' Compute the regression estimates for a blb dataset.
#'
#'
#' @param X independent variables of the data.
#' @param y dependent variables of the data.
#' @param n number of rows of the data.
#' @param fast_option optional to use Rcpp for faster computation. By default, `fast_option = FALSE`.
#'
#' @return a list of coefficients and sigma for the blb dataset.
lm1 <- function(X, y, n, fast_option = FALSE) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  if (fast_option) {
    fit <- fast_wlm(as.matrix(X), as.matrix(y), as.matrix(freqs))
    list(coef = fit$coefficients, sigma = mean(fit$sigma))
  }
  else {
    fit <- lm.wfit(X, y, freqs)
    list(coef = blbcoef(fit), sigma = blbsigma(fit))
  }
}


#' Computing Coefficients
#'
#' Compute the coefficients from fit.
#'
#'
#' @param fit a list containing information about the fitted model.
#'
#' @return coefficients of the model.
blbcoef <- function(fit) {
  coef(fit)
}


#' Computing Sigma
#'
#' Compute sigma from fit.
#'
#'
#' @param fit a list containing information about the fitted model.
#'
#' @return the sigma of the model.
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Printing Regression Model
#'
#' Print the formula for the blblm model.
#'
#'
#' @param x a blblm model.
#' @param ... additional arguments to be passed.
#'
#' @return formula.
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' Getting Sigma of Fit
#'
#' Get the estimates of the standard deviation for the model. Can be used to generate a confidence interval for the sigma.
#'
#'
#' @param object a blblm model.
#' @param confidence logicals. If `TRUE` the confidence interval will be calculated. By default, `confidence = FALSE`.
#' @param level the significance level for the confidence interval.
#' @param ... additional arguments to be passed.
#'
#' @return estimated sigma or a confidence interval for it.
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Getting Coefficients of Fit
#'
#' Get the coefficients for the blblm model.
#'
#'
#' @param object a blblm model.
#' @param ... additional arguments to be pased.
#'
#' @return coefficients of the model.
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' Getting the Confidence Interval
#'
#' Get the confidence interval for the coeficients.
#'
#'
#' @param object a blblm model.
#' @param parm a list of parameters used to construct confidence interval.
#' @param level the significance level for the confidence interval.
#' @param ... additional arguments to be passed.
#'
#' @return a matrix or vector with columns giving lower and upper confidence limits for each parameter.
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' Model Prediction
#'
#' Predict from the results of various blblm model fitting functions for new data.
#'
#'
#' @param object a blblm model.
#' @param new_data a set of new data used for prediction
#' @param confidence logicals. If `TRUE` the confidence interval will be calculated. By default, `confidence = FALSE`.
#' @param level the significance level for the confidence interval.
#' @param ... additional arguments to be passed.
#'
#' @return a numeric of vector for predicted value (and its confidence interval if applicable).
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}


map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}


map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}


map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
