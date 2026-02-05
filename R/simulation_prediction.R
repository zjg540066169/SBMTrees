#' @title Simulate Continuous Longitudinal Data for Prediction
#'
#' @description Generates synthetic longitudinal data with continuous outcomes, specifically designed for
#' evaluating prediction models. The function creates a population of subjects with correlated covariates
#' and outcomes, then splits them into training and testing sets. It offers flexible options for
#' simulating non-normal random effects (e.g., skewed, mixtures, t-distributions) and residuals,
#' as well as nonlinear relationships.
#'
#' @param train_prop A numeric value between 0 and 1 indicating the proportion of the population to be used
#' for the training set. Default: \code{0.7}.
#' @param n_subject An integer specifying the total number of subjects in the population. Default: \code{1000}.
#' @param n_obs_per_sub An integer specifying the number of observations per subject. Default: \code{5}.
#' @param seed An optional integer for setting the random seed to ensure reproducibility. Default: \code{NULL}.
#' @param nonlinear A logical value. If \code{TRUE}, the outcome \code{Y} is generated using a complex
#' nonlinear function of the covariates. If \code{FALSE}, \code{Y} is a linear combination of covariates.
#' Default: \code{FALSE}.
#' @param residual A character string specifying the distribution of the residual errors added to the training outcome.
#' Options are:
#' \itemize{
#'   \item \code{"normal"}: Standard normal distribution.
#'   \item \code{"normal_mixture"}: Mixture of two normal distributions.
#'   \item \code{"skewed_normal"}: Skew-normal distribution.
#'   \item \code{"t3"}: Student's t-distribution with 3 degrees of freedom.
#'   \item \code{"t2"}: Student's t-distribution with 2 degrees of freedom.
#' }
#' @param randeff A character string specifying the distribution of the random effects.
#' Options are:
#' \itemize{
#'   \item \code{"MVN"}: Multivariate Normal distribution.
#'   \item \code{"MVN_mixture"}: Mixture of Multivariate Normal distributions.
#'   \item \code{"skewed_MVN"}: Multivariate Skew-normal distribution.
#'   \item \code{"MVT3"}: Multivariate t-distribution with 3 degrees of freedom.
#'   \item \code{"MVT2"}: Multivariate t-distribution with 2 degrees of freedom.
#' }
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{subject_id_train}{A vector of subject IDs for the training set.}
#'   \item{Z_train}{A matrix of random predictors (time/intercept) for the training set.}
#'   \item{X_train}{A matrix of covariates for the training set.}
#'   \item{Y_train}{A vector of observed outcomes for the training set (Signal + Random Effects + Residual Error).}
#'   \item{subject_id_test}{A vector of subject IDs for the testing set.}
#'   \item{Z_test}{A matrix of random predictors for the testing set.}
#'   \item{X_test}{A matrix of covariates for the testing set.}
#'   \item{Y_test}{A vector of "true" outcomes for the testing set (Signal + Random Effects), without residual error.}
#'   \item{X_pop}{A matrix of covariates for the entire population.}
#'   \item{y_pop}{A vector of "true" outcomes for the entire population (Signal + Random Effects).}
#'   \item{I}{A logical vector indicating which observations belong to the training set.}
#'   \item{X_src}{Duplicate of \code{X_train}, provided for convenience.}
#'   \item{Y_src}{Duplicate of \code{Y_train}, provided for convenience.}
#' }
#'
#' @details
#' The function first simulates correlated covariates \code{X} using a multivariate normal distribution,
#' adding subject-specific random variations. The outcome \code{Y} is then constructed based on \code{X}
#' (either linearly or nonlinearly) and combined with random effects \code{Z * Bi} drawn from the
#' specified \code{randeff} distribution.
#'
#' The data is split into training and testing sets based on \code{train_prop}. Crucially, residual noise
#' (specified by \code{residual}) is added **only** to \code{Y_train}. The \code{Y_test} values represent
#' the conditional mean (Fixed + Random Effects) and serve as the ground truth for prediction tasks
#' aiming to recover the de-noised signal.
#'
#' @examples
#' sim_data <- simulation_prediction_conti(
#'   train_prop = 0.7,
#'   n_subject = 200,
#'   n_obs_per_sub = 5,
#'   nonlinear = TRUE,
#'   residual = "normal",
#'   randeff = "skewed_MVN",
#'   seed = 123
#' )
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS mvrnorm
#' @importFrom sn rmsn rmst rsn
#' @importFrom dplyr if_else case_when
simulation_prediction_conti = function(train_prop = 0.7, n_subject = 1000, n_obs_per_sub = 5, seed = NULL, nonlinear = FALSE, residual = c("normal", "normal_mixture", "skewed_normal", "t3", "t2"), randeff = c("MVN", "MVN_mixture", "skewed_MVN", "MVT3", "MVT2")){


  if(!is.null(seed))
    set.seed(seed)
  n_obs_per_sub = sapply(1:n_subject, function(x) n_obs_per_sub)
  subject_id = c(unlist(sapply(1:n_subject, function(x) rep(x, n_obs_per_sub[x]))))
  n_obs = length(subject_id)
  Z = c(unlist(sapply(1:n_subject, function(x) 1:n_obs_per_sub[x])))
  trajectory = cbind(Z)
  Z = cbind(Z)
  Z = apply(Z, 2, scale)
  Z = cbind(1, Z)
  Z_O = cbind(Z)
  n = dim(Z)[1]


  X = mvtnorm::rmvnorm(n_obs, c(0, 2, 0, 3, 0, -1, 0))

  X[,1:7] = apply(cbind(X[,1:7]), 2, function(X_v){
    Bi = mvtnorm::rmvnorm(n_subject, rep(0, dim(Z)[2]), 0.25 * diag(ncol(Z)))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    #re
    X_v + re
  })

  if(nonlinear){
    Y = 1.5 * X[,3] - 1 * X[,4] +        # Interaction1 (scaled down)
      -0.7 * (X[,7] + 0.5)^2 +        # Interaction2 (scaled down)
      X[,1] * X[,6] +
      3 * log(3 * (X[,5])^2 + 0.5) +        # Simple addition interaction (scaled down)
      0.5 * X[,2]^2
  }else{
    Y = X[,1] - X[,2] + 3 * X[,3] - 2 * X[,4] + 0.5 * X[,5] - X[,6] + 0.5 * X[,7]
  }


  Omega <- matrix(c(1.0, 0.3,
                    0.3, 1.0),
                  nrow = 2, byrow = TRUE)
  xi <- c(0, 0)  # Location vector


  Bi = dplyr::case_when(
    randeff == "MVN" ~ MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = Omega),
    randeff == "MVN_mixture" ~ dplyr::if_else(rbinom(n_subject, 1, 0.5) == 1, MASS::mvrnorm(n = n_subject, mu = c(2,2), Sigma = Omega), MASS::mvrnorm(n = n_subject, mu = c(-2,-2), Sigma = Omega)),
    randeff == "skewed_MVN" ~ sn::rmsn(n = n_subject, xi = xi, Omega = Omega, alpha = c(2, -8)),
    randeff == "MVT3" ~ sn::rmst(n = n_subject, xi = xi, Omega = Omega, alpha = c(0, 0), nu = 3),
    randeff == "MVT2" ~ sn::rmst(n = n_subject, xi = xi, Omega = Omega, alpha = c(0, 0), nu = 2)
  )
  re = sapply(1:length(subject_id), function(x){
    Z[x,] %*% Bi[subject_id[x],]
  })
  Y = Y + re

  I = sample(1:n, size = round(n * train_prop), replace = F)
  I = replace(rep(FALSE, n), I, TRUE)

  Y_train = Y[I == TRUE]
  subject_id_train = subject_id[I == TRUE]

  res = dplyr::case_when(
    residual == "normal" ~ rnorm(n = sum(I), 0, sd = 5),
    residual == "normal_mixture" ~ dplyr::if_else(rbinom(sum(I), 1, 0.5) == 1, rnorm(n = sum(I), 5, sd = 2), rnorm(n = sum(I), -5, sd = 2)),
    residual == "skewed_normal" ~ sn::rsn(n = sum(I), xi = 0, omega = 5, alpha=8, tau = 3.5),
    residual == "t3" ~ rt(sum(I), df = 3),
    residual == "t2" ~ rt(sum(I), df = 2)
  )

  return(list(subject_id_train = subject_id_train, Z_train = Z[I,], X_train = X[I,], Y_train = Y_train + res, X_pop = X, y_pop = Y, I = I, subject_id_test = subject_id[!I],  Z_test = Z[!I,], X_test = X[!I,], Y_test = Y[!I], X_src = X[I,], Y_src = Y[I]))
}


#' @title Simulate Binary Longitudinal Data for Prediction
#'
#' @description Generates synthetic longitudinal data with binary outcomes, designed for evaluating
#' classification and prediction models. The function creates a latent continuous variable based on
#' covariates and random effects, then converts it into binary outcomes using various link functions
#' (corresponding to the \code{residual} argument).
#'
#' @param train_prop A numeric value between 0 and 1 indicating the proportion of the population to be used
#' for the training set. Default: \code{0.7}.
#' @param n_subject An integer specifying the total number of subjects in the population. Default: \code{1000}.
#' @param n_obs_per_sub An integer specifying the number of observations per subject. Default: \code{5}.
#' @param seed An optional integer for setting the random seed to ensure reproducibility. Default: \code{NULL}.
#' @param nonlinear A logical value. If \code{TRUE}, the latent variable is generated using a complex
#' nonlinear function of the covariates. If \code{FALSE}, it is a linear combination. Default: \code{FALSE}.
#' @param residual A character string specifying the link function (CDF) used to generate probabilities from the latent variable.
#' This effectively acts as the error distribution assumption in a Generalized Linear Mixed Model (GLMM) context:
#' \itemize{
#'   \item \code{"normal"}: Uses the standard normal CDF (Probit link).
#'   \item \code{"logistic"}: Uses the logistic CDF (Logit link).
#'   \item \code{"t3"}: Uses the Student's t (df=3) CDF.
#'   \item \code{"t2"}: Uses the Student's t (df=2) CDF.
#' }
#' @param randeff A character string specifying the distribution of the random effects added to the latent variable.
#' Options are:
#' \itemize{
#'   \item \code{"MVN"}: Multivariate Normal distribution.
#'   \item \code{"MVN_mixture"}: Mixture of Multivariate Normal distributions.
#'   \item \code{"skewed_MVN"}: Multivariate Skew-normal distribution.
#'   \item \code{"MVT3"}: Multivariate t-distribution with 3 degrees of freedom.
#'   \item \code{"MVT2"}: Multivariate t-distribution with 2 degrees of freedom.
#' }
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{subject_id_train}{A vector of subject IDs for the training set.}
#'   \item{Z_train}{A matrix of random predictors (time/intercept) for the training set.}
#'   \item{X_train}{A matrix of covariates for the training set.}
#'   \item{Y_train}{A vector of **observed binary outcomes** (0 or 1) for the training set.}
#'   \item{subject_id_test}{A vector of subject IDs for the testing set.}
#'   \item{Z_test}{A matrix of random predictors for the testing set.}
#'   \item{X_test}{A matrix of covariates for the testing set.}
#'   \item{Y_test}{A vector of **true probabilities** for the testing set. These represent the ground truth propensity scores (0 to 1) used for evaluation.}
#'   \item{X_pop}{A matrix of covariates for the entire population.}
#'   \item{y_pop}{A vector of true probabilities for the entire population.}
#'   \item{I}{A logical vector indicating which observations belong to the training set.}
#'   \item{X_src}{Duplicate of \code{X_train}, provided for convenience.}
#'   \item{Y_src}{Vector of true probabilities for the training set (unlike \code{Y_train} which is binary).}
#' }
#'
#' @details
#' The function simulates a latent continuous variable \eqn{Y^*} based on fixed effects (linear or nonlinear \code{X})
#' and random effects (\code{Z * Bi}). This latent variable is scaled and then transformed into a probability \eqn{p}
#' using the CDF specified by \code{residual}.
#'
#' For the training set, the observed outcome \code{Y_train} is sampled from a Bernoulli distribution
#' with probability \eqn{p}. For the testing set, the function returns the probability \eqn{p} itself (\code{Y_test}),
#' allowing for precise evaluation of the model's ability to estimate propensity scores or risk.
#'
#' @examples
#' # Simulate data with logistic link (Logit) and mixture of normal random effects
#' sim_bin <- simulation_prediction_binary(
#'   train_prop = 0.7,
#'   n_subject = 500,
#'   residual = "logistic",
#'   randeff = "MVN_mixture",
#'   seed = 123
#' )
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS mvrnorm
#' @importFrom sn rmsn rmst
#' @importFrom dplyr if_else case_when
#' @importFrom arm invlogit
#' @importFrom stats pnorm pt rbinom rnorm
simulation_prediction_binary = function(train_prop = 0.7, n_subject = 1000, n_obs_per_sub = 5, seed = NULL, nonlinear = FALSE, residual = c("normal", "logistic", "t3", "t2"), randeff = c("MVN", "MVN_mixture", "skewed_MVN", "MVT3", "MVT2")){


  if(!is.null(seed))
    set.seed(seed)
  n_obs_per_sub = sapply(1:n_subject, function(x) n_obs_per_sub)
  subject_id = c(unlist(sapply(1:n_subject, function(x) rep(x, n_obs_per_sub[x]))))
  n_obs = length(subject_id)
  Z = c(unlist(sapply(1:n_subject, function(x) 1:n_obs_per_sub[x])))
  trajectory = cbind(Z)
  Z = cbind(Z)
  Z = apply(Z, 2, scale)
  Z = cbind(1, Z)
  Z_O = cbind(Z)
  n = dim(Z)[1]


  X = mvtnorm::rmvnorm(n_obs, c(0, 2, 0, 3, 0, -1, 0))

  X[,1:7] = apply(cbind(X[,1:7]), 2, function(X_v){
    Bi = mvtnorm::rmvnorm(n_subject, rep(0, dim(Z)[2]), 0.25 * diag(ncol(Z)))
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    #re
    X_v + re
  })

  if(nonlinear){
    Y = 1.5 * X[,3] - 1 * X[,4] +        # Interaction1 (scaled down)
      -0.7 * (X[,7] + 0.5)^2 +        # Interaction2 (scaled down)
      X[,1] * X[,6] +
      3 * log(3 * (X[,5])^2 + 0.5) +        # Simple addition interaction (scaled down)
      0.5 * X[,2]^2
  }else{
    Y = X[,1] - X[,2] + 3 * X[,3] - 2 * X[,4] + 0.5 * X[,5] - X[,6] + 0.5 * X[,7]
  }


  Omega <- matrix(c(1.0, 0.3,
                    0.3, 1.0),
                  nrow = 2, byrow = TRUE)
  xi <- c(0, 0)  # Location vector


  Bi = dplyr::case_when(
    randeff == "MVN" ~ MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = Omega),
    randeff == "MVN_mixture" ~ dplyr::if_else(rbinom(n_subject, 1, 0.5) == 1, MASS::mvrnorm(n = n_subject, mu = c(2,2), Sigma = Omega), MASS::mvrnorm(n = n_subject, mu = c(-2,-2), Sigma = Omega)),
    randeff == "skewed_MVN" ~ sn::rmsn(n = n_subject, xi = xi, Omega = Omega, alpha = c(2, -8)),
    randeff == "MVT3" ~ sn::rmst(n = n_subject, xi = xi, Omega = Omega, alpha = c(0, 0), nu = 3),
    randeff == "MVT2" ~ sn::rmst(n = n_subject, xi = xi, Omega = Omega, alpha = c(0, 0), nu = 2)
  )
  re = sapply(1:length(subject_id), function(x){
    Z[x,] %*% Bi[subject_id[x],]
  })
  Y = Y + re
  Y = Y / 5

  I = sample(1:n, size = round(n * train_prop), replace = F)
  I = replace(rep(FALSE, n), I, TRUE)

  Y_train = Y[I == TRUE]
  subject_id_train = subject_id[I == TRUE]

  if(residual == "normal"){
    Y_train = rbinom(length(Y_train), 1, pnorm(Y_train, sd = 1))
    y_pop =  pnorm(Y, sd = 5)
    Y_test = pnorm(Y[!I], sd = 5)
    Y_src = pnorm(Y[I], sd = 5)
  }else if(residual == "logistic"){
    Y_train = rbinom(length(Y_train), 1, arm::invlogit(Y_train))
    y_pop =  arm::invlogit(Y)
    Y_test = arm::invlogit(Y[!I])
    Y_src = arm::invlogit(Y[I])
  }else if(residual == "t3"){
    Y_train = rbinom(length(Y_train), 1, pt(Y_train, df = 3))
    y_pop =  pt(Y, df = 3)
    Y_test = pt(Y[!I], df = 3)
    Y_src = pt(Y[I], df = 3)
  }else if(residual == "t2"){
    Y_train = rbinom(length(Y_train), 1, pt(Y_train, df = 2))
    y_pop =  pt(Y, df = 2)
    Y_test = pt(Y[!I], df = 2)
    Y_src = pt(Y[I], df = 2)
  }
  return(list(subject_id_train = subject_id_train, Z_train = Z[I,], X_train = X[I,], Y_train = Y_train, X_pop = X, y_pop = y_pop, I = I, subject_id_test = subject_id[!I],  Z_test = Z[!I,], X_test = X[!I,], Y_test = Y_test, X_src = X[I,], Y_src = Y_src))
}
