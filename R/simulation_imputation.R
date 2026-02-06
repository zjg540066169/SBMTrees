#' @title Simulate Longitudinal Data with Missing Values for Imputation
#'
#' @description Generates synthetic longitudinal data specifically designed to evaluate
#' missing data imputation methods. The function creates a complex dataset with:
#' \itemize{
#'   \item **Time-varying covariates** with autoregressive structures and random effects.
#'   \item **Non-linear relationships** and interactions between covariates.
#'   \item **Mixed data types** (continuous and binary/logical).
#'   \item **Non-normal Distributions** (optional) for both random effects and residuals (Skew-t, t-distribution).
#'   \item **Missing Data Mechanisms**:
#'   \itemize{
#'      \item \emph{Intermittent Missingness}: Generated via logistic models conditioned on outcomes and other covariates.
#'      \item \emph{Loss to Follow-up (LTFU)}: Simulates subject dropout starting from time point 4 based on values at time point 3.
#'   }
#' }
#'
#' @param NNY A logical value. If \code{TRUE}, the outcome \code{Y} is generated using non-normal distributions
#' (Skew-t random effects, t-distribution residuals). If \code{FALSE}, it uses standard Normal distributions.
#' Default: \code{TRUE}.
#' @param NNX A logical value. If \code{TRUE}, the covariates \code{X_7} through \code{X_12} are generated using
#' non-normal distributions (Mixture models, Skew-t random effects). If \code{FALSE}, they use standard Normal distributions.
#' Default: \code{TRUE}.
#' @param n_subject An integer specifying the number of subjects. Default: \code{1000}.
#' @param seed An optional integer for setting the random seed to ensure reproducibility. Default: \code{NULL}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{data_E}{A data frame of the **complete** data (ground truth) without any missing values.}
#'   \item{data_M}{A data frame of the **incomplete** data, containing \code{NA}s introduced by intermittent missingness and dropout.}
#'   \item{data_O}{A duplicate of \code{data_E} used internally for generating missingness probabilities.}
#'   \item{Z}{A matrix of random predictors (intercept and time slopes) used in generation.}
#'   \item{pair}{A matrix summarizing the missing data pattern (generated via \code{mice::md.pattern}).}
#' }
#'
#' @details
#' The simulation process creates 12 covariates (\code{X_1} to \code{X_12}):
#' \itemize{
#'   \item \code{X_1} to \code{X_6}: Base covariates generated via multivariate normal distributions with autoregressive sigma. \code{X_4, X_5, X_6} are converted to binary.
#'   \item \code{X_7} to \code{X_12}: Derived covariates dependent on the base set, involving non-linear transformations (squares, logs, interactions).
#' }
#'
#' Missingness is introduced in two stages:
#' 1. **Intermittent Missingness**: For variables \code{X_7} to \code{X_12}, missingness indicators are drawn from Bernoulli distributions where the probability depends on the outcome \code{Y} and other covariates.
#' 2. **Dropout**: A "Loss to Follow-up" indicator is generated based on data at time point 3. If a subject drops out, all values for time points 4 and 5 become \code{NA}.
#'
#' @examples
#' # Simulate data with non-normal errors and random effects
#' sim_data <- simulation_imputation(NNY = TRUE, NNX = TRUE, n_subject = 10, seed = 123)
#'
#' # View missing data pattern
#' sim_data$pair
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS mvrnorm
#' @importFrom sn rmsn rmst
#' @importFrom dplyr if_else tibble
#' @importFrom arm invlogit
#' @importFrom mice md.pattern
#' @importFrom stats rnorm rbinom rt rlogis
simulation_imputation = function(NNY = TRUE, NNX = TRUE, n_subject = 1000, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  n_obs_per_sub = 5
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


  Omega <- matrix(c(1.0, 0.3,
                    0.3, 1.0),
                  nrow = 2, byrow = TRUE)

  X = mvtnorm::rmvnorm(n_obs, rep(0, 6), sigma = AR(0.5, 6) * 0.5)

  X[,1:6] = apply(cbind(X[,1:6]), 2, function(X_v){
    Bi <- MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  X[,4:6] = (X[,4:6] >= 0)

  X_7 = cbind(-0.5 * (X[,1] + X[,4])^2 + 1 * X[,4])
  X_7 = apply(X_7, 2, function(X_v){
    if(NNX){
      Bi = dplyr::if_else(rbinom(n_subject, 1, 0.5) == 1, MASS::mvrnorm(n = n_subject, mu = c(1,1), Sigma = 0.25 * Omega), MASS::mvrnorm(n = n_subject, mu = c(-1,-1), Sigma = 0.25 * Omega))
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })


  X_8 = cbind(0.8 * X[,2] * X[,3] - X[,3])
  X_8 = apply(X_8, 2, function(X_v){
    if(NNX){
      Bi = dplyr::if_else(rbinom(n_subject, 1, 0.5) == 1, MASS::mvrnorm(n = n_subject, mu = c(2,2), Sigma = 0.25 * Omega), MASS::mvrnorm(n = n_subject, mu = c(-2,-2), Sigma = 0.25 * Omega))
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  X_9 = cbind(log(1e-80 + (X[,3] - X[,6] - 0.5)^2))
  X_9 = apply(X_9, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 3)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })

  if(NNX){
    X_7 = X_7 + dplyr::if_else(rbinom(n_obs, 1, 0.5) == 1, rnorm(n = n_obs, 2, sd = 0.5), rnorm(n = n_obs, -2, sd = 0.5))
    X_8 = X_8 + rt(n_obs, df = 3)
    X_9 = X_9 + rt(n_obs, df = 2)
  }else{
    X_7 = X_7 + rnorm(n_obs, 0, 1)
    X_8 = X_8 + rnorm(n_obs, 0, 1)
    X_9 = X_9 + rnorm(n_obs, 0, 1)
  }


  X_10 = cbind(0.5 * X[,2] - X[,5] - 0.8 * X_7)
  X_10 = apply(X_10, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 3)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })

  X_11 = cbind(X[,3] * X[,4] - 0.5 * X[,6] * X_9)
  X_11 = apply(X_11, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 2)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })


  X_12 = cbind(0.8 *X[,6] * X_8 - X[,1])
  X_12 = apply(X_12, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 2)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })



  if(NNX){
    X_10 = (X_10 + rlogis(n_obs)) >= 0
    X_11 = (X_11 + rt(n_obs, df = 3)) >= 0
    X_12 = (X_12 + rt(n_obs, df = 2)) >= 0
  }else{
    X_10 = (X_10 + rnorm(n_obs, 0, 1)) >= 0
    X_11 = (X_11 + rnorm(n_obs, 0, 1)) >= 0
    X_12 = (X_12 + rnorm(n_obs, 0, 1)) >= 0
  }


  Y = cbind(10 + 1* X[,1] +1* X[,2] + 1*X[,3] + 1 *X[,4] + 1 *X[,5] + 1 *X[,6] + 1 *X_7 + 1 *X_8 + 1 *X_9 + 1*X_10 + 1*X_11 + 1*X_12)
  Y = apply(Y, 2, function(X_v){
    if(NNY){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = 0.5 * Omega, alpha = c(0, 0), nu = 3)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 1 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })



  if(NNY){
    Y = Y + rt(n_obs, df = 2)
  }else{
    Y = Y + rnorm(n_obs, 0, 2)
  }




  data_E = dplyr::tibble(
    subject_id = subject_id,
    trajectory = trajectory[,1],
    X_1 = X[,1],
    X_2 = X[,2],
    X_3 = X[,3],
    X_4 = as.logical(X[,4]),
    X_5 = as.logical(X[,5]),
    X_6 = as.logical(X[,6]),
    X_7 = X_7[,1],
    X_8 = X_8[,1],
    X_9 = X_9[,1],
    X_10 = X_10[,1],
    X_11 = X_11[,1],
    X_12 = X_12[,1],
    Y = Y[,1]
  )

  data_O = dplyr::tibble(
    subject_id = subject_id,
    trajectory = trajectory[,1],
    X_1 = X[,1],
    X_2 = X[,2],
    X_3 = X[,3],
    X_4 = as.logical(X[,4]),
    X_5 = as.logical(X[,5]),
    X_6 = as.logical(X[,6]),
    X_7 = X_7[,1],
    X_8 = X_8[,1],
    X_9 = X_9[,1],
    X_10 = X_10[,1],
    X_11 = X_11[,1],
    X_12 = X_12[,1],
    Y = Y[,1]
  )

  R_7 = logical(n_obs)
  R_8 = logical(n_obs)
  R_9 = logical(n_obs)
  R_10 = logical(n_obs)
  R_11 = logical(n_obs)
  R_12 = logical(n_obs)
  #R_Y = logical(n_obs)

  alpha = -5 #-1
  for(j in 1:n_obs_per_sub[1]){
    if(j == 1){
      R_7[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == 1]) - 5.5 * data_O$X_6[data_O$trajectory == 1] + 1.5))
      R_8[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == 1]) - 3 * data_O$X_5[data_O$trajectory == 1] + 1.5))
      R_9[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == 1]) + 3 * data_O$X_4[data_O$trajectory == 1] - 1.5))
      R_10[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == 1]) + 2 * data_O$X_3[data_O$trajectory == 1] - 0.5))
      R_11[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == 1]) - 5 * (data_O$X_2[data_O$trajectory == 1] + 1)^2 + 3))
      R_12[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == 1]) - 3 * (data_O$X_1[data_O$trajectory == 1] + 1)^2 + 2.5))
      #R_Y[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(alpha + -1.5 * (data_O$Y[data_O$trajectory == 1]) + 1*data_O$X_1[data_O$trajectory == 1] + 1*data_O$X_2[data_O$trajectory == 1] + 1*data_O$X_3[data_O$trajectory == 1] + 1*data_O$X_4[data_O$trajectory == 1] + 1*data_O$X_5[data_O$trajectory == 1] + 1*data_O$X_6[data_O$trajectory == 1] ))
    }else{
      R_7[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == j]) - 5.5 * data_O$X_6[data_O$trajectory == j] + 1 * data_O$X_7[data_O$trajectory == (j-1)] + 1.5))
      R_8[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == j]) - 3 * data_O$X_5[data_O$trajectory == j] + 1 * data_O$X_8[data_O$trajectory == (j-1)] + 1.5))
      R_9[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == j]) +3 * data_O$X_4[data_O$trajectory == j] + 1 * data_O$X_9[data_O$trajectory == (j-1)] - 1.5))
      R_10[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == j]) + 2 * data_O$X_3[data_O$trajectory == j] + 1 * data_O$X_10[data_O$trajectory == (j-1)]  - 0.5))
      R_11[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == j]) - 5 * (data_O$X_2[data_O$trajectory == j] + 1)^2 + 1 * data_O$X_11[data_O$trajectory == (j-1)] + 3))
      R_12[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-0.5 * (data_O$Y[data_O$trajectory == j]) - 3 * (data_O$X_1[data_O$trajectory == j] + 1)^2 + 1 * data_O$X_12[data_O$trajectory == (j-1)] + 2.5))


      #R_Y[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(alpha + -1.5 * (data_O$Y[data_O$trajectory == j]) + 1*data_O$X_1[data_O$trajectory == j] + 1*data_O$X_2[data_O$trajectory == j] + 1*data_O$X_3[data_O$trajectory == j] + 1*data_O$X_4[data_O$trajectory == j] + 1*data_O$X_5[data_O$trajectory == j] + 1*data_O$X_6[data_O$trajectory == j] - 1 * data_O$Y[data_O$trajectory == (j-1)]))

      R_7[data_O$trajectory == j] = R_7[data_O$trajectory == j] * (1 - R_7[data_O$trajectory == (j-1)])
      R_8[data_O$trajectory == j] = R_8[data_O$trajectory == j] * (1 - R_8[data_O$trajectory == (j-1)])
      R_9[data_O$trajectory == j] = R_9[data_O$trajectory == j] * (1 - R_9[data_O$trajectory == (j-1)])
      R_10[data_O$trajectory == j] = R_10[data_O$trajectory == j] * (1 - R_10[data_O$trajectory == (j-1)])
      R_11[data_O$trajectory == j] = R_11[data_O$trajectory == j] * (1 - R_11[data_O$trajectory == (j-1)])
      R_12[data_O$trajectory == j] = R_12[data_O$trajectory == j] * (1 - R_12[data_O$trajectory == (j-1)])

      #R_Y[data_O$trajectory == j] = R_Y[data_O$trajectory == j] * (1 - R_Y[data_O$trajectory == (j-1)])

    }
  }

  data_M = dplyr::tibble(
    subject_id = subject_id,
    trajectory = trajectory[,1],
    X_1 = X[,1],
    X_2 = X[,2],
    X_3 = X[,3],
    X_4 = as.logical(X[,4]),
    X_5 = as.logical(X[,5]),
    X_6 = as.logical(X[,6]),
    X_7 = X_7[,1],
    X_8 = X_8[,1],
    X_9 = X_9[,1],
    X_10 = X_10[,1],
    X_11 = X_11[,1],
    X_12 = X_12[,1],
    Y = Y[,1]
  )
  data_M$X_7[R_7 == 1] = NA
  data_M$X_8[R_8 == 1] = NA
  data_M$X_9[R_9 == 1] = NA
  data_M$X_10[R_10 == 1] = NA
  data_M$X_11[R_11 == 1] = NA
  data_M$X_12[R_12 == 1] = NA


  R_ltfu = rbinom(n = n_subject, size = 1, prob = arm::invlogit(-9.5 + 0.5 * (data_O$Y[data_O$trajectory == 3]) + 2*data_O$X_1[data_O$trajectory == 3] + 2*data_O$X_2[data_O$trajectory == 3] + 2*data_O$X_3[data_O$trajectory == 3]))
  for (i in 1:n_subject) {
    if(R_ltfu[i] == 1){
      data_M[data_M$subject_id == i & data_M$trajectory %in% c(4,5), 3:(dim(data_M)[2])] = NA
    }
  }
  return(list(data_E = data_E, data_M = data_M, data_O = data_O, Z = Z, pair = mice::md.pattern(data_M[,3:dim(data_M)[2]], plot = FALSE)))
}











#' @title Simulate Longitudinal Data with Loss to Follow-up (LTFU) for Imputation
#'
#' @description Generates synthetic longitudinal data specifically designed to stress-test
#' imputation methods against **Loss to Follow-up (Dropout)**. While it includes intermittent
#' missingness, the parameters are tuned to simulate scenarios where subjects permanently leave
#' the study based on their characteristics at specific time points.
#'
#' @param NNY A logical value. If \code{TRUE}, the outcome \code{Y} is generated using non-normal distributions
#' (Skew-t random effects, t-distribution residuals). If \code{FALSE}, it uses standard Normal distributions.
#' Default: \code{TRUE}.
#' @param NNX A logical value. If \code{TRUE}, the covariates \code{X_7} through \code{X_12} are generated using
#' non-normal distributions (Mixture models, Skew-t random effects). If \code{FALSE}, they use standard Normal distributions.
#' Default: \code{TRUE}.
#' @param n_subject An integer specifying the number of subjects. Default: \code{1000}.
#' @param seed An optional integer for setting the random seed to ensure reproducibility. Default: \code{NULL}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{data_E}{A data frame of the **complete** data (ground truth) without any missing values.}
#'   \item{data_M}{A data frame of the **incomplete** data, containing \code{NA}s introduced by intermittent missingness and significant LTFU.}
#'   \item{data_O}{A duplicate of \code{data_E} used internally for generating missingness probabilities.}
#'   \item{Z}{A matrix of random predictors (intercept and time slopes) used in generation.}
#'   \item{pair}{A matrix summarizing the missing data pattern (generated via \code{mice::md.pattern}).}
#' }
#'
#' @details
#' The data generation process mirrors \code{\link{simulation_imputation}} regarding covariate structure (time-varying, non-linear, mixed types),
#' but utilizes specific coefficients to drive the missingness mechanisms:
#'
#' **1. Loss to Follow-up (LTFU):**
#' Dropout is simulated based on the subject's state at **time point 3**. A logistic model determines the probability of dropout using:
#' \itemize{
#'   \item The outcome \code{Y} at time 3.
#'   \item Covariates \code{X_1}, \code{X_2}, and \code{X_3} at time 3.
#' }
#' If a subject is selected for LTFU, all their observations for **time points 4 and 5** are set to \code{NA}.
#'
#' **2. Intermittent Missingness:**
#' Variable-specific missingness is applied to \code{X_7} through \code{X_12} using logistic models that depend on the concurrent outcome \code{Y},
#' other covariates, and the previous value of the variable itself (autoregressive missingness).
#'
#' @examples
#' lt_data <- simulation_imputation_LTFU(NNY = TRUE, NNX = TRUE, n_subject = 10, seed = 42)
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS mvrnorm
#' @importFrom sn rmsn rmst
#' @importFrom dplyr if_else tibble
#' @importFrom arm invlogit
#' @importFrom mice md.pattern
#' @importFrom stats rnorm rbinom rt rlogis
simulation_imputation_LTFU = function(NNY = TRUE, NNX = TRUE, n_subject = 1000, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  n_obs_per_sub = 5
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


  Omega <- matrix(c(1.0, 0.3,
                    0.3, 1.0),
                  nrow = 2, byrow = TRUE)

  X = mvtnorm::rmvnorm(n_obs, rep(0, 6), sigma = AR(0.5, 6) * 0.5)

  X[,1:6] = apply(cbind(X[,1:6]), 2, function(X_v){
    Bi <- MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  X[,4:6] = (X[,4:6] >= 0)

  X_7 = cbind(-0.5 * (X[,1] + X[,4])^2 + 1 * X[,4])
  X_7 = apply(X_7, 2, function(X_v){
    if(NNX){
      Bi = dplyr::if_else(rbinom(n_subject, 1, 0.5) == 1, MASS::mvrnorm(n = n_subject, mu = c(1,1), Sigma = 0.25 * Omega), MASS::mvrnorm(n = n_subject, mu = c(-1,-1), Sigma = 0.25 * Omega))
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })


  X_8 = cbind(0.8 * X[,2] * X[,3] - X[,3])
  X_8 = apply(X_8, 2, function(X_v){
    if(NNX){
      Bi = dplyr::if_else(rbinom(n_subject, 1, 0.5) == 1, MASS::mvrnorm(n = n_subject, mu = c(2,2), Sigma = 0.25 * Omega), MASS::mvrnorm(n = n_subject, mu = c(-2,-2), Sigma = 0.25 * Omega))
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })
  X_9 = cbind(log(1e-80 + (X[,3] - X[,6] - 0.5)^2))
  X_9 = apply(X_9, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 3)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })

  if(NNX){
    X_7 = X_7 + dplyr::if_else(rbinom(n_obs, 1, 0.5) == 1, rnorm(n = n_obs, 2, sd = 0.5), rnorm(n = n_obs, -2, sd = 0.5))
    X_8 = X_8 + rt(n_obs, df = 3)
    X_9 = X_9 + rt(n_obs, df = 2)
  }else{
    X_7 = X_7 + rnorm(n_obs, 0, 1)
    X_8 = X_8 + rnorm(n_obs, 0, 1)
    X_9 = X_9 + rnorm(n_obs, 0, 1)
  }


  X_10 = cbind(0.5 * X[,2] - X[,5] - 0.8 * X_7)
  X_10 = apply(X_10, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 3)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })

  X_11 = cbind(X[,3] * X[,4] - 0.5 * X[,6] * X_9)
  X_11 = apply(X_11, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 2)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })


  X_12 = cbind(0.8 *X[,6] * X_8 - X[,1])
  X_12 = apply(X_12, 2, function(X_v){
    if(NNX){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.25 , alpha = c(0, 0), nu = 2)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 0.25 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })



  if(NNX){
    X_10 = (X_10 + rlogis(n_obs)) >= 0
    X_11 = (X_11 + rt(n_obs, df = 3)) >= 0
    X_12 = (X_12 + rt(n_obs, df = 2)) >= 0
  }else{
    X_10 = (X_10 + rnorm(n_obs, 0, 1)) >= 0
    X_11 = (X_11 + rnorm(n_obs, 0, 1)) >= 0
    X_12 = (X_12 + rnorm(n_obs, 0, 1)) >= 0
  }


  Y = cbind(10 + 1* X[,1] +1* X[,2] + 1*X[,3] + 1 *X[,4] + 1 *X[,5] + 1 *X[,6] + 1 *X_7 + 1 *X_8 + 1 *X_9 + 1*X_10 + 1*X_11 + 1*X_12)
  Y = apply(Y, 2, function(X_v){
    if(NNY){
      Bi = sn::rmst(n = n_subject, xi = c(0, 0), Omega = Omega * 0.5, alpha = c(0, 0), nu = 3)
    }else{
      Bi = MASS::mvrnorm(n = n_subject, mu = c(0,0), Sigma = 1 * Omega)
    }
    re = sapply(1:length(subject_id), function(x){
      Z[x,] %*% Bi[subject_id[x],]
    })
    X_v + re
  })



  if(NNY){
    Y = Y + rt(n_obs, df = 2) #
  }else{
    Y = Y + rnorm(n_obs, 0, 2)
  }




  data_E = dplyr::tibble(
    subject_id = subject_id,
    trajectory = trajectory[,1],
    X_1 = X[,1],
    X_2 = X[,2],
    X_3 = X[,3],
    X_4 = as.logical(X[,4]),
    X_5 = as.logical(X[,5]),
    X_6 = as.logical(X[,6]),
    X_7 = X_7[,1],
    X_8 = X_8[,1],
    X_9 = X_9[,1],
    X_10 = X_10[,1],
    X_11 = X_11[,1],
    X_12 = X_12[,1],
    Y = Y[,1]
  )


  data_O = dplyr::tibble(
    subject_id = subject_id,
    trajectory = trajectory[,1],
    X_1 = X[,1],
    X_2 = X[,2],
    X_3 = X[,3],
    X_4 = as.logical(X[,4]),
    X_5 = as.logical(X[,5]),
    X_6 = as.logical(X[,6]),
    X_7 = X_7[,1],
    X_8 = X_8[,1],
    X_9 = X_9[,1],
    X_10 = X_10[,1],
    X_11 = X_11[,1],
    X_12 = X_12[,1],
    Y = Y[,1]
  )

  R_7 = logical(n_obs)
  R_8 = logical(n_obs)
  R_9 = logical(n_obs)
  R_10 = logical(n_obs)
  R_11 = logical(n_obs)
  R_12 = logical(n_obs)
  #R_Y = logical(n_obs)

  alpha = -5 #-1
  for(j in 1:n_obs_per_sub[1]){
    if(j == 1){
      R_7[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == 1]) - 5.5 * data_O$X_6[data_O$trajectory == 1] - 3))
      R_8[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == 1]) - 3 * data_O$X_5[data_O$trajectory == 1] - 3))
      R_9[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == 1]) + 3 * data_O$X_4[data_O$trajectory == 1] -6))
      R_10[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == 1]) + 2 * data_O$X_3[data_O$trajectory == 1] - 4))
      R_11[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == 1]) - 5 * (data_O$X_2[data_O$trajectory == 1] + 1)^2 - 0))
      R_12[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == 1]) - 3 * (data_O$X_1[data_O$trajectory == 1] + 1)^2 - 2))
      #R_Y[data_O$trajectory == 1] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(alpha + -1.5 * (data_O$Y[data_O$trajectory == 1]) + 1*data_O$X_1[data_O$trajectory == 1] + 1*data_O$X_2[data_O$trajectory == 1] + 1*data_O$X_3[data_O$trajectory == 1] + 1*data_O$X_4[data_O$trajectory == 1] + 1*data_O$X_5[data_O$trajectory == 1] + 1*data_O$X_6[data_O$trajectory == 1] ))
    }else{
      R_7[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == j]) - 5.5 * data_O$X_6[data_O$trajectory == j] + 1 * data_O$X_7[data_O$trajectory == (j-1)] -3))
      R_8[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == j]) - 3 * data_O$X_5[data_O$trajectory == j] + 1 * data_O$X_8[data_O$trajectory == (j-1)] - 3))
      R_9[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == j]) +3 * data_O$X_4[data_O$trajectory == j] + 1 * data_O$X_9[data_O$trajectory == (j-1)] -6))
      R_10[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == j]) + 2 * data_O$X_3[data_O$trajectory == j] + 1 * data_O$X_10[data_O$trajectory == (j-1)]  -4))
      R_11[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == j]) - 5 * (data_O$X_2[data_O$trajectory == j] + 1)^2 + 1 * data_O$X_11[data_O$trajectory == (j-1)] -0))
      R_12[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(-1 * (data_O$Y[data_O$trajectory == j]) - 3 * (data_O$X_1[data_O$trajectory == j] + 1)^2 + 1 * data_O$X_12[data_O$trajectory == (j-1)] -2))


      #R_Y[data_O$trajectory == j] = rbinom(n = sum(data_O$trajectory == 1), size = 1, prob = arm::invlogit(alpha + -1.5 * (data_O$Y[data_O$trajectory == j]) + 1*data_O$X_1[data_O$trajectory == j] + 1*data_O$X_2[data_O$trajectory == j] + 1*data_O$X_3[data_O$trajectory == j] + 1*data_O$X_4[data_O$trajectory == j] + 1*data_O$X_5[data_O$trajectory == j] + 1*data_O$X_6[data_O$trajectory == j] - 1 * data_O$Y[data_O$trajectory == (j-1)]))

      R_7[data_O$trajectory == j] = R_7[data_O$trajectory == j] * (1 - R_7[data_O$trajectory == (j-1)])
      R_8[data_O$trajectory == j] = R_8[data_O$trajectory == j] * (1 - R_8[data_O$trajectory == (j-1)])
      R_9[data_O$trajectory == j] = R_9[data_O$trajectory == j] * (1 - R_9[data_O$trajectory == (j-1)])
      R_10[data_O$trajectory == j] = R_10[data_O$trajectory == j] * (1 - R_10[data_O$trajectory == (j-1)])
      R_11[data_O$trajectory == j] = R_11[data_O$trajectory == j] * (1 - R_11[data_O$trajectory == (j-1)])
      R_12[data_O$trajectory == j] = R_12[data_O$trajectory == j] * (1 - R_12[data_O$trajectory == (j-1)])

      #R_Y[data_O$trajectory == j] = R_Y[data_O$trajectory == j] * (1 - R_Y[data_O$trajectory == (j-1)])

    }
  }

  data_M = dplyr::tibble(
    subject_id = subject_id,
    trajectory = trajectory[,1],
    X_1 = X[,1],
    X_2 = X[,2],
    X_3 = X[,3],
    X_4 = as.logical(X[,4]),
    X_5 = as.logical(X[,5]),
    X_6 = as.logical(X[,6]),
    X_7 = X_7[,1],
    X_8 = X_8[,1],
    X_9 = X_9[,1],
    X_10 = X_10[,1],
    X_11 = X_11[,1],
    X_12 = X_12[,1],
    Y = Y[,1]
  )
  data_M$X_7[R_7 == 1] = NA
  data_M$X_8[R_8 == 1] = NA
  data_M$X_9[R_9 == 1] = NA
  data_M$X_10[R_10 == 1] = NA
  data_M$X_11[R_11 == 1] = NA
  data_M$X_12[R_12 == 1] = NA
  #data_M$Y[R_Y == 1] = NA


  R_ltfu = rbinom(n = n_subject, size = 1, prob = arm::invlogit(-8.5 + 1 * (data_O$Y[data_O$trajectory == 3]) + 2*data_O$X_1[data_O$trajectory == 3] + 2*data_O$X_2[data_O$trajectory == 3] + 2*data_O$X_3[data_O$trajectory == 3]))
  for (i in 1:n_subject) {
    if(R_ltfu[i] == 1){
      data_M[data_M$subject_id == i & data_M$trajectory %in% c(4,5), 3:(dim(data_M)[2])] = NA
    }
  }
  print(mean(R_ltfu))
  #mice::md.pattern(data_M[,3:dim(data_M)[2]])
  return(list(data_E = data_E, data_M = data_M, data_O = data_O, Z = Z, pair = mice::md.pattern(data_M[,3:dim(data_M)[2]], plot = FALSE)))
}
