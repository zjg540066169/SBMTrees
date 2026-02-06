#' @title Initialize Missing Values using LOCF and NOCB
#'
#' @description Imputes missing values in longitudinal data using a hierarchical three-step strategy
#' to ensure complete data for model initialization. The process prioritizes within-subject information
#' using Last Observation Carried Forward (LOCF) and Next Observation Carried Backward (NOCB),
#' falling back to cross-sectional summary statistics (mean or mode) only when a subject has absolutely
#' no observed data for a specific variable.
#'
#' @param X A data.frame or matrix containing the variables to be imputed. Columns correspond to variables.
#' @param subject_id A vector of subject identifiers with length equal to \code{nrow(X)}.
#' @param is_binary A vector of length \code{ncol(X)} indicating the type of each variable.
#' Values can be \code{TRUE}/\code{1} (for binary variables) or \code{FALSE}/\code{0} (for continuous variables).
#'
#' @return A data.frame with the same dimensions as \code{X} but with all missing values imputed.
#'
#' @details
#' **Pre-requisite:** The rows of \code{X} must be ordered by time within each subject prior to calling this function.
#'
#' The imputation proceeds in three specific stages:
#' \enumerate{
#'   \item **Subject-wise LOCF**: For each subject, missing values are filled using the immediately preceding observed value (forward fill). This handles gaps in the middle or end of a subject's timeline.
#'   \item **Subject-wise NOCB**: For each subject, any remaining missing values (typically at the start of the timeline, before the first observation) are filled using the next available observed value (backward fill).
#'   \item **Global Fallback**: If a subject has \emph{no} observed data for a specific variable (i.e., the entire column is \code{NA} for that \code{subject_id}), the function imputes these values using the global statistics calculated from the rest of the population:
#'   \itemize{
#'      \item **Continuous variables**: Imputed with the global mean.
#'      \item **Binary variables**: Imputed with the global mode (ties default to 0).
#'   }
#' }
#'
#' @examples
#' # Create a toy dataset with missing values
#' X <- data.frame(
#'   cont = c(NA, 5, NA,   NA, NA, NA),  # Subj 1: Gap/Lead/Trail, Subj 2: All NA
#'   bin  = c(0, NA, 1,    1, 1, 0)      # Subj 1: Gap,            Subj 2: Complete
#' )
#' subject_id <- c(1, 1, 1, 2, 2, 2)
#' is_binary <- c(FALSE, TRUE)
#'
#' # Run imputation
#' X_imputed <- apply_locf_nocb(X, subject_id, is_binary)
#'
#' @export
apply_locf_nocb <- function(X, subject_id, is_binary) {
  # X: data.frame or matrix (columns = variables to impute)
  # subject_id: vector same length as nrow(X)
  # is_binary: 0/1 (or TRUE/FALSE) vector of length ncol(X); 1= binary, 0 = continuous
  # NOTE: rows must already be ordered by time within subject

  # --- helpers ---
  locf <- function(x) {
    if (length(x) <= 1) return(x)
    for (i in 2:length(x)) if (is.na(x[i])) x[i] <- x[i - 1]
    x
  }
  nocb <- function(x) {
    if (length(x) <= 1) return(x)
    for (i in (length(x) - 1):1) if (is.na(x[i])) x[i] <- x[i + 1]
    x
  }
  mode01 <- function(v) {
    v <- as.numeric(v)
    v <- v[!is.na(v)]
    if (!length(v)) return(0)           # fallback if all NA
    n1 <- sum(v == 1); n0 <- sum(v == 0)
    ifelse(n1 > n0, 1, 0)               # ties -> 0
  }

  # --- checks ---
  X <- as.data.frame(X)
  if (ncol(X) == 0) return(X)
  if (length(is_binary) != ncol(X))
    stop("`is_binary` must have length equal to ncol(X).")
  is_binary <- as.logical(is_binary)

  # --- LOCF/NOCB per subject ---
  uniq_ids <- unique(subject_id)
  for (id in uniq_ids) {
    rows <- which(subject_id == id)
    if (length(rows) < 2) next
    # LOCF across all columns for this subject
    X[rows, ] <- lapply(X[rows, ], locf)
  }
  for (id in uniq_ids) {
    rows <- which(subject_id == id)
    if (length(rows) < 2) next
    # NOCB across all columns for this subject
    X[rows, ] <- lapply(X[rows, ], nocb)
  }

  # --- final columnwise fill (mean for continuous, mode for binary) ---
  for (j in seq_len(ncol(X))) {
    col <- X[[j]]
    if (!anyNA(col)) next

    if (is_binary[j]) {
      fill_val <- mode01(col)
      col[is.na(col)] <- fill_val
    } else {
      col_num <- suppressWarnings(as.numeric(col))
      # mean; if all NA, fallback to 0 (change if you prefer another default)
      m <- mean(col_num, na.rm = TRUE)
      if (is.nan(m)) m <- 0
      col[is.na(col)] <- m
    }
    X[[j]] <- col
  }

  X
}




#' @title Longitudinal Sequential Imputation for Longitudinal Missing Data
#'
#' @description Implements sequential imputation for missing covariates and outcomes in longitudinal data.
#' The function uses a Bayesian non-parametric framework with mixed-effects models to handle both
#' normal and non-normal random effects and errors. It sequentially imputes missing values by constructing
#' univariate models in a fixed order, initializing with LOCF/NOCB, and ensuring consistency with a valid joint distribution.
#'
#' @param X A matrix of missing covariates.
#' @param Y A vector of missing outcomes (numeric or logical).
#' @param Z A matrix of complete random predictors. Default: \code{NULL}.
#' @param subject_id A vector of subject IDs corresponding to the rows of \code{X} and \code{Y}. Can be integer, factor, or character.
#' @param type A vector indicating whether each covariate in \code{X} is binary (1) or continuous (0).
#' @param binary_outcome A logical value indicating whether the outcome \code{Y} is binary. Default: \code{FALSE}.
#' @param model A character vector specifying the imputation model for the covariates. Options are \code{"BMTrees"} (default),
#' \code{"BMTrees_R"} (residual DP), \code{"BMTrees_RE"} (random effect DP), and \code{"mixedBART"}.
#' @param outcome_model A character vector specifying the model used for the outcome. Options are \code{"BMTrees"} (default) or
#' \code{"BMLM"} (Bayesian Mixed Linear Model). If \code{"BMLM"} is selected, posterior estimates for beta and sigma are returned.
#' @param nburn An integer specifying the number of burn-in iterations. Default: \code{0}.
#' @param npost An integer specifying the number of sampling iterations. Default: \code{3}.
#' @param skip An integer specifying the interval for keeping samples in the sampling phase. Default: \code{1}.
#' @param verbose A logical value indicating whether to display progress and MCMC information. Default: \code{TRUE}.
#' @param seed A random seed for reproducibility. Default: \code{NULL}.
#' @param tol A small numerical tolerance to prevent numerical overflow or underflow in the model. Default: \code{1e-20}.
#' @param k A numeric value for the BART prior parameter controlling the standard deviation of the terminal node values. Default: \code{2.0}.
#' @param ntrees An integer specifying the number of trees in BART. Default: \code{200}.
#' @param reordering A logical value indicating whether to apply a reordering strategy for sorting covariates based on missingness. Default: \code{TRUE}.
#' @param pi_DP A value between 0 and 1 for calculating the empirical prior in the DP prior. Default: \code{0.99}.
#'
#' @return A list containing:
#' \item{imputed_data}{A three-dimensional array of imputed data with dimensions \code{(npost / skip, N, p + 1)}, where \code{N} is the number of observations and \code{p} is the number of covariates. The last column represents the outcome \code{Y}.}
#' \item{posterior_sigma}{(Only if \code{outcome_model = "BMLM"}) A vector of posterior samples for the error standard deviation.}
#' \item{posterior_beta}{(Only if \code{outcome_model = "BMLM"}) A matrix of posterior samples for the regression coefficients.}
#'
#' @details The function builds on the Bayesian Trees Mixed-Effects Model (BMTrees), which extends Mixed-Effects
#' BART by using centralized Dirichlet Process Normal Mixture priors. This framework handles non-normal
#' random effects and errors, addresses model misspecification, and captures complex relationships.
#'
#' The algorithm initializes missing values using Last Observation Carried Forward (LOCF) and Next Observation Carried Backward (NOCB)
#' before starting the MCMC sequential imputation process.
#'
#' @examples
#' data <- simulation_imputation(NNY = TRUE, NNX = TRUE, n_subject = 10, seed = 123)
#' BMTrees <- sequential_imputation(X = data$data_M[,3:5], Y = data$data_M$Y, Z = data$Z,
#'   subject_id = data$data_M$subject_id, type = c(0, 0, 0),
#'   outcome_model = "BMLM", binary_outcome = FALSE, model = "BMTrees", nburn = 0,
#'   npost = 1, skip = 1, verbose = FALSE, seed = 123)
#'
#' # Access imputed data
#' dim(BMTrees$imputed_data)
#' @rdname sequential_imputation
#' @note
#' This function utilizes modified C++ code originally derived from the
#' BART3 package (Bayesian Additive Regression Trees). The original package
#' was developed by Rodney Sparapani and is licensed
#' under GPL-2. Modifications were made by Jungang Zou, 2024.
#' @references
#' For more information about the original BART3 package, see:
#' https://github.com/rsparapa/bnptools/tree/master/BART3
#' @export
#' @useDynLib SBMTrees, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats var
sequential_imputation <- function(X, Y,  Z = NULL, subject_id, type, binary_outcome = FALSE, model = c("BMTrees", "BMTrees_R", "BMTrees_RE", "mixedBART"), outcome_model = c("BMTrees", "BMLM"), nburn = 0L, npost = 3L, skip = 1L, verbose = TRUE, seed = NULL, tol = 1e-20, k = 2.0, ntrees = 200, reordering = TRUE, pi_DP = 0.99) {
  model = match.arg(model)
  if(is.null(dim(X))){
    stop("More than one covariate is needed!")
  }
  if(!is.null(seed)){
    set.seed(seed)
  }

  outcome_BMTrees = ifelse(outcome_model == "BMTrees", TRUE, FALSE)
  constant = any(apply(X, 2, function(x) stats::var(x, na.rm = TRUE)) == 0)
  trajectory = c(sapply(table(subject_id), function(i) 1:i))
  subject_id = as.factor(subject_id)
  p = dim(X)[2]
  if(!constant){
    #if(sum(mis_num == -Inf) == 0 | correct_prob > 0){
    X = cbind("intercept" = 1, X)
    type = c(0, type)
  }


  mis_num = colSums(is.na(X))
  if(reordering == TRUE){
    mis_num[mis_num ==0] = -Inf
    mis_num[mis_num > 0 & type == 1] = mis_num[mis_num > 0 & type == 1] - base::max(mis_num)
    mis_order = order(mis_num)
    re_order = seq(ncol(is.na(X)))
    type = type[mis_order]
    X[,re_order] = X[,mis_order]
    colnames(X) = colnames(X)[mis_order]
    message("reordering: new covariates order is ", paste(colnames(X), collapse = " "), "\n")
  }


  R = is.na(X)
  R = cbind(R, "Y" = is.na(Y))
  message("Start to initialize imputed missing data by LOCF and NOCB. ")
  X_locf_nocb <- apply_locf_nocb(cbind(X, Y), subject_id, c(type, 0))
  #return(X_locf_nocb)
  Y = X_locf_nocb[,dim(X_locf_nocb)[2]]
  X = X_locf_nocb[,1:(dim(X_locf_nocb)[2] - 1)]
  #print(apply(X[,type == 1], 2, unique))


  message("Completed.\n")


  message("Start to impute using Longitudinal Sequential Imputation with: ")

  if(model == "BMTrees_R"){
    message("BMTrees_R\n")
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), outcome_BMTrees = outcome_BMTrees, binary_outcome = binary_outcome, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = TRUE, CDP_re = FALSE, seed = seed, ncores = 0, ntrees = ntrees, k = k, pi_CDP = pi_DP, correct_prob = 0)
  }
  else if(model == "BMTrees_RE"){
    message("BMTrees_RE\n")
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), outcome_BMTrees = outcome_BMTrees, binary_outcome = binary_outcome, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = FALSE, CDP_re = TRUE, seed = seed, ncores = 0, ntrees = ntrees, k = k, pi_CDP = pi_DP, correct_prob = 0)
  }
  else if(model == "BMTrees"){
    message("BMTrees\n")
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), outcome_BMTrees = outcome_BMTrees, binary_outcome = binary_outcome, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = TRUE, CDP_re = TRUE, seed = seed, ncores = 0, ntrees = ntrees, k = k, pi_CDP = pi_DP, correct_prob = 0)
  }
  else if(model == "mixedBART"){
    message("mixedBART\n")
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), outcome_BMTrees = outcome_BMTrees, binary_outcome = binary_outcome, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = FALSE, CDP_re = FALSE, seed = seed, ncores = 0,  ntrees = ntrees, k = k, pi_CDP = pi_DP, correct_prob = 0)
  }
  else{
    message("BMTrees\n")
    imputation_X_DP = sequential_imputation_cpp(as.matrix(X), as.numeric(Y), as.logical(type), as.matrix(Z), as.character(subject_id), as.matrix(R), outcome_BMTrees = outcome_BMTrees, binary_outcome = binary_outcome, nburn = nburn, npost = npost, skip = skip, verbose = verbose, CDP_residual = TRUE, CDP_re = TRUE, seed = seed, ncores = 0, ntrees = ntrees, k = k, pi_CDP = pi_DP, correct_prob = 0)
  }
  #return(imputation_X_DP)
  posterior_sigma = imputation_X_DP$posterior_sigma
  posterior_beta = imputation_X_DP$posterior_beta
  #return(imputation_X_DP)
  imputation_Y = t(do.call(cbind, imputation_X_DP$imputation_Y_DP))
  if(reordering == TRUE){
    imputation_X_DP = lapply(imputation_X_DP$imputation_X_DP, function(x){
      #if(sum(mis_num == -Inf) == 0){
      x[,mis_order] = x[,re_order]
      colnames(x)[mis_order] = colnames(x)[re_order]
      if(!constant){
        x = x[,-1]
      }
      x
    })

    posterior_beta[, mis_order] = posterior_beta[, re_order]
    colnames(posterior_beta)[mis_order] = colnames(posterior_beta)[re_order]
  }else{
    imputation_X_DP = lapply(imputation_X_DP$imputation_X_DP, function(x){
      #if(sum(mis_num == -Inf) == 0){
      if(!constant){
        x = x[,-1]
      }
      x
    })
  }

  imputation_X = array(NA, dim = c(length(imputation_X_DP), dim(imputation_X_DP[[1]])[1], dim(imputation_X_DP[[1]])[2] + 1))
  for (i in 1:length(imputation_X_DP)) {
    imputation_X[i,,] = cbind(imputation_X_DP[[i]], imputation_Y[i,])
  }
  message("\n")
  message("Finish imputation with ", length(imputation_X_DP), " imputed sets\n")
  if(outcome_BMTrees){
    return(list(imputed_data = imputation_X))
  }else{
    return(list(imputed_data = imputation_X, posterior_sigma = posterior_sigma, posterior_beta = posterior_beta))
  }

}
