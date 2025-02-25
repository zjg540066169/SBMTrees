\name{SBMTrees-package}
\alias{SBMTrees-package}
\alias{SBMTrees}
\docType{package}
\title{
  Sequential Imputation with Bayesian Trees Mixed-Effects Models
}
\description{
  The SBMTrees package implements a Bayesian non-parametric framework for imputing missing covariates and outcomes in longitudinal data under the Missing at Random (MAR) assumption. Its core model, the Bayesian Trees Mixed-Effects Model (BMTrees), extends Mixed-Effects BART by employing centralized Dirichlet Process (CDP) Normal Mixture priors. This allows handling non-normal random effects and errors, addressing model misspecification, and capturing complex relationships.
}
\details{
  SBMTrees offers tools for predicting and imputing missing values in longitudinal data using Bayesian Trees Mixed-Effects Models. The package supports various semiparametric variants, including BMTrees_R and BMTrees_RE, and integrates mixedBART as a baseline model. Key functionalities include:

  - \code{BMTrees_prediction}: Predicts longitudinal outcomes based on mixed-effects models.
  
  - \code{sequential_imputation}: Imputes missing covariates and outcomes sequentially in longitudinal datasets.

  The package supports flexibility in specifying priors for random effects and errors, making it suitable for diverse longitudinal data settings. Core computations leverage efficient Gibbs samplers implemented in C++.
  
  This package modifies and extends C++ code originally derived from the 
  BART3 package, developed by Rodney Sparapani, 
  which is licensed under the GNU General Public License version 2 (GPL-2).

  The modified code is redistributed in accordance with the GPL-2 license. 
  For more details on the modifications, see the package's documentation.
}
\author{
  Jungang Zou

  Maintainer: Jungang Zou <jungang.zou@gmail.com>
}
\keyword{SBMTrees, longitudinal missing data, Bayesian non-parametric methods, sequential imputation, mixed-effects models, multiple imputation}
\seealso{
  \code{\link[=BMTrees_prediction]{BMTrees_prediction}}, \code{\link[=sequential_imputation]{sequential_imputation}}
}
\note{
  This package and all associated documentation are licensed under the GNU General Public License version 2 (GPL-2). See the \code{LICENSE} file for the full text of the license.
}
\references{
  BART3 package: https://github.com/rsparapa/bnptools/tree/master, originally developed by Rodney Sparapani. 
}
\examples{
  \dontrun{
     ## Example of predicting longitudinal outcomes
     data <- simulation_prediction(n_subject = 800, seed = 123, nonlinear = TRUE,
             nonrandeff = TRUE, nonresidual = TRUE)
     X_train <- data$X_train
     Y_train <- data$Y_train
     Z_train <- data$Z_train
     subject_id_train <- data$subject_id_train

     X_test <- data$X_test
     Z_test <- data$Z_test
     subject_id_test <- data$subject_id_test

     model <- BMTrees_prediction(X_train, Y_train, Z_train, subject_id_train, 
              X_test, Z_test, subject_id_test, model = "BMTrees")
     model$post_predictive_y_test
     
     
     
     data2 = simulation_imputation(n_subject = 800, seed = 123, nonrandeff = TRUE, 
             nonresidual = TRUE, alligned = F) 
     X_mis = data2$X_mis # get missing covariates
     Y_mis = data2$Y_mis # get missing outcomes
     Z = data2$Z # get random predictors
     subject_id = data2$subject_id  # get subject id

     model2 = sequential_imputation(X_mis, Y_mis, Z, subject_id, rep(0, 9), F, 
             model = "BMTrees", nburn = 3000L, npost = 4000L, skip = 200L, 
             verbose = TRUE, seed = 123)
     model2$imputed_data
     model2$imputed_data[,,10]
  }
}

