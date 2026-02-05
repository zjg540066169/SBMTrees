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
  Jungang Zou <jungang.zou@gmail.com>
}
\concept{SBMTrees}
\concept{longitudinal missing data}
\concept{Bayesian non-parametric methods}
\concept{sequential imputation}
\seealso{
  \code{\link[=BMTrees_prediction]{BMTrees_prediction}}, \code{\link[=sequential_imputation]{sequential_imputation}}
}
\note{
  This package and all associated documentation are licensed under the GNU General Public License version 2 (GPL-2). See the \code{LICENSE} file for the full text of the license.
}
\references{
  BART3 package: https://github.com/rsparapa/bnptools/tree/master, originally developed by Rodney Sparapani.
}

