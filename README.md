# SBMtrees 

The R package **SBMtrees** (Sequential imputation with Bayesian Trees Mixed-Effects models) implements a Bayesian non-parametric framework for imputing missing covariates and outcomes in longitudinal data under the Missing at Random (MAR) assumption. Its core model, the Bayesian Trees Mixed-Effects Model (BMTrees), extends Mixed-Effects BART by employing centralized Dirichlet Process (CDP) Normal Mixture priors, allowing it to handle non-normal random effects and errors, address model misspecification, and capture complex relationships. The package also includes two semiparametric variants, BMTrees_R and BMTrees_RE. Built on BMTrees, the longitudinal sequential imputation framework employs a Metropolis-Hastings (M-H) MCMC method to sequentially impute missing values by constructing univariate models in a fixed order, ensuring both simplicity and consistency with a valid joint distribution.

For more details on these models and their applications, please consult the following paper: "Nonparametric Bayesian Additive Regression Trees for Predicting and Handling Missing Covariates and Outcomes in Longitudinal Data".


## Installation
This package is based on `Rcpp`, `RcppArmadillo`, and `RcppDist`, please make sure these three packages can be installed.

Right now, the package has not been uploaded to Rcran yet, so please install this package from Github:
```
require("devtools")
install_github("https://github.com/zjg540066169/SBMtrees")
library(SBMTrees)
```

## Models
This package is based on the mixed-effects model for longitudinal data: <img src="https://latex.codecogs.com/gif.latex?Y_{ij}=BART(X_{ij})+Z_{ij}b_i+\epsilon_{ij}" /> 

Different models impose different prior distributions on <img src="https://latex.codecogs.com/gif.latex?b_i"/> and <img src="https://latex.codecogs.com/gif.latex?\epsilon_{ij}" />. We also include the existing model Mixed-Effects BART (mixedBART) in this package.

<table>
   <tr>
      <th align="center">Models</th>
      <th align="center">Prior on random effects <img src="https://latex.codecogs.com/gif.latex?b_i" /> </th>
      <th align="center">Prior on random errors <img src="https://latex.codecogs.com/gif.latex?\epsilon_{ij}" /> </th>
   </tr>
   <tr>
      <td style="text-align:center" align="center" rowspan="1" colspan="1">BMTrees</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">CDP Multivariate Normal Mixture</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">CDP Normal Mixture</td>
   </tr>
   <tr>
      <td style="text-align:center" align="center" rowspan="1" colspan="1">BMTrees_R</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">Multivariate Normal</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">CDP Normal Mixture</td>
   </tr>
   <tr>
      <td style="text-align:center" align="center" rowspan="1" colspan="1">BMTrees_RE</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">CDP Multivariate Normal Mixture</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">Normal</td>
   </tr>
   <tr>
      <td style="text-align:center" align="center" rowspan="1" colspan="1">mixedBART</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">Multivariate Normal</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">Normal</td>
   </tr>
</table>

The inference is done with posterior samples by Gibbs samplers in C++. 


## Usage
There are two main functions in this package. `BMTrees_prediction` is employed to estimate and predict longitudinal outcomes. `sequential_imputation` is used to multiply-impute longitudinal missing covariates and outcomes. 

### Prediction
We first generate a data with some individuals, each has 6 follow-up time points. As described in paper, we can specify if the linear/nonlinear associations, normal/non-normal random effects and random error. For each subject, one to three time points were randomly chosen to form the testing dataset, while the remaining time points constituted the training dataset. The testing dataset accounted for roughly 40% of the total data.

This can be achieved by running the function `simulation_prediction(n_subject = 800, seed, nonlinear, nonrandeff, nonresidual)`. The five parameters are:
* n_subject (integer): positive number to indicate the total number of individuals. The total number of observations is n_subject * 6. In training set, there are about n_subject * 3.6 observations, while the number of observations is n_subject * 2.4 in testing set.
* seed: random seed for replicates. The population data is fixed. The random seed ensures different training-testing splits for the fixed population data.
* nonlinear (binary): indicate if the outcome model is nonlinear or not. See more details the paper.
* nonrandeff (binary): indicate if the random effects is non-normal or not. See more details the paper.
* nonresidual (binary): indicate if the random errors is non-normal or not. See more details the paper.

Here is an example:
```
library(SBMTrees)
data = simulation_prediction(n_subject = 800, seed = 123, nonlinear = TRUE, nonrandeff = TRUE, nonresidual = TRUE) 
X_train = data$X_train # get predictors in training set
Y_train = data$Y_train # get outcomes in training set
Z_train = data$Z_train # get random predictors in training set
subject_id_train = data$subject_id_train # get subject id in training set

X_test = data$X_test # get predictors in testing set
Y_test = data$Y_test # get outcomes in testing set
Z_test = data$Z_test # get random predictors in testing set
subject_id_test = data$subject_id_test # get subject id in testing set
```

After we get data, we can run the prediction model based on function `BMTrees_prediction(X_train, Y_train, Z_train, subject_id_train, X_test, Z_test, subject_id_test, model = c("BMTrees", "BMTrees_R", "BMTrees_RE", "mixedBART"), binary = FALSE, nburn = 3000L, npost = 4000L, skip = 1L, verbose = TRUE, seed, tol = 1e-20, resample = 5, ntrees = 200, pi_CDP = 0.99)`. These parameters are:
* X_train (matrix): covariates matrix in training set.
* Y_train (numeric or logical): outcomes in training set.
* Z_train (matrix): matrix of random predictors in training set.
* subject_id_train (character): subject id in training set.
* X_test (matrix): covariates matrix in testing set.
* Y_test (numeric or logical): outcomes in testing set.
* Z_test (matrix): matrix of random predictors in testing set.
* subject_id_test (character): subject id in testing set.
* model (character): choose the predictive model.
* binary (logical): indicate if the outcome is binary or continuous.
* nburn (integer): number of burn-in phase.
* npost (integer): number of sampling phase.
* skip (integer): keep one sample for every skip draws in sampling phase.
* verbose (logical): if output the MCMC information. If verbose is FALSE, we will run a progress-bar.
* seed: random seed for prediction.
* tol: a small numerical tolenrance in model to prevent potential numerical overflow and underflow.
* resample (integer): the number of resamping number for CDP prior. See details in paper.
* ntrees (integer): number of trees.
* pi_CDP: a value between 0 and 1, to calculate the empirical prior in CDP prior. See details in paper.

Here is an example to run the predictive model.
```
model = BMTrees_prediction(X_train, Y_train, Z_train, subject_id_train, X_test, Z_test, subject_id_test, model = "BMTrees", binary = FALSE, nburn = 3000L, npost = 4000L, skip = 1L, verbose = TRUE, seed = 123)
model$posterior_Y_test
model$posterior_sigma
```
The users can get the posterior predictive samples for Y_test and posterior draws of other parameters.


### Imputation
For imputation, we first generate a dataset comprising individuals, each with six follow-up time points. As described in the paper, we can specify whether the random effects and random errors follow normal or non-normal distributions. Different missingness mechanisms are applied to create MAR missing values, resulting in approximately 35% of the observations having missing data.

The data with missingness is generated by running the function `simulation_imputation = function(n_subject = 800, seed, nonrandeff, nonresidual, alligned = F)`. The five parameters are:
* n_subject (integer): positive number to indicate the total number of individuals. The total number of observations is n_subject * 6.
* seed: random seed for replicates. The random seed ensures different covariates, outcomes and missingness are generated.
* nonrandeff (binary): indicate if the random effects is non-normal or not. See more details the paper.
* nonresidual (binary): indicate if the random errors is non-normal or not. See more details the paper.
* alligned (binary): if alligned is TRUE, the returned missing covariates order is (X7, X8, X9). If alligned is FALSE, we will shuffle the order as (X9, X8, X7).
  
Here is an example:
```
library(SBMTrees)
data = simulation_imputation(n_subject = 800, seed = 123, nonrandeff = TRUE, nonresidual = TRUE, alligned = F) 
X_mis = data$X_mis # get missing covariates
Y_mis = data$Y_mis # get missing outcomes
Z = data$Z # get random predictors
subject_id = data$subject_id  # get subject id
time = data$time # get time point

X_O = data$X_O # get the original covariates matrix which is complete, for evaluation
Y_O = data$Y_O # get the original outcome which is complete, for evaluation
```

After we get data, we can run the prediction model based on function `sequential_imputation(X_mis, Y_mis, Z, subject_id, type, binary_outcome, model = c("BMTrees", "BMTrees_R", "BMTrees_RE", "mixedBART"), nburn = 3000L, npost = 4000L, skip = 200L, verbose = TRUE, seed, tol = 1e-20, resample = 5, ntrees = 200, reordering = T, pi_CDP = 0.99)`. These parameters are:
* X_mis (matrix): missing covariates matrix.
* Y_mis (numeric or logical): missing outcomes.
* Z (matrix): matrix of random predictors, which is complete
* subject_id (character): subject id.
* type (logical): a vector indicate if the covariates are binary (1) or continuous (0). The length should be same as the number of covariates in X_mis.
* binary_outcome (logical): a logical value indicate if the outcome are binary (1) or continuous (0).
* model (character): choose the imputation model.
* nburn (integer): number of burn-in phase.
* npost (integer): number of sampling phase.
* skip (integer): keep one sample for every skip draws in sampling phase.
* verbose (logical): if output the MCMC information. If verbose is FALSE, we will run a progress-bar.
* seed: random seed for imputation.
* tol: a small numerical tolenrance in model to prevent potential numerical overflow and underflow.
* resample (integer): the number of resamping number for CDP prior. See details in paper.
* ntrees (integer): number of trees.
* reordering (logical): a logical value to indicate if the reordering strategy will be applied to sort covariates order. See details in paer.
* pi_CDP: a value between 0 and 1, to calculate the empirical prior in CDP prior. See details in paper.

Here is an example to run the predictive model.
```
model = sequential_imputation(X_mis, Y_mis, Z, subject_id, rep(0, 9), F, model = "BMTrees", nburn = 3000L, npost = 4000L, skip = 200L, verbose = TRUE, seed = 123)
model$imputed_data
model$imputed_data[,,10] # get imputed outcomes.
```
The returned `imputed_data` is a three-array, whose dimension is (npost / skip, N, p + 1). N is the number of observations. p is the number of covariates.


## Attribution

This package includes code derived from the [BART3](https://github.com/rsparapa/bnptools/tree/master) package, originally developed by Rodney Sparapani. 

The original source code, licensed under the [GNU General Public License version 2 (GPL-2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html), has been modified as follows:
- We include part of the C++ code in BART3, primarily about functions about `wbart` and `cpwart`. We also modify some files to make sure our package can be successfully compiled.
- Modifications were made by Jungang Zou, 2024.

### Licensing

- The original BART3 package is licensed under the GNU General Public License version 2 (GPL-2).
- This package, as a derived work, is also licensed under the GNU General Public License version 2 (GPL-2) to comply with the licensing terms.

For more details, please see the `LICENSE` file included in this repository.

## Contact

If you find there is any bug, please contact me: jungang.zou@gmail.com.
