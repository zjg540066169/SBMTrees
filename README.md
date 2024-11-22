# SBMtrees 

The R package **SBMtrees** (Sequential imputation with Bayesian Trees Mixed-Effects models) implements a Bayesian non-parametric framework for imputing missing covariates and outcomes in longitudinal data under the Missing at Random (MAR) assumption. Its core model, the Bayesian Trees Mixed-Effects Model (BMTrees), extends Mixed-Effects BART by employing centralized Dirichlet Process (CDP) Normal Mixture priors, allowing it to handle non-normal random effects and errors, address model misspecification, and capture complex relationships. The package also includes two semiparametric variants, BMTrees_R and BMTrees_RE. Built on BMTrees, the longitudinal sequential imputation framework employs a Metropolis-Hastings (M-H) MCMC method to sequentially impute missing values by constructing univariate models in a fixed order, ensuring both simplicity and consistency with a valid joint distribution.

For more details on these models and their applications, please consult the following paper: "Nonparametric Bayesian Additive Regression Trees for Predicting and Handling Missing Covariates and Outcomes in Longitudinal Data".


## Installation
This package is based on `Rcpp`, `RcppArmadillo`, and `RcppDist`, please make sure these three packages can be installed.

Right now, the package has not been uploaded to Rcran yet, so please install this package from Github:
```
require("devtools")
install_github("https://github.com/zjg540066169/SBMtrees")
library(SBMtrees)
```

## Models
This package is based on the mixed-effects model for longitudinal data: <img src="https://latex.codecogs.com/gif.latex?Y_{ij}=BART(X_{ij})+Z_{ij}b_i+\epsilon_{ij}" /> 

Different models impose different prior distributions on <img src="https://latex.codecogs.com/gif.latex?b_i"/> and <img src="https://latex.codecogs.com/gif.latex?\epsilon_{ij}" />. 

<table>
   <tr>
      <th align="center">Models</th>
      <th align="center">Prior on random effects <img src="https://latex.codecogs.com/gif.latex?b_i" /> </th>
      <th align="center">Prior on random errors <img src="https://latex.codecogs.com/gif.latex?\epsilon_{ij}" /> </th>
      <th align="center">Other priors</th>
   </tr>
   <tr>
      <td style="text-align:center" align="center" rowspan="1" colspan="1">BMTrees</td>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">CDP Normal Mixture <img src="https://latex.codecogs.com/gif.latex?MVN(\bmu_i, \mathbf{\Sigma})" /></td>
      <td style="text-align:center" align="center" colspan="1">r</td>
      <td style="text-align:center" align="center" colspan="1">2</td>
   </tr>
   <tr>
    <td style="text-align:center" align="center" colspan="1">s</td>
      <td style="text-align:center"  align="center" colspan="1">15</td>
   </tr>
   <tr>
      <td style="text-align:center" align="center" colspan="1" rowspan="1">ARD(Ridge)</td>
      <td style="text-align:center" align="center" colspan="2">No hyper-parameters</td>
   </tr>
   <tr>
    <td style="text-align:center" align="center" colspan="1">Horseshoe</td>
      <td style="text-align:center" align="center" colspan="2">No hyper-parameters</td>
   </tr>
   <tr>
      <td style="text-align:center" align="center" colspan="1" rowspan="5">Discrete Mixture Model</td>
      <td style="text-align:center" align="center" rowspan="2">Spike-Ridge</td>
      <td style="text-align:center" align="center" rowspan="1">p0</td>
      <td style="text-align:center" align="center" rowspan="1">0.5</td>
   <tr>
     <td style="text-align:center" align="center" colspan="1">v0</td>
      <td style="text-align:center" align="center" colspan="1">4</td>
   </tr>
   <tr>
      <td style="text-align:center" align="center" colspan="1" rowspan="3">Spike-Laplace</td>
      <td style="text-align:center" align="center" colspan="1">lambda</td>
      <td style="text-align:center" align="center" colspan="1">6/11</td>
   </tr>
   <tr>
      <td style="text-align:center" align="center" colspan="1">a</td>
      <td style="text-align:center" align="center" colspan="1">1</td>
   </tr>
      <tr>
      <td style="text-align:center" align="center" colspan="1">b</td>
         <td style="text-align:center" align="center" colspan="1">1</td>
   </tr>
</tr>
</table>

The inference is done with posterior samples by running MCMC. 
