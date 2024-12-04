# NEWS for SBMTrees

## Version 1.1 (2024-12-04) - Initial Release

This release introduces the **SBMTrees** package, which implements a **Bayesian non-parametric framework** for imputing missing covariates and outcomes in longitudinal data under the **Missing at Random (MAR)** assumption. The core model, **Bayesian Trees Mixed-Effects Model (BMTrees)**, extends **Mixed-Effects BART** by using **centralized Dirichlet Process (CDP) Normal Mixture priors**, allowing it to handle non-normal random effects and errors, address model misspecification, and capture complex relationships in longitudinal data.

### Features:
- **BMTrees Prediction**: Predicts longitudinal outcomes based on mixed-effects models.
- **Sequential Imputation**: Imputes missing covariates and outcomes in longitudinal datasets with flexibility in specifying priors for random effects and errors.
- **Efficient Gibbs sampling**: Implemented in **C++** for improved performance.
- **Integration with mixedBART**: Includes **mixedBART** as a baseline model for additional flexibility in handling data.

### Key Functionalities:
- **BMTrees_prediction**: Predicts longitudinal outcomes based on mixed-effects models.
- **sequential_imputation**: Handles missing data imputation sequentially in longitudinal datasets.

### License:
This package is licensed under the **GNU General Public License version 2 (GPL-2)**. The core computations leverage code derived from the **BART3** package, originally developed by Rodney Sparapani, which is also under **GPL-2**.

### Known Issues:
- No significant issues at the moment. For more details on potential limitations, please refer to the documentation.
