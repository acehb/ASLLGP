This repository contains MATLAB code for the paper titled "Robust Parameter Design Based on Time-Consuming Bi-fidelity Simulations With Autoregressive Shifted Log Loss Gaussian Process Model."

Robust parameter design is a quality improvement method that mitigates the effect of input noise by minimizing the expected quality loss (EQL). To reduce simulation time in RPD based on bi-fidelity simulations, the autoregressive shifted log loss Gaussian process (ASLLGP) model is proposed to estimate the high-fidelity (HF) EQL using HF and low-fidelity (LF) simulation data. The excellent performance of the proposed model is illustrated with three examples. Descriptions on how the codes should be run to reproduce the results for the three examples are given below. 
## Requirements: MATLAB R2020a, Global Optimization Toolbox, Parallel Computing Toolbox.

## Note: Please put all MATLAB scripts and data files in this repository in one folder on your computer.

## Instructions for reproducing the results in Example 1:

Please open and run the BarCompare3methods.m script to get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models that are presented in Table 2 and the paired sample t-statistics in Table 2. The code also gives the results in Table 3, which summarizes the EQL prediction performance for three single-fidelity GP models, and Table 4, which summarizes the robust optimization performance of the ASLLGP, ALL, and AGPL models.

### `BarGridOptimizationEQL.m`

Run this script to construct Figure 3. The figure plots the posterior means and upper and lower 95% credible interval limits for the EQL given by the ASLLGP, ALL, and AGPL models fitted using one nested design, together with the true EQL.

## Instructions for Reproducing Example 2

### `PiezoCompare3methods.m`

Run this script to obtain the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models presented in Table 6, together with the paired-sample t-statistics. The script also summarizes EQL prediction performance for three single-fidelity GP models and provides Table 7, which summarizes robust-optimization performance for the ASLLGP, ALL, and AGPL models.

### `PiezoContourPlotsofEQL.m`

Run this script to construct Figure 5. The figure is a contour plot of the true EQL and the posterior mean EQLs given by the ASLLGP, ALL, and AGPL models fitted using one nested design.

### `PiezoBoundedLossCompare3methods.m`

Run this script to obtain the summaries of the five EQL prediction performance measures and paired-sample t-statistics for the ASLLGP, ALL, and AGPL models using the bounded loss function. The script also summarizes EQL prediction performance for three single-fidelity GP models and robust-optimization performance under the bounded loss function.

## Data Files

### Example 1

- `BarOutput.mat`: Maximum von Mises stress outputs of the HF and LF simulators, stored in `stress`, and mass outputs of the small bar, stored in `mass`, on the interpolation grid.
- `Designs for the bar example.mat`: 500 maximin nested Latin hypercube designs for fitting the ASLLGP, ALL, and AGPL models; 500 maximin Latin hypercube designs for fitting the single-fidelity GP models; and the corresponding outputs at the design points.
- `200 test points for the bar example.mat`: The 200 test points for Example 1 and their true EQL values.

### Example 2

- `Designs for the piezoelectric actuator example.mat`: 500 maximin nested Latin hypercube designs for fitting the ASLLGP, ALL, and AGPL models, together with the corresponding outputs at the design points.
- `400 test points for the piezoelectric actuator example.mat`: The 400 Latin hypercube test points for Example 2 and their true EQL values under the unbounded and bounded loss functions.
- `Piezo True EQL on grid points.mat`: The true EQL values at all control-factor settings on the grid $x_c \in \{0, 0.005, \ldots, 1\}^2$, for both the unbounded and bounded loss functions.

## Supporting MATLAB Scripts

- `AGPLexpectedloss.m`: Computes the posterior mean EQL and 95% upper and lower credible interval limits for the AGPL model.
- `ALLexpectedloss.m`: Computes the posterior mean EQL and 95% upper and lower credible interval limits for the ALL model.
- `ASLLexpectedloss.m`: Computes the posterior mean EQL and 95% upper and lower credible interval limits for the ASLLGP model.
- `correlax.m`: Computes the correlation function.
- `gpfit1level.m`: Finds maximum likelihood estimates of the parameters of the one-level GP model.
- `gpfit2level.m`: Finds maximum likelihood estimates of the parameters of the two-level GP model and the autoregressive coefficient.
- `gpfitASLL1level.m`: Estimates the low-fidelity parameter vector for the proposed model.
- `gpfitASLL2level.m`: Estimates the high-fidelity parameter vector for the proposed model.
- `gpfitGPL.m`: Finds maximum likelihood estimates for the GP-for-loss model.
- `gpfitLL.m`: Finds maximum likelihood estimates for the lognormal loss model.
- `gpfitSLL.m`: Finds maximum likelihood estimates for the shifted log loss GP model.
- `gppredict.m`: Computes posterior mean and covariance functions for the HF and LF responses in the autoregressive GP model.
- `invandlogdet.m`: Computes the inverse and log determinant of a positive definite matrix.
- `lgwt.m`: Provides Gauss-Legendre quadrature weights and nodes. Download it from [MATLAB File Exchange](https://ww2.mathworks.cn/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes) and place it in the same folder as the other scripts.
- `lossBar.m`: Computes the loss function used in Example 1.
- `lossPiezo.m`: Computes the unbounded and bounded loss functions used in Example 2.
- `PiezoDesignPoints.m`: Generates 500 maximin nested Latin hypercube designs and 500 maximin Latin hypercube designs, and computes outputs at the design points for Example 2.
- `PiezoelectricActuator.m`: Contains the HF (`fidelity = 2`) and LF (`fidelity = 1`) simulators for Example 2.
- `PiezoTestPoints.m`: Generates 400 test points for Example 2 and computes their true EQL values under the unbounded and bounded loss functions.
- `PiezoTrueEQL.m`: Computes the true EQL under the unbounded or bounded loss function for Example 2.
- `PiezoTrueEQLonGridPoints.m`: Computes the true EQL at all control-factor settings on the grid for Example 2.
- `BarDesignpoints.m`: Generates 500 maximin nested Latin hypercube designs and 500 maximin Latin hypercube designs, and computes outputs at the design points for Example 1.
- `BarInterp.m`: Interpolates the bar mass output and the HF (`fidelity = 2`) and LF (`fidelity = 1`) maximum von Mises stress outputs for Example 1.
- `BarTestpoints.m`: Generates 200 test points for Example 1 and computes their true EQL values.
- `BarTrueEQL.m`: Computes the true EQL for Example 1.
