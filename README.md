# Robust Parameter Design Based on Time-Consuming Bi-fidelity Simulations with Autoregressive Shifted Log Loss Gaussian Process Model

This repository contains the source code for the paper “Robust Parameter Design Based on Time-Consuming Bi-fidelity Simulations with Autoregressive Shifted Log Loss Gaussian Process Model.”

Robust parameter design is a quality improvement method to mitigate the effect of input noise by minimizing the expected quality loss (EQL). To reduce simulation time in RPD based on bi-fidelity simulations, the autoregressive shifted log loss Gaussian process (ASLLGP) model is proposed to estimate the high-fidelity (HF) EQL using HF and low-fidelity (LF) simulation data. The excellent performance of the proposed model is illustrated with two examples. Descriptions on how to run the codes to reproduce the results for the two examples are presented as follows. 

**Requirements: Matlab R2020a, Global Optimization Toolbox, Parallel Computing Toolbox.**

**Instructions for reproducing the results in Example 1:**

Please open and run the `ShaftCompare3methods.m` script to get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models as given in Table 2. The code also gives the results in Table 3.

Run the `ShaftGridOptimizationEQL.m` script to construct Figure 3, i.e., which is a plot of the posterior mean, and lower and upper 95% credible limits for the EQL given by the ASLLGP, ALL and AGPL models fitted with data from one nested design versus $`x_c`$.

**Instructions for reproducing the results in Example 2:**

Please open and run the `PiezoCompare3methods.m` script to get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models as given in Table 5. The code also gives the results in Table 6.

Run the `PiezoContourPlotsofEQL.m` script to construct Figure 5, which is a contour plot of the posterior mean, and lower and upper 95% credible limits for the EQL given by the ASLLGP, ALL and AGPL models fitted with data from one nested design versus $`x_c=(x_1,x_2)`$.

**Data Files**

**1.	Example 1**

The `ShaftOutput.mat` file gives the maximum Von Mises stress outputs of the HF and LF simulators (stored in variable stress) and the mass of the small shaft output of both simulators (stored in variable mass) on the $`{11}^2`$ grid $`{\{0,0.1,…,1}\}^2`$.

The `Designs for the shaft example.mat` file gives the 100 maximin nested Latin hypercube designs and 100 maximin Latin hypercube designs for 3 single-fidelity models.

The `200 test points for the shaft example.mat` file gives the 200 test points, i.e., the grid {0,1/199,…,1}.

**2.	Example 2**

The `Designs for the piezoelectric actuator example.mat` file gives the 100 maximin nested Latin hypercube designs.

The `400 test points for the piezoelectric actuator example.mat` file gives the 400 test points from a Latin hypercube design.

The `Piezo True EQL on grid points.mat` file gives the true EQL at the points in $`{\{0,0.05,…,1}\}^2`$.

**Other Matlab scripts**

`AGPLexpectedloss.m`: Function for computing the posterior mean of the EQL and 95% upper and lower credible interval limits for the EQL given by the AGPL model.

`ALLexpectedloss.m`: Function for computing the posterior mean of the EQL and 95% upper and lower credible interval limits for the EQL given by the ALL model.

`ASLLexpectedloss.m`: Function for computing the posterior mean of the EQL and 95% upper and lower credible interval limits for the EQL given by the ASLLGP model.

`correlax.m`: Function for computing the prior correlations between points in the input domain.

`gpfit1level.m`: Function used to find the MLEs of the parameters of the GP $`G_l (∙)`$ in Section 3.1.

`gpfit2level.m`: Function used to find the MLEs of the parameters of the GP $`\Delta(∙)`$ in Section 3.1 and also the MLE of $`\rho`$.

`gpfitASLL1level.m`: Function used to find the estimate of the vector of parameters $`\mathbf{\Theta_l}`$ for the proposed model in Section 4.

`gpfitASLL2level.m`: Function used to find the estimate of the vector of parameters $`\mathbf{Θ_h}`$ for the proposed model in Section 4.

`gppredict.m`: Function used to compute point and interval predictions of the ASLLGP/ALL/ AGPL emulator.  

`invandlogdet.m`: Function used to compute the inverse and the log determinant of a positive definite matrix.

`lgwt.m`:	Please	download	this	file	from [https://ww2.mathworks.cn/matlabcentral/fileexchange/4540-legendre-gauss-quadrature- weights-and-nodes](url), and place it in the same folder as the other scripts.

`lossShaft.m`: Function for computing the loss function used in Example 1.

`lossPiezo.m`: Function for computing the loss function used in Example 2.

`NestedLHD.m`: Function for generating maximin nested Latin hypercube designs.

`PiezoDesignPoints.m`: Generates 100 designs and outputs at the design points for Example 2 (the results from one execution of this script are stored in `Designs for the piezoelectric actuator example.mat`, and the data in this .mat file are used to obtain the results presented in Example 2)

`PiezoelectricActuator.m`: Function for constructing the HF (fidelity=2) and LF (fidelity=1) simulators for Example 2.

`PiezoTestPoints.m`: Script used to generate 400 test points Example 2 (used to generate the points stored in `400 test points for the piezoelectric actuator example.mat`).

`PiezoTrueEQL.m`: Function for computing the true EQL for Example 2.

`PiezoTrueEQLonGridPoints.m`: Script used to compute the true EQL on the points in the grid $`{\{0,0.005,…,1}\}^2`$ for Example 2 (used to compute the data stored in `Piezo True EQL on grid points.mat`).

`ShaftDesignpoints.m`: Generates 200 designs and outputs at the design points for Example 1 (the results from one execution of this script are stored in `Designs for the shaft example.mat`, and the data in this .mat file are used to obtain the results presented in Example 1).

`ShaftInterp.m`: Interpolators for the mass of the small shaft output and the HF (fidelity=2) and LF (fidelity=1) maximum von Mises stress outputs on the grid $`{\{0,0.1,…,1}\}^2`$ for Example 1.

`ShaftTestpoints.m`: Script used to generate the 200 test points for Example 1 (used to generate the points stored in `200 test points for the shaft example.mat`).

`ShaftTrueEQL.m`: Function for computing the true EQL for Example 1.
