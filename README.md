# Robust Parameter Design Based on Time-Consuming Bi-fidelity Simulations With Autoregressive Shifted Log Loss Gaussian Process Model

This repository contains the MATLAB codes for the paper titled “Robust Parameter Design Based on Time-Consuming Bi-fidelity Simulations With Autoregressive Shifted Log Loss Gaussian Process Model.”

Robust parameter design is a quality improvement method that mitigates the effect of input noise by minimizing the expected quality loss (EQL). To reduce simulation time in RPD based on bi-fidelity simulations, the autoregressive shifted log loss Gaussian process (ASLLGP) model is proposed to estimate the high-fidelity (HF) EQL using HF and low-fidelity (LF) simulation data. The excellent performance of the proposed model is illustrated with two examples. Descriptions on how the codes should be run to reproduce the results for the two examples are given below. 

**Requirements: MATLAB R2020a, Global Optimization Toolbox, Parallel Computing Toolbox.**

**Note: Please put all MATLAB scripts and data files in this repository in one folder on your computer.**

**Instructions for reproducing the results in Example 1:**

Please open and run the `BarCompare3methods.m` script to get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models that are presented in Table 2 and the paired sample t-statistics in Table 2. The code also gives the results in Table 3, which summarizes the EQL prediction performance for three single-fidelity GP models, and Table 4, which summarizes the robust optimization performance of the ASLLGP, ALL, and AGPL models.

Run the `BarGridOptimizationEQL.m` script to construct Figure 3, i.e., which is a plot of the posterior means, and upper and lower 95% credible interval limits for the EQL given by the ASLLGP, ALL and AGPL models fitted with the data from one nested design versus $x_c$. The true EQL is also plotted in the figure.

**Instructions for reproducing the results in Example 2:**

Please open and run the `PiezoCompare3methods.m` script to get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models that are presented in Table 6 and the paired sample t-statistics in Table 6. The code also gives the results in Table 7, which summarizes the robust optimization performance of the ASLLGP, ALL, and AGPL models.

Run the `PiezoContourPlotsofEQL.m` script to construct Figure 5, which is a contour plot of the posterior means of the EQL given by the ASLLGP, ALL and AGPL models fitted with the data from one nested design versus $\mathbf{x}_c=(x_1,x_2)$. The true EQL is also plotted in the figure.

**Data Files**

**1.	Example 1**

The `BarOutput.mat` file gives the maximum Von Mises stress outputs of the HF and LF simulators (stored in variable stress) and the mass of the small bar outputs of both simulators (stored in variable mass) on the $`{11}^2`$ grid $`{\{0,0.1,…,1}\}^2`$ of $`(x_c,x_e)`$ values.

The `Designs for the bar example.mat` file gives the 100 maximin nested Latin hypercube designs for fitting the ASLLGP, ALL, and AGPL models, the 100 maximin Latin hypercube designs for fitting the single-fidelity GP models, and the corresponding outputs at all design points.

The `200 test points for the bar example.mat` file  gives the 200 test points for Example 1, and the true EQLs at those points.

**2.	Example 2**

The `Designs for the piezoelectric actuator example.mat` file gives the 100 maximin nested Latin hypercube designs for fitting the ASLLGP, ALL, and AGPL models, and the corresponding outputs at all design points.

The `400 test points for the piezoelectric actuator example.mat` file gives the 400 test points from a Latin hypercube design for Example 2, and the true EQLs at those points.

The `Piezo True EQL on grid points.mat` file gives the true EQL at each $`\mathbf{x}_c\in{\{0,0.005,…,1}\}^2`$.

**MATLAB scripts other than the four main .m files mentioned in the instructions above**

`AGPLexpectedloss.m`: Function for computing the posterior mean of the EQL and 95% upper and lower credible interval limits for the EQL given by the AGPL model.

`ALLexpectedloss.m`: Function for computing the posterior mean of the EQL and 95% upper and lower credible interval limits for the EQL given by the ALL model.

`ASLLexpectedloss.m`: Function for computing the posterior mean of the EQL and 95% upper and lower credible interval limits for the EQL given by the ASLLGP model.

`correlax.m`: Function for computing the correlation function (5).

`gpfit1level.m`: Function used to find the MLEs of the parameters of the GP $`G_l (∙)`$ in Section 3.1.

`gpfit2level.m`: Function used to find the MLEs of the parameters of the GP $`\Delta(∙)`$ in Section 3.1 and also the MLE of $`\rho`$.

`gpfitASLL1level.m`: Function used to find the estimate of the vector of parameters $`\mathbf{\Theta}_l`$ for the proposed model in Section 4.

`gpfitASLL2level.m`: Function used to find the estimate of the vector of parameters $`\mathbf{\Theta}_h`$ for the proposed model in Section 4.

`gppredict.m`: Function used to compute the posterior mean functions and posterior covariance functions of the HF and LF responses for the autoregressive GP model.  
  
`invandlogdet.m`: Function used to compute the inverse and the log determinant of a positive definite matrix.

`lgwt.m`:	Please	download	this	file	from https://ww2.mathworks.cn/matlabcentral/fileexchange/4540-legendre-gauss-quadrature- weights-and-nodes, and place it in the same folder as the other scripts.

`lossBar.m`: Function for computing the loss function  $`L (∙)`$ used in Example 1.

`lossPiezo.m`: Function for computing the loss function $`L (∙)`$ used in Example 2.

`NestedLHD.m`: Function for generating maximin nested Latin hypercube designs.

`PiezoDesignPoints.m`: Generates 100 maximin nested Latin hypercube designs and computes the outputs at the design points for Example 2 (the results from one execution of this script are stored in the `Designs for the piezoelectric actuator example.mat` file, and the data in this .mat file are used to obtain the results presented in Example 2).

`PiezoelectricActuator.m`: Function that contains the HF (fidelity=2) and LF (fidelity=1) simulators for Example 2.

`PiezoTestPoints.m`: Script used to generate 400 test points for Example 2 and compute the true EQLs at those points (used to generate the data stored in the `400 test points for the piezoelectric actuator example.mat` file).

`PiezoTrueEQL.m`: Function for computing the true EQL for Example 2.

`PiezoTrueEQLonGridPoints.m`: Script used to compute the true EQL at each $`\mathbf{x}_c\in{\{0,0.005,…,1}\}^2`$ for Example 2 (used to compute the data stored in the `Piezo True EQL on grid points.mat` file).

`BarDesignpoints.m`: Generates 100 maximin nested Latin hypercube designs and 100 maximin Latin hypercube designs, and computes the outputs at all design points for Example 1 (the results from one execution of this script are stored in the `Designs for the bar example.mat` file, and the data in this .mat file are used to obtain the results presented in Example 1).

`BarInterp.m`:  Interpolators for the mass of the small bar output (which is the same for both HF and LF simulators) on the grid $`{\{0,0.1,…,1}\}`$ of $`x_c`$ values and the HF (fidelity=2) and LF (fidelity=1) maximum von Mises stress outputs on the grid $`{\{0,0.1,…,1}\}^2`$ of $`(x_c,x_e)`$ values for Example 1.

`BarTestpoints.m`: Script used to generate the 200 test points for Example 1 and compute the true EQLs at those points (used to generate the data stored in the `200 test points for the bar example.mat` file).

`BarTrueEQL.m`: Function for computing the true EQL for Example 1.
