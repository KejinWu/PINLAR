# PINLAR
The prediction inference of non-linear time series based on bootstrap. This repository provides the necessary to reproduce the simulation and empirical studies of the paper **Bootstrap prediction inference of nonlinear autoregressive models** [Link](https://onlinelibrary.wiley.com/doi/full/10.1111/jtsa.12739)

## Usage of this repository
The file *Simulation_example.R* is self-contained to perform prediction inference of a specific non-linear AR model in the format *X_{t+1} = a + log(b + abs(X_t))*. It will compute the bootstrap-based L_1 and L_2 point predictions; the QPI with fitted or predictive residuals, and the PPI centered at L_2 and L_1 optimal point with fitted or predictive residuals predictions. 


