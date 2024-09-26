# Prediction inference of Non-linear Autoregressive Model (PINLAR)
We consider the prediction inference of the non-linear time series based on bootstrap. This repository provides the necessary code to reproduce the simulation and empirical studies in the paper **Bootstrap prediction inference of nonlinear autoregressive models** [Link](https://onlinelibrary.wiley.com/doi/full/10.1111/jtsa.12739).

## Usage of this repository
- The file **Simulation_example.R** is self-contained to perform prediction inference of a specific non-linear AR model in the format **X_{t+1} = a + log(b + abs(X_t))**. It will compute the bootstrap-based L_1 and L_2 point predictions; the QPI with fitted or predictive residuals, and the PPI centered at L_2 and L_1 optimal point predictions with fitted or predictive residuals.
- Other simulation studies can be performed similarly by applying some auxiliary functions from file **All_simulated_data.R**
- The file **Empirical_data_example.R** is another self-contained source to reproduce the empirical study with UnempRate data in the original paper. The algorithm is almost identical to the one in the file **Simulation_example.R**. The difference is that we do not need to simulate data and the model estimation method is changed.

## Note
I would like to recommend reading **Simulation_example.R** first and then **Empirical_data_example.R**. 

