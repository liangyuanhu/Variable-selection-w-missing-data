# Variable-selection-w-missing-data
A general variable selection approach in the presence of missing data in both covariates and outcomes. 
This approach exploits the flexibility of machine learning modeling techniques and bootstrap imputation, 
which is amenable to nonparametric methods in which the effect sizes of predictor variables are not naturally
defined as in parametric models. Six methods are considered and compared: XGBoost, Random Forests, 
Conditional Random Forests, Bayesian Additive Regression Trees, lasso, stepwise backward selection. 
Details can be found in Hu et al. (2021). 
Variable selection with missing data in both covariates and outcomes: Imputation and machine learning. Statistical Methods in Medical Research 30 (12), 2651-2671. 
