# SparsePortfolioSelection

## Description
SparsePortfolioSelection is a Julia package which computes sparse L2-regularized efficient Markowitz portfolios. Sparsity is enforced explicitly via a cardinality constraint (or by penalizing cardinality, if you edit one line of the code), which is implemented by decomposing the problem as a mixed-integer linear programming and imposing cuts via outer-approximation. The algorithm is warm-started using a quasi-ADMM method with random-restarts and subsequently solves quadratic programs at each iteration to generate cuts.

## Development Notes
This package should be used at your own risk.
If you run into any issues, please file them [here](https://github.com/ryancorywright/SparsePortfolioSelection.jl/issues/new).

## A High-Level Overview
I gave a high-level talk on this package at INFORMS 2018. The slides are available [here](https://ryancorywright.github.io/pdf/A_scalable_algorithm_talk.pdf)

## Reproducing output from the paper
You can reproduce Tables 3-9 by running the script "createTablexResults.jl", where x is the relevant table number. Note that you will need to manually change the cardinality constraint right-hand-side value (i.e. run each script seperately for k=5, 10, 20), and change the file path to point to the right file on your system. Note that all experiments were run using Julia 1.2, JuMP.jl 0.18.5, CPLEX 12.8.0 and Mosek version 9.0. More up to-date versions of Julia, CPLEX or Mosek should be fine, but the code will not work with JuMP.jl 0.19 or higher.

## Citing SparsePortfolioSelection.jl

If you use SparsePortfolioSelection.jl, we ask that you please cite the following [paper](https://ryancorywright.github.io/pdf/Bertsimas_CoryWright_SparsePortfolios_2018.pdf):
```
@article{BCw_sparse2019,
	title = {A scalable algorithm for sparse portfolio selection},
	journal = {	arXiv preprint arXiv:1811.00138},
	author = {Bertsimas, Dimitris and Cory-Wright, Ryan},
	year = {2018}
}
```

Another paper which could potentially be useful to look at is the following [paper](https://ryancorywright.github.io/pdf/UnifiedFrameworkforMIO.pdf)


[![Build Status](https://travis-ci.org/ryancorywright/SparsePortfolioSelection.jl.svg?branch=master)](https://travis-ci.org/ryancorywright/SparsePortfolioSelection.jl)

## Related packages:

My colleague Jean Pauphilet has created a package which computes sparse L2-regularized efficient regressors. It is available [here](https://github.com/jeanpauphilet/SubsetSelectionCIO.jl). 