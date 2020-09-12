# Scalable-Target-Marketing
Scalable Target Marketing: Distributed MCMC for Hierarchical Models

JMR_code.R: R application file that demonstrates the proposed algorithm

The following files must be compiled as part of the bayesm R package:
1.	Stage one of proposed algorithm
 -	rhierMnlRwMixtureParallel_rcpp.R
-	rhierMnlRwMixtureParallel_rcpp_loop.cpp

2.	Stage two of proposed algorithm
-	IndepMetropParallel_rcpp.R
-	IndepMetropParallel_rcpp_loop.cpp

3.	Utility file
-	llmnlParallel_rcpp.cpp
