# Subject-specific Dirichlet-Multinomial regression
Building upon the Dirichlet-multinomial regression framework, we develop an high-dimensional Bayesian hierarchical model for microbiota analysis that exploits subject-specific regression coefficients to simultaneously borrow strength across districts and include complex interactions between diet and clinical factors if supported by the data. 

For details regarding posterior computation see references section.

It is implemented in C++ through the use of [Rcpp](http://www.rcpp.org/) and [RcppArmadillo](https://dirk.eddelbuettel.com/code/rcpp.armadillo.html).

Authors: Matteo Pedone, Francesco Stingo.

Maintainer: Matteo Pedone.

## References

TBA

## Example

For a quick demonstration of the use of the function please check the script file in `demo/` directory.

## Installation

You can install the package from this GitHub repository, using the **devtools** package, running the following command:

```
devtools::install("mattpedone/subject-specific-dm-reg")
```
