An R Package for Causal Effect Estimation in Graphical Models with
Unmeasured Variables
================

- [Installation](#installation)
- [Quick Start](#quick-start)

This package is built for estimating the Average Causal Effect (ACE) in
graphical models with unmeasured variables. This package is an
implementation of the proposed estimators by , based on the theory of
influence functions and targeted minimum loss based estimation (TMLE).

If you find this package useful, please cite:

``` r
print("placeholder")
```

Graphical models with unmeasured variables can be depicted via the
Acyclic Directed Mixed Graphs (ADMG). For example, consider the
following ADMG, where $A$ is the treatment variable and $Y$ is the
outcome variable: ![](ADMG.png)

## Installation

To install, run the following code in terminal:

``` bash
# install the devtools package first if it's not yet installed
devtools::install_github("annaguo-bios/ADMGtmle")
```

The source code for `ADMGtmle` package is available on GitHub at
[ADMEtmle](https://github.com/annaguo-bios/ADMGtmle/tree/main).

## Quick Start

The main function in this package is `ADMGtmle()`, which estimates the
Average Causal Effect (ACE) using both TMLE esetimator and onestep
estimator in graphical models with unmeasured variables. To get a flavor
of how to use this package, we provide a quick example below:

``` r
library(ADMGtmle) # load the package
head(data_fig_4a) # take a glance of the data, which is a simulated dataset under above figure (a).
```

    ##            X          U A        M1         M2          L         Y
    ## 1 0.98890930  3.4036578 0 0.5423635 -1.6361458  0.7632388  3.257711
    ## 2 0.39774545 -0.7055857 0 2.4330827 -0.6538274  3.9498004  5.338658
    ## 3 0.11569778  1.0320105 1 4.5009622  2.0672577 11.4239744 19.819490
    ## 4 0.06974868  1.6994103 1 2.8610542 -1.1488686  5.0942460 10.803918
    ## 5 0.24374939  2.9995114 1 2.6677580 -1.1273291  3.2955105  8.205794
    ## 6 0.79201043  2.9700555 1 3.3190450  3.5674934 10.1107529 21.442111

``` r
est <- ADMGtmle(a=c(1,0),data=data_fig_4a, vertices=c('A','M','L','Y','X'),
                bi_edges=list(c('A','Y')),
                di_edges=list(c('X','A'), c('X','M'), c('X','L'),c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')),
                treatment='A', outcome='Y',
                multivariate.variables = list(M=c('M1','M2')))
```

    ## [1] "The treatment is primal fixable in the provided graph."
    ## [1] "The graph is nonparametrically saturated."
    ## TMLE estimated ACE: 1.83; 95% CI: (0.7, 2.97) 
    ## Onestep estimated ACE: 1.81; 95% CI: (0.65, 2.96)

The code above estimates the ACE of treatment $A$ on outcome $Y$,
defined as $E(Y^1)-E(Y^0)$, using the data `data_fig_4a` generated based
on Figure (a). The function `ADMGtmle()` takes the following arguments:

- `a`: a vector of length 2, specifying the values of treatment $A$ to
  compare. For example, `a=c(1,0)` compares causal effect under $A=1$
  verse $A=0$.
- `data`: a data frame containing the data.
- `vertices`: a vector of strings, specifying the names of all vertices
  in the ADMG.
- `bi_edges`: a list of bi-directed edges in the ADMG.
- `di_edges`: a list of directed edges in the ADMG. For example,
  `c('X','A')` specifies an edge from $X$ to $A$.
- `treatment`: a string, specifying the name of the treatment variable.
- `outcome`: a string, specifying the name of the outcome variable.
- `multivariate.variables`: a list of variables that are multivariate.
  For example, `list(M=c('M1','M2'))` specifies that $M$ is a
  multivariate variable with components $M1$ and $M2$.
