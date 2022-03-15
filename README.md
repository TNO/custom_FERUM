# Customized version of FERUM 4.1

To perfrom structural reliability analysis.

Work in progress, many changes are experimental, expect nothing ;)


## Changes

- additional distribution functions:
    - Generalized extreme value distribution (GEV)
    - 3-parameter lognormal distribution
    - sample based distribution - fitting kernel to the sample
    - non-parametric distribution - vector based definition
    - powered distribution type option, `F_n(x) = F_1(x)^n`

- new methods:
    - lambda method,  inspired by [doi.org/10.1016/j.strusafe.2015.09.004](doi.org/10.1016/j.strusafe.2015.09.004)
    - directional sampling with adaptive response surface (ADIS), inspired by [doi:10.1016/j.probengmech.2010.11.002](doi:10.1016/j.probengmech.2010.11.002)
    
- solvers:
    - `COBYLA` for `FORM` (relying on [NLopt](https://github.com/stevengj/nlopt))
- general:
    - small changes in the code: too many to list and were not properly recorded...
 
License: the original license applies.   
    
## Acknowledgements

The latest version of FERUM is available from https://www.sigma-clermont.fr/en/ferum

Originally developed at the University of California, Berkeley. Initiated by Terje Haukaas and Armen Der Kiureghian. http://projects.ce.berkeley.edu/ferum/

The scripts in this repo have been developed during my time at:
* Structural Reliability Department, Building and Construction Research,  [TNO](https://www.tno.nl/en/). **(2016 - 2022)**
* [Department of Structural Engineering](http://hsz.bme.hu/?language=en), Budapest University of Technology and Economics  . **(2014 - 2016)**