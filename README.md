# InvGaussianQuad-R  

This the set of R codes used for the numerical examples of __"Inverse Gaussian quadrature and finite normal-mixture approximation of the generalized hyperbolic distribution"__ paper by Jaehyuk Choi([@jaehyukchoi](https://github.com/jaehyukchoi)), Yeda Du, Qingshuo Song([@songqsh](https://github.com/songqsh)).

## Paper Information
### Title 
Inverse Gaussian quadrature and finite normal-mixture approximation of the generalized hyperbolic distribution

### Abstract
In this study, a numerical quadrature for the generalized inverse Gaussian distribution is derived from the Gauss-Hermite quadrature by exploiting its relationship with the normal distribution. The proposed quadrature is not Gaussian, but it exactly integrates the polynomials of both positive and negative orders. Using the quadrature, the generalized hyperbolic distribution is efficiently approximated as a finite normal variance-mean mixture. Therefore, the expectations under the distribution, such as cumulative distribution function and European option price, are accurately computed as weighted sums of those under normal distributions. The generalized hyperbolic random variates are also sampled in a straightforward manner. The accuracy of the methods is illustrated with numerical examples.

### Links
[DOI](https://doi.org/10.1016/j.cam.2020.113302) | [arXiv](https://arxiv.org/abs/1810.01116) | [SSRN](http://ssrn.com/abstract=3259013)

## Files
* [InvGaussianQuad/igquad.R](InvGaussianQuad/igquad.R): the collection of functions sourced in the other R files
* [InvGaussianQuad/Fig1-IG-Moments.R](InvGaussianQuad/Fig1-IG-Moments.R): Figure 1
* [InvGaussianQuad/Fig2-GIG-MGF.R](InvGaussianQuad/Fig2-GIG-MGF.R): Figure 2
* [InvGaussianQuad/Fig3-Table1-GH-CDF-Sets.R](InvGaussianQuad/Fig3-Table1-GH-CDF-Sets.R): Table 1 (parameter sets) and Figure 3 
* [InvGaussianQuad/Fig4-GH-CDF-Param.R](InvGaussianQuad/Fig4-GH-CDF-Param.R): Figure 4
* [InvGaussianQuad/Table2-GH-CDF-Time.R](InvGaussianQuad/Table2-GH-CDF-Time.R): Table 2
* [InvGaussianQuad/Table3-GH-CDF-Extreme.R](InvGaussianQuad/Table3-GH-CDF-Extreme.R): Table 3
* [InvGaussianQuad/Table4-GH-RV.R](InvGaussianQuad/Table4-GH-RV.R): Table 4

## Citation
* __Choi, J., Du, Y., & Song, Q.__ (2021). Inverse Gaussian quadrature and finite normal-mixture approximation of the generalized hyperbolic distribution. *Journal of Computational and Applied Mathematics*, 388, 113302. https://doi.org/10.1016/j.cam.2020.113302
