
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DANTsR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jeffduda/DANTsR.svg?branch=master)](https://travis-ci.org/jeffduda/DANTsR)
[![muschellij2 Travis build
status](https://travis-ci.com/muschellij2/DANTsR.svg?branch=master)](https://travis-ci.com/muschellij2/DANTsR)

<!-- badges: end -->

DANTsR provides diffusion imaging extensions to ANTsR.

“Tensors? Richard had no idea what a tensor was, but he had noticed that
when math geeks started throwing the word around, it meant that they
were headed in the general direction of actually getting something
done.”

  - Neal Stephenson, Reamde (2011).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jeffduda/DANTsR")
```

## Example

``` r
library(DANTsR)
img = antsImageRead("yourfile.nii.gz")
fa = dtiAnisotrpyImage(img, "FractionalAnisotropy")
ra = dtiAnisotropyImage(img, "RelativeAnisotropy")
eig = dtiEigenSystem(img)
dec = dtiColorMap(img)
seeds = labelsToPoints( fa > 0.1 )
```
