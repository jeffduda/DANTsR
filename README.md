
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DANTsR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jeffduda/DANTsR.svg?branch=master)](https://travis-ci.org/jeffduda/DANTsR)
[![muschellij2 Travis build
status](https://travis-ci.com/muschellij2/DANTsR.svg?branch=master)](https://travis-ci.com/muschellij2/DANTsR)

<!-- badges: end -->

The goal of DANTsR is to provide diffusion imaging extensions to ANTsR.

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
