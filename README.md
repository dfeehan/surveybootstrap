
[![Travis-CI Build Status](https://travis-ci.org/dfeehan/surveybootstrap.svg?branch=master)](https://travis-ci.org/dfeehan/surveybootstrap)

Summary
================

The goal of the `surveybootstrap` package is to help people use the bootstrap
to estimate sampling uncertainty from surveys, including surveys with complex
sample designs.

For more information on the bootstrap for surveys, see
- Rao, JNK and Wu, CFJ (1988). Resampling inference with complex survey data. Journal of the American Statistical Association, 83(401):231–241.
- Rao, J., Wu, C., and Yue, K. (1992). Some recent work on resampling methods for complex surveys. Survey Methodology, 18(2):209–217.
- Rust, K. and Rao, J. (1996). [Variance estimation for complex surveys using replication techniques.](http://dx.doi.org/10.1177/096228029600500305) Statistical Methods in Medical Research, 5(3):283-310.
- Efron, B. and Tibshirani, R. J. (1993). An introduction to the bootstrap. Chapman and Hall/CRC.

The development of this software was supported by a grant from the National Institutes of Health (R01-HD075666).

Installing
-----------

You can install:

* the latest released version from CRAN with

    ```R
    install.packages("surveybootstrap")
    ```
            
* the latest development version from github with

    ```R
    install.packages("devtools")
    devtools::install_github("dfeehan/surveybootstrap")
    ```

            
Issues
---------
If you would like to suggest a feature or report a bug, please create an [issue](https://github.com/dfeehan/surveybootstrap/issues)

Citation
-----------

If you use our package for your research, please cite it so that we can continue to develop it.

- Feehan, Dennis M. and Salganik, Matthew J. (2016) "The surveybootstrap package." http://cran.r-project.org/package=surveybootstrap

Wishlist
--------

* update to use the `purrr` package, instead of `plyr` (but we'll probably need to wait until `purrr` supports parallelization)
