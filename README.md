
<!-- README.md is generated from README.Rmd. Please edit that file -->

# poset

The **poset** package provides simple and efficient statistical routines
for partially ordered data, such as multivariate ordinal response under
consensus or prioritized order. The current version focuses on the win
ratio/net benefit approach (Mao 2024) via generalized pairwise
comparisons (Buyse 2010).

## Installation

Install **poset** from CRAN with:

``` r
install.packages("poset")
```

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("lmaowisc/poset")
```

## Examples

Here is a basic example for two-sample testing and regression.

``` r
library(poset)
## data example
head(liver)
#>   R1NASH R2NASH Sex    AF Steatosis  SSF2  LSN
#> 1      3      2   M FALSE        30  0.21 2.33
#> 2      1      1   F FALSE         5  0.38 2.86
#> 3      4      2   M FALSE        70  0.58 3.65
#> 4      4      4   F  TRUE        30 -0.08 2.73
#> 5      4      3   M  TRUE        70 -0.04 2.53
#> 6      3      3   M FALSE        10  0.02 2.88
```

### Compare bivariate ratings by fibrosis stage

``` r
Y1 <- liver[liver$AF, c("R1NASH", "R2NASH")] # advanced
Y0 <- liver[!liver$AF, c("R1NASH", "R2NASH")] # not advanced
wrtest(Y1, Y0)
#> Call:
#> wrtest(Y1 = Y1, Y0 = Y0)
#> 
#> Two-sample (Y1 vs Y0) win ratio/net benefit analysis
#> 
#> Number of pairs: N1 x N0 =  69 x 116  =  8004 
#>   Win: 4251 (53.1%)
#>   Loss: 2392 (29.9%)
#>   Tie: 1361 (17%)
#> 
#> Win ratio (95% CI): 1.78 (1.16, 2.73), p-value = 0.00856547
#> Net benefit (95% CI): 0.232 (0.065, 0.4), p-value = 0.006577537
```

### Regression analysis

``` r
Y <- 5 - liver[, c("R1NASH", "R2NASH")] # lower score is better
Z <- cbind("Female" = liver$Sex == "F",
           liver[, c("AF", "Steatosis",   "SSF2",  "LSN")]) # covariates
obj <- wreg(Y, Z) # fit model
obj
#> Call:
#> wreg(Y = Y, Z = Z)
#> 
#> n = 154 subjects with complete data
#> Comparable (win/loss) pairs: 9548/11781 = 81%
#> 
#>    Female         AF   Steatosis         SSF2        LSN
#>  -0.18956 -0.9660827 -0.02779146 -0.007926333 -0.1029914
```

``` r
summary(obj)
#> Call:
#> wreg(Y = Y, Z = Z)
#> 
#> n = 154 subjects with complete data
#> Comparable (win/loss) pairs: 9548/11781 = 81%
#> 
#> Newton-Raphson algoritm converged in 7 iterations
#> 
#>                coef exp(coef)  se(coef)      z Pr(>|z|)    
#> Female    -0.189560    0.8273  0.259988 -0.729 0.465934    
#> AF        -0.966083    0.3806  0.280313 -3.446 0.000568 ***
#> Steatosis -0.027791    0.9726  0.005281 -5.262 1.42e-07 ***
#> SSF2      -0.007926    0.9921  0.003953 -2.005 0.044953 *  
#> LSN       -0.102991    0.9021  0.125718 -0.819 0.412657    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>           exp(coef) exp(-coef) lower .95 upper .95
#> Female      0.82732    1.20872   0.49702    1.3771
#> AF          0.38057    2.62763   0.21970    0.6592
#> Steatosis   0.97259    1.02818   0.96258    0.9827
#> SSF2        0.99210    1.00796   0.98445    0.9998
#> LSN         0.90213    1.10848   0.70512    1.1542
#> 
#> Overall Wald test = 79.129 on 5 df,  p = 1.221245e-15
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-buyse2010" class="csl-entry">

Buyse, Marc. 2010. “Generalized Pairwise Comparisons of Prioritized
Outcomes in the Two-Sample Problem.” *Statistics in Medicine* 29 (30):
3245–57. <https://doi.org/10.1002/sim.3923>.

</div>

<div id="ref-mao2024" class="csl-entry">

Mao, Lu. 2024. “Win Ratio for Partially Ordered Data,” Under revision.

</div>

</div>
