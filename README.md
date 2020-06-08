
<!-- README.md is generated from README.Rmd. Please edit that file -->
GECKO
=====

Accurate genetic and environmental covariance estimation with composite likelihood in genome-wide association studies

Installation
------------

To install the GECKO, it's easiest to use the 'devtools' package.

``` r
#install.packages("devtools")
library(devtools)
install_github("borangao/GECKO")
```

Usage
-----

The following help page will provide quick references for GECKO package and the example command lines:

``` r
library(GECKO)
package?GECKO
###Read in the summary statistics
Z1<-fread(summary_name_1)
Z2<-fread(summary_name_2)
###Need to specify the number of individuals within each study: n1,n2, and number of the overlapping individuals in the two studies:ns
###if the two samples are from the separate studies, nsin = 0, and Fix_Vein = 1
###if the two samples are from the same study, nsin need to be specified
n1in<-round(mean(Z1$N))
n2in<-round(mean(Z2$N))

#nsin number of overlapped samples
LDscorein<-read.table("/eur_chr_all.l2.ldscore.gz",header = T)
Weightin = T
Fix_Vein = T #(if two studies have non-overlapped samples, otherwise false)
Test = T
###Use GECKO function for the GECKO program, Z1, Z2 are the summary statistic files 
###n1, n2, ns are the number of individuals from the two studies, and the number of overlapped individuals
### LDscore is the LDscore output from the LD score software. We suggest to use 1000 genome reference panel for the LD score calculation
### Weightin = T is always set for GECKO to improve the efficiency
### Test = T is you want to test if the two trait are significantly correlated, otherwise you could specify it to be false
Result<-GECKO(Z1,Z2,n1in,n2in,nsin,LDscore,Weightin,Fix_Vein,Test)
```

Development
-----------

This package is developed and maintained by Boran Gao (<borang@umich.edu>).
