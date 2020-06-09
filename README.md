
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
Z1<-data(sumstat_1)
Z2<-data(sumstat_2)
####summary_stat_study_1 and summary_stat_study_1 are the file names of the summary statistics with column name chr, bp, SNP, A1, A2, N, Z, P representing chromosome, base pair position, SNP iD, major allele, minor allele, number of individuals in the study, Z score, P value
###Need to specify the number of individuals within each study: n1,n2, and number of the overlapping individuals in the two studies:ns
###if the two samples are from the separate studies, nsin = 0, and Fix_Vein = 1
###if the two samples are from the same study, nsin need to be specified
n1in<-round(mean(sumstat_1$N))
n2in<-round(mean(sumstat_2$N))


#nsin number of overlapped samples

Weightin = T
Fix_Vein = T #(if two studies have non-overlapped samples, otherwise false)
Test = T
###Use GECKO function for the GECKO program, Z1, Z2 are the summary statistic files 
###n1, n2, ns are the number of individuals from the two studies, and the number of overlapped individuals
### LDscore is the LDscore output from the LD score software. We suggest to use 1000 genome reference panel for the LD score calculation
### Weightin = T is always set for GECKO to improve the efficiency
### Test = T is you want to test if the two trait are significantly correlated, otherwise you could specify it to be false
Result<-GECKO_R(sumstat_1,sumstat_2,n1in,n2in,0,ldscore,Weightin,T,Test)
```

Development
-----------

This package is developed and maintained by Boran Gao (<borang@umich.edu>).
