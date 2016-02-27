# CAMer

**CAMer** (Continuous Admixture Modeler) is an R package for Continuous Admixture Modeling (CAM), based on the result of **MALDmef** (available on http://www.picb.ac.cn/PGG/resource.php).

##Installation

One way is to use `devtools`

```r
library(devtools)
install_github("david940408/CAMer")
```

Another way is to download the repo from the website of population genetic group: http://www.picb.ac.cn/PGG/resource.php or https://github.com/david940408/CAMer, and run `R CMD INSTALL path/to/CAMer` in the shell/command line.

##Notes for Windows Users

Windows users need to install a suitable version of Rtools, which can be downloaded from [CRAN](https://cran.r-project.org/). They might also need to change system path when installing Rtools and add paths to the binaries of the installed R into the system variable named PATH. See the requirement of `Rcpp` and `RcppArmadillo` packages for more information.

##Vignette

See [An Introduction to CAMer package](https://github.com/david940408/CAMer/blob/master/inst/doc/intro.md), the .md file under inst/doc/.
