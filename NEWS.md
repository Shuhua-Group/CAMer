##v0.3

###Implementation
* Re-implemented part of `singleCAM`, `CAM` and `HI` in C++ to achieve faster speed. Hence this package now depends on `Rcpp` package (and `RcppArmadillo` for linear algebra as well).
* Removed compatibility with older versions of R. Now only support R (>= 3.0.0).
* Removed dependency of `snow` package (since newer versions of R has `parallel` package installed by default).
* Included notes for Windows users in README.md.
* Updated documentations accordingly, i.e. remove the parts explaining `parallel` and `snow` packages.

###Miscellany
* Reduced the size of data.
* Corrected a typo in README.md (comand->command).




##v0.2

###Statistics
* Used psedusovalues and signed-rank test to select best model(s).
* Used medians to determine the best-fit model.

###Miscellany
* Forced the levels of Model in data frames to be in the order in manuscript.
* Set starting time of HI to be the same as its ending time.
