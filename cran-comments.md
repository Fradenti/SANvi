Please note this is the first submission of this package to CRAN.

## R CMD check results

- Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces no errors, no warnings, nor notes.

- Running `devtools::check(remote = TRUE, manual = TRUE)` produces

0 errors | 0 warnings | 1 note

> checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Francesco Denti <francescodenti.personal@gmail.com>’
  
  New submission
  
- Running `devtools::check_win_devel()` produces

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Francesco Denti <francescodenti.personal@gmail.com>'
  
  New submission

- Finally, this package, in its current state, also passes all the standard checks performed via *GitHub actions*.

## Downstream dependencies

There are currently no downstream dependencies for this package.
