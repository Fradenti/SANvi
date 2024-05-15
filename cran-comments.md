
We now report, for the fixed version of the package, the following

## R CMD check results

- Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces no errors, warnings, or notes.

- Running `devtools::check(remote = TRUE, manual = TRUE)` produces

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

Found the following (possibly) invalid DOIs:
    DOI: 10.1080/01621459.2021.1933499
      From: DESCRIPTION
      Status: Forbidden
      Message: 403

However, the links work when building the documentation.
  
- Running `devtools::check_win_devel()` produces

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Francesco Denti <francescodenti.personal@gmail.com>'
  
  New submission

- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.

## Downstream dependencies

There are currently no downstream dependencies for this package.
