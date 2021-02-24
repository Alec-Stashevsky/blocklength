## Test environments
* local R installation, R 4.0.3
* ubuntu 20.04 (devel and release on GitHub Actions), R 4.0.3
* macOS-latest (GitHub Actions)
* win-builder (devel and release)

## Resubmission

CHANGES: 
  * Added references as doi hyperlinks in description section of DESCRIPTION
  * added \value to .Rd files for plot.hhj.Rd and plot.pwsd.Rd
  * changed all references from \dontrun to \donttest
    * NOTE: Examples wrapped in \donttest are for hhj.R calls which take >5s to run


### Previous Resubmissions
CHANGES: Updated .Rd files to remove unicode characters


## R CMD check results

* This is a new release.


### Local Machine

There were no ERRORs, WARNINGs or NOTEs


### Win-builder

No ERRORs or WARNINGs

One NOTE:
  1. This is a new submission

### R-hub

No ERRORs or WARNINGs.

Two NOTEs:

  1. This is a new submission
  2. Package unavailable to check Rd xrefs: 'np'
     - There is a reference to package {np} in the .Rd file for pwsd, but there are no other dependencies on this package.


## Downstream dependencies

* None
