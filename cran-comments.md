## Test environments
* local R installation, R 4.1.2
* Ubuntu 20.04 (devel and release on GitHub Actions)
* macOS-latest (release on GitHub Actions)
* windows-latest (release on GitHub Actions)
* win-builder (devel and release)
* Fedora (devel on Rhub)
* Ubuntu (release on Rhub)
* Windows Server 2008 R2 (devel on Rhub)

## Resubmission - Version 0.1.5
This is resubmission upgrading to version 0.1.5 from 0.1.4.

Version Changes: 
  * Addressed CRAN reviewer comments for incorrect DOI link in `man/hhj.Rd`
  * Addressed CRAN review comments to change https://codecov.io/gh/Alec-Stashevsky/blocklength to https://app.codecov.io/gh/Alec-Stashevsky/blocklength

  
NOTE: Examples wrapped in \donttest are for hhj.R calls which take >5s to run


## R CMD check results

### Local Machine

No ERRORs, WARNINGs or NOTEs


### Win-builder

No ERRORs or WARNINGs or NOTEs


### R-hub

No ERRORs or WARNINGs.

One NOTEs:

  1. Package unavailable to check Rd xrefs: 'np'
     - There is a reference to package {np} in the .Rd file for pwsd, but there are no other dependencies on this package.


## Downstream dependencies

* None
