## Test environments
* local R installation, R 4.0.3
* ubuntu 20.04 (devel and release on GitHub Actions), R 4.0.3
* macOS-latest (GitHub Actions)
* win-builder (devel and release)


## R CMD check results

* This is a new release.


### Local Machine

There were no ERRORs, WARNINGs or NOTEs


### Win-builder

1 ERROR and 1 WARNING due to failed pdf latex compilation. This is a well known issue and I believe it is due to the pdflatex engine used on the win-builder server, because there is unicode (hyphens only) contained in .Rd files. It would be nice to be able to keep the hyphens in there for readability but I can remove if needed. Output included below:

```
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! Package inputenc Error: Unicode char ‐ (U+2010)
(inputenc)                not set up for use with LaTeX.
```

### R-hub

No ERRORs or WARNINGs.

Two NOTEs:

  1. This is a new submission
  2. Package unavailable to check Rd xrefs: 'np'
     - There is a reference to package {np} in the .Rd file for pwsd, but there are no other dependencies on this package.


## Downstream dependencies

* None
