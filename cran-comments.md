## Test environments
* local R installation, R 4.0.3
* ubuntu 20.04 (devel and release on GitHub Actions), R 4.0.3
* macOS-latest (GitHub Actions)
* win-builder (devel and release)


## R CMD check results


### Local Machine

There were no ERRORs, WARNINGs or NOTEs

* This is a new release.
* This is a first time submission.


### Win-builder

I did get two NOTEs about the length of time to run the example code for hhj.R
I have reduced the time as much as possible to compute the example for `hhj()`,
but the function is somewhat intensive inherently. 

I also got a WARNING about building documentation, pasted below:

```
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! Package inputenc Error: Unicode char ‐ (U+2010)
(inputenc)                not set up for use with LaTeX.
```

### R-hub

No ERRORs or WARNINGs, 3 NOTEs

  1. This is a new submisions
  2. This run time on example for `hhj()` is >5s
  3. Package unavailable to check Rd xrefs: 'np'
     - There is a reference to package {np} in the .Rd file for pwsd, but there are no other dependencies on this package.

## Downstream dependencies

* There are no downstream dependencies
