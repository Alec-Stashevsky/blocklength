# blocklength 0.1.3.9004

* A new vignette has been included on tuning and diagnosing problematic output from the selection functions

* the argument `sub_block_length = ` of `hhj()` has been changed to `sub_sample = ` to avoid confusion with the other tuning parameter, `pilot_block_length()`

* `rho.k.critical` has been changed to snake case `rho_k_critical` in the `$parameters` matrix from output of class 'pwsd' objects from `pwsd()`

* if `subsample = ` of `hhj()` is set directly, it is now rounded to the nearest whole number.

* `plot.pwsd()` now includes "darkmagenta" significance lines to match `pwsd()`

* `plot.pwsd()` now explicity includes an option to customize title with `main = `

# blocklength 0.1.3

* This is the first submission accepted to CRAN
* Minor changes to documentation and testing of examples


# blocklength 0.1.0

* This is the first submission of blocklength package.
