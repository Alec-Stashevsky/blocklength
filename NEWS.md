# blocklength 0.1.5

* Minor adjustments to documentation per CRAN requests.


# blocklength 0.1.4


## Bug Fixes

* In `pwsd()`, `m_hat` was indexing the first insignificant lag of the correlation structure rather than the first *significant* lag prior to the consecutive run of insignificant lags.

* `pwsd()` and `plot.pwsd()` now output significance bands that match `rho_k_critical` exactly. Prior to this, the significance bands on the correlogram output were generated using the `plot.acf(ci = )` argument which led to misleading graphical representations of the implied hypothesis test's critical value in the PWSD method.

* In `hhj()`, `sub_block_length` has been changed to `sub_sample` to avoid confusion with the other tuning parameter, `pilot_block_length`


## Minor Changes

* A new vignette has been included on tuning and diagnosing problematic output from the selection functions!

* `pwsd()` now includes a new argument to override the implied hypothesis test by setting `m_hat` directly.

* In `pwsd()`, `rho.k.critical` has been changed to snake case `rho_k_critical` in the `$parameters` matrix from output of class 'pwsd' objects from `pwsd().`

* In `hhj()`, if `subsample = ` is set directly it is now rounded to the nearest whole number.

* `plot.pwsd()` now includes "darkmagenta" significance lines to match `pwsd().`

* `plot.pwsd()` now explicitly includes an option to customize title with `main = .`

* `hhj()` now includes a warning message if the supplied iteration limit `n_iter` is reached but still outputs an object of class 'hhj.'


# blocklength 0.1.3

* This is the first submission accepted to CRAN
* Minor changes to documentation and testing of examples


# blocklength 0.1.0

* This is the first submission of blocklength package.
