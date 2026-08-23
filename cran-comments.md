## Test environments

* Local: Ubuntu 26.04, R 4.6.1 -- `R CMD check --as-cran`
* (Please also run win-builder / R-hub before submission, per the checklist.)

## R CMD check results

0 errors | 0 warnings | 0 notes, on a machine with `mrfDepth` installed and
a full LaTeX/HTML-tidy toolchain.

Notes seen only in restricted check environments, and why they do not apply
on CRAN:

* GitHub URL "Not Found" notes disappear once the public repository exists.
* "checking PDF version of manual" needs the `inconsolata` LaTeX package;
  the manual builds cleanly where a full TeX installation is present.
* HTML-validation / math-rendering notes require `tidy` and `V8`, which are
  present on CRAN's check machines.

## New submission

This is the first submission of raiseR.

The robust methods depend on the `mrfDepth` package (declared in Imports)
for the projection outlyingness measure used in the published methodology.
Examples and tests that require it are guarded with `\donttest{}` and
`testthat::skip_if_not_installed("mrfDepth")` respectively, so a check
environment without it still passes; with it installed they run in full.
