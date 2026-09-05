## Resubmission

This is a resubmission addressing the points raised by Konstanze Lauseker.

Changes made:

1. All acronyms are now explained on first use in the Description text:
   Simultaneous Raise Regression (SRR), Sequential Variance Inflation
   Factor (SVIF), Variance Inflation Factor (VIF), Robust Variance
   Inflation Factor (RVIF), Median Absolute Deviation Normalized (MADN).
   MM-estimates is the name coined by Yohai (1987) for the estimator and
   is cited accordingly rather than expanded.

2. The Title has been shortened to under 65 characters:
   "Raise Regression and Robust Methods for Multicollinearity" (57 chars).

3. \value tags have been added to all exported methods that were missing
   them (the six plot methods and the raiseReg-diagnostics page),
   describing the class and meaning of each return value.

4. All \dontrun{} wrappers in examples have been replaced with \donttest{}.
   These examples require the 'mrfDepth' package (declared in Imports) and
   run correctly when it is installed.

5. The package no longer writes to .GlobalEnv. The seed handling in the
   internal projection-outlyingness function now uses withr::with_seed(),
   which saves and restores the RNG state without modifying the global
   environment. 'withr' has been added to Imports.

## Test environments

* Local: Ubuntu 26.04, R 4.6.1 -- 0 errors, 0 warnings, 0 notes
* Local: Windows 11, R 4.6.1 -- 0 errors, 0 warnings, 0 notes

## R CMD check results

0 errors | 0 warnings | 0 notes (with mrfDepth installed and a full
LaTeX/HTML-tidy toolchain).

The only remaining note in restricted environments is the standard
"New submission" note.
