# Top-level test harness for `R CMD check` and `testthat::test_local()`.
# `test_check()` discovers `tests/testthat/test-*.R` files and runs them
# against the installed copy of the package (the one R CMD check sets up).

library(testthat)
library(germlinevaR)

test_check("germlinevaR")
