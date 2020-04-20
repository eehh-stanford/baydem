library(testthat)
library(baydem)

test_check("baydem")
list.files("analysis",
  full.names = TRUE,
  pattern = "FINAL"
)
