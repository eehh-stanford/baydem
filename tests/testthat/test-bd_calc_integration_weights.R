test_that("test bd_calc_integration_widths", {
  expect_equal(bd_calc_integration_weights(c(-1.5,2,3,4,7)), c(1.75,2.25,1,2,1.5))
})

