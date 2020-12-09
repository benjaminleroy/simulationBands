test_that("cut_to_numeic basic tests (incomplete)", {
  my_string <- factor("(-Inf,1.092006e-05]") # scientific notation and infinity

  upper <- cut_to_numeric(my_string, .lower = FALSE)
  lower <- cut_to_numeric(my_string, .lower = TRUE)

  testthat::expect_equal(upper, 1.092006e-05)
  testthat::expect_equal(lower, -Inf)
})
