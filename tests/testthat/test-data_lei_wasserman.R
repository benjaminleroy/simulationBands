test_that("inner_discrete_mass_cde tests (basic)", {
  cde_vec <- c(.11,.19,.09,.1,.3,.21)
  delta_y <- 1
  mass_expected <- c(.3,.49,.09,.19, 1,.7)
  mass <- inner_discrete_mass_cde(cde_vec, delta_y)
  testthat::expect_equal(mass, mass_expected)

  reorder <- sample(6)
  cde_vec2 <- cde_vec[reorder]
  mass_expected2 <-  mass_expected[reorder]
  mass2 <- inner_discrete_mass_cde(cde_vec2, delta_y)
  testthat::expect_equal(mass2, mass_expected2)
})
