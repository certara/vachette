test_that("covariates are correctly extracted if missing reference values", {

  #DF
  indivsam.obs <- read.csv(system.file(package = "vachette", "examples", "oral-two-cov-obs.csv"))
  vachette.covs <- vachette:::.process_covariates(vachette.covs = c("WT", "AGE"), indivsam.obs)
  testthat::expect_equal(vachette.covs, c("WT" = "70", "AGE" = "30"))

  #Test Tibble
  mtcars$cyl <- paste0("cyl_", mtcars$cyl)
  vachette.covs <- vachette:::.process_covariates(vachette.covs = c("mpg", "cyl"), mtcars)
  testthat::expect_equal(vachette.covs, c("mpg" = "19.2", "cyl" = "cyl_8"))

  vachette.covs <- vachette:::.process_covariates(vachette.covs = c("mpg" = 18.7, "cyl" = "cyl_6"), mtcars)
  testthat::expect_equal(vachette.covs, c("mpg" = "18.7", "cyl" = "cyl_6"))

  vachette.covs <- vachette:::.process_covariates(vachette.covs = c("cyl", "hp" = 123), mtcars)
  testthat::expect_equal(vachette.covs, c("cyl" = "cyl_8", "hp" = "123"))
})

