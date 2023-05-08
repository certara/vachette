test_that("vachette curves are accurate for iv", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "iv-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "iv-typ.csv")),
    covariates = c("vachette.cov1" = 70),
    ref.dosenr = 1
  )

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "vachette-curves-iv.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_identical(vd$typ.data, reference)

})


test_that("vachette curves are accurate for oral-two-cov", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-typ.csv")),
    covariates =  c("vachette.cov1" = 70, "vachette.cov2" = 30),
    ref.dosenr = 1
  )

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "vachette-curves-oral-two-cov-v1.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_identical(vd$typ.data, reference)

})
