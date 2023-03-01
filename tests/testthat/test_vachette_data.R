test_that("vachette curves are accurate for iv", {

  set.seed(121)

  vd <- vachette_data(
    indivsam.obs = read.csv(system.file(package = "vachette", "examples", "iv-obs.csv")),
    output.typ = read.csv(system.file(package = "vachette", "examples", "iv-typ.csv")),
    vachette.covs = c("vachette.cov1" = 70),
    ref.dose = 1
  )

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "vachette-curves-iv.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_identical(vd$output.typ, reference)

})


test_that("vachette curves are accurate for oral-two-cov", {

  set.seed(121)

  vd <- vachette_data(
    indivsam.obs = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-v1-obs.csv")),
    output.typ = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-v1-typ.csv")),
    vachette.covs =  c("vachette.cov1" = 70, "vachette.cov2" = 30),
    ref.dose = 1
  )

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "vachette-curves-oral-two-cov-v1.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_identical(vd$output.typ, reference)

})
