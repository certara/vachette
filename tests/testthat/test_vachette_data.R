test_that("vachette curves are accurate for iv", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "iv-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "iv-typ.csv")),
    covariates = c(WT = 70),
    mappings = c(OBS = "DV",
                 x = "time"),
    ref.dosenr = 1
  )

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "obs-curves-iv.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$obs.orig, tidyr::as_tibble(reference))

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "typ-curves-iv.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$typ.orig, tidyr::as_tibble(reference))

})


test_that("vachette curves are accurate for oral-two-cov", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-typ.csv")),
    covariates =  c(WT = 70, AGE = 30),
    mappings = c(OBS = "DV",
                 x = "time"),
    ref.dosenr = 1
  )

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "obs-curves-oral-two-cov.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$obs.orig, tidyr::as_tibble(reference))

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "typ-curves-oral-two-cov.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$typ.orig, tidyr::as_tibble(reference))

})
