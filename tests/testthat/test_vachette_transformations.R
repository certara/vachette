test_that("vachette transformations are accurate for iv", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "iv-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "iv-typ.csv")),
    covariates = c(WT = 70),
    mappings = c(OBS = "DV",
                 x = "time"),
    ref.dosenr = 1
  )

  vd <- apply_transformations(vd)

  #write.csv(vd$obs.all, "inst/test_files/iv_obs_all.csv", row.names = FALSE)

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "iv_obs_all.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$obs.all, tidyr::as_tibble(reference))

})

test_that("vachette transformations are accurate for oral absorption", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "oral-absorption-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "oral-absorption-typ.csv")),
    covariates = c(WT = 70),
    mappings = c(OBS = "DV",
                 x = "time"),
    ref.dosenr = 1
  )

  vd <- apply_transformations(vd)

  #write.csv(vd$obs.all, "inst/test_files/oral_absorption_obs_all.csv", row.names = FALSE)

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "oral_absorption_obs_all.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$obs.all, tidyr::as_tibble(reference))

})

test_that("vachette transformations are accurate for oral absorption", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "oral-two-cov-typ.csv")),
    covariates = c(WT = 70, AGE = 30),
    mappings = c(OBS = "DV",
                 x = "time"),
    ref.dosenr = 1
  )

  vd <- apply_transformations(vd)

  #write.csv(vd$obs.all, "inst/test_files/oral_two_cov_obs_all.csv", row.names = FALSE)

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "oral_two_cov_obs_all.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$obs.all, tidyr::as_tibble(reference))

})


test_that("vachette transformations are accurate for sigmoid", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "sigmoid-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "sigmoid-typ.csv")),
    covariates = c(WT = 70),
    mappings = c(OBS = "DV",
                 x = "bmx"),
    ref.dosenr = 1,
    log.x = TRUE
  )

  vd <- vd |>
    apply_transformations(asymptote_right = TRUE,
                          asymptote_left = TRUE,
                          zero_asymptote_right = FALSE,
                          zero_asymptote_left = TRUE)

  #write.csv(vd$obs.all, "inst/test_files/sigmoid_obs_all.csv", row.names = FALSE)

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "sigmoid_obs_all.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  #look into why vd$obs.all is not a tbl
  #testthat::expect_equal(vd$obs.all, tidyr::as_tibble(reference))
  testthat::expect_equal(vd$obs.all, reference)


})


test_that("vachette transformations are accurate for indirect response", {

  set.seed(121)

  vd <- vachette_data(
    obs.data = read.csv(system.file(package = "vachette", "examples", "indirect-response-obs.csv")),
    typ.data = read.csv(system.file(package = "vachette", "examples", "indirect-response-typ.csv")),
    covariates = c(WT = 30),
    mappings = c(OBS = "DV",
                 x = "time"),
    ref.dosenr = 1
  ) |>
    apply_transformations(asymptote_right = TRUE,
                          asymptote_left = TRUE,
                          zero_asymptote_right = FALSE,
                          zero_asymptote_left = TRUE)

  #write.csv(vd$obs.all, "inst/test_files/indirect_response_obs_all.csv", row.names = FALSE)

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "indirect_response_obs_all.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$obs.all, tidyr::as_tibble(reference))

})
