# Regression test: peak/valley landmark refinement must tolerate typical curves
# that yield multiple stationary points or a degenerate search window.
#
# Trigger: when two subjects share a covariate combination the typical
# simulation produces duplicate, overlapping curves (duplicate x-values). This
# collapses the peak-refinement window onto the whole curve, so `extremes()`
# returns several roots. Previously the unguarded `if (is.na(newpeak.refined))`
# check errored with "the condition has length > 1", leaving `obs.all` NULL.

test_that("apply_transformations completes when a covariate group has duplicate typical curves", {
  obs <- read.csv(system.file(package = "vachette", "test_files", "duplicate_typ_obs.csv"))
  typ <- read.csv(system.file(package = "vachette", "test_files", "duplicate_typ_typ.csv"))

  # wt = 50 / female is represented by two subjects -> duplicate typical curves
  dup <- typ[typ$wt == 50 & typ$sex == "female", ]
  testthat::expect_gt(dplyr::n_distinct(dup$ID), 1)

  vd <- suppressWarnings(suppressMessages(
    vachette_data(
      obs.data   = obs,
      typ.data   = typ,
      covariates = c(wt = 78, sex = "male"),
      model.name = "duplicate-typ"
    ) |>
      apply_transformations()
  ))

  # The transformation must complete (no length > 1 condition error)
  testthat::expect_false(is.null(vd$obs.all))
  # and retain all input subjects
  testthat::expect_setequal(unique(vd$obs.all$ID), unique(obs$ID))
})


test_that("extremes() reduces multiple minmax roots to the refined landmark closest to the detected peak", {
  # A search window that spans a full biexponential curve makes the degree-6
  # polynomial fit yield several stationary points. get.x.multi.landmarks must
  # still return a single, finite 'max' landmark rather than a length > 1 value.
  x <- c(0, seq(0.1, 132, length.out = 199))
  ka <- 1; ke <- 0.134 / 7.8; dose <- 100; v <- 7.8
  y <- (dose / v) * (ka / (ka - ke)) * (exp(-ke * x) - exp(-ka * x))
  y[1] <- 0

  # Duplicate the curve so each x appears twice (mimics duplicate typical curves)
  lm <- suppressWarnings(suppressMessages(
    vachette:::get.x.multi.landmarks(c(x, x), c(y, y))
  ))

  max_x <- lm$x[lm$type == "max"]
  testthat::expect_length(max_x, 1L)
  testthat::expect_true(all(is.finite(lm$x)))
})
