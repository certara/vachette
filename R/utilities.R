vachette_env <- new.env()

vachette_env$summary_out <- list(errors = list(),
                    warnings = list(),
                    values = list()
)

vachette_env$warning_out <- NULL

log_output <- function(message) {

  if (exists("vachette_log_file", envir = vachette_env)){
    log_file <- get("vachette_log_file", envir = vachette_env)
  } else {
    log_file <- NULL
  }

  if (!is.null(log_file)) {
    cat(message, file = log_file, append = TRUE, sep = "\n")
  } else {
    cat(message, sep = "\n")
  }
}



#' Summarise \code{vachette_data}
#'
#' Summary generic used to return information about \code{vachette_data} object
#'
#' @param x An \code{vachette_data}.
#' @param log_file character; File path to direct console output e.g., \code{"log.txt"}
#' @param ... Additional args.
#' @return Returns \code{x} invisibly.
#' @export
summary.vachette_data <- function(x, log_file = NULL, ...) {
  stopifnot(inherits(x, "vachette_data"))

  if (!is.null(log_file)) {
    if (file.exists(log_file))  {
      res <-
        askYesNo(msg = paste0(log_file, " exists and will be overwritten"))
      if (res) {
        unlink(log_file)
      } else {
        stop("Cannot delete log file. Set argument `log_file`=NULL to ignore.")
      }
    }
    conn <- file(log_file, open = "at")
    assign("vachette_log_file", value = conn, envir = vachette_env)
    on.exit(close(conn))
    on.exit(assign("vachette_log_file", value = NULL, envir = vachette_env),
            add = TRUE)
  }

  log_output(paste0(sprintf("Model name:\t\t%s", x$model.name), "\n"))
  log_output(paste0(sprintf("Covariate reference:\t%s", paste0(
    paste0(names(x$covariates), "=", x$covariates), collapse = " , "
  )), "\n"))
  log_output(paste0(sprintf(
    "Unique covariate combinations:\t\t%s",
    x$summary$values$n_iter
  ),
  "\n")
  )

  log_output("\n--------------------------------------------\n")
  log_output("vachette transformations")
  log_output("\n--------------------------------------------\n")

  log_output("\nAsymptotes:\n")
  log_output(paste0(sprintf("\tAsymptote right:\t%s", x$summary$values$asymptote_right),
             "\n"))
  log_output(paste0(sprintf("\tAsymptote left:\t\t%s", x$summary$values$asymptote_left),
             "\n"))
  log_output(
    paste0(sprintf(
      "\tZero asymptote right:\t%s",
      x$summary$values$zero_asymptote_right
    ),
    "\n"
  )
  )
  log_output(
    paste0(sprintf(
      "\tZero asymptote left:\t%s",
      x$summary$values$zero_asymptote_left
    ),
    "\n"
  )
  )


  ucov <- x$tab.ucov

  errors <- x$summary$errors

  if (length(errors) > 0) {
    log_output("\nErrors:")
    for (i in seq_along(errors)) {
      emsg <- errors[[i]]
      if (any(nzchar(emsg))) {
        log_output(paste0("\n\t---i.ucov = ", i, "---", "\n"))
        if (!is.null(log_file)) {
          log_output(capture.output(print.data.frame(ucov[i, ])))
        } else {
          print.data.frame(ucov[i, ])
        }
        log_output(emsg)
      }
    }
  } else {
    log_output("\nObservations:\n")
    log_output(paste0(sprintf(
      "\tN transformed:\t\t%s",
      nrow(x$obs.all) - nrow(x$obs.excluded)
    ), "\n"))
    log_output(paste0(sprintf("\tN excluded:\t\t%s", nrow(x$obs.excluded)), "\n"))
    log_output("\nObservations Excluded:\n")
    if (!is.null(log_file)) {
      log_output(capture.output(print.data.frame(x$obs.excluded %>% dplyr::select(-c("REP", "COV", "PRED", "exclude")))))
    } else {
      print.data.frame(x$obs.excluded %>% dplyr::select(-c("REP", "COV", "PRED", "exclude")))
    }
  }

  warnings <- x$summary$warnings

  if (length(warnings) > 0) {
    log_output("\n\nWarnings:")
    for (i in seq_along(warnings)) {
      wmsg <- warnings[[i]]
      if (any(nzchar(wmsg))) {
        log_output(paste0("\n\t---i.ucov = ", i, "---", "\n"))
        if (!is.null(log_file)) {
          log_output(capture.output(print.data.frame(ucov[i, ])))
        } else {
          print.data.frame(ucov[i, ])
        }
        log_output(wmsg)
      }
    }
  }
  # output error type
  invisible(x)
}
