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
    message(message, appendLF = TRUE)
  }
}



#' Summarise \code{vachette_data}
#'
#' Summary generic used to return information about \code{vachette_data} object
#'
#' @param object \code{vachette_data} object.
#' @param verbose logical; Set to \code{TRUE} to view additional warnings output
#' @param trim logical; If \code{TRUE}, only the first 20 rows are printed for \code{data.frame} output. Set to \code{FALSE} to display all rows
#' @param log_file character; File path to direct console output e.g., \code{"log.txt"}
#' @param ... additional args.
#' @return \code{object} invisibly.
#' @export
summary.vachette_data <-
  function(object,
           verbose = FALSE,
           trim = TRUE,
           log_file = NULL,
           ...) {
    stopifnot(inherits(object, "vachette_data"))

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
    n_ucov <- nrow(object$tab.ucov)
    n_obs_excl <-  nrow(object$obs.excluded)
    if (!trim) {
      old_max_print <- getOption("max.print")
      options(max.print = 1000000)
      on.exit(options(max.print = old_max_print), add = TRUE)
    }

    log_output(paste0(sprintf("Model name:\t\t%s", object$model.name), "\n"))
    log_output(paste0(sprintf(
      "Covariate reference:\t%s", paste0(paste0(names(object$covariates), "=", object$covariates), collapse = " , ")
    ), "\n"))
    log_output(paste0(
      sprintf("Unique covariate combinations:\t\t%s",
              n_ucov),
      "\n"
    ))

    if (trim) {
      ucov_data <- head(object$tab.ucov, 20)
    } else {
      ucov_data <- object$tab.ucov
    }

    if (!is.null(log_file)) {
      log_output(capture.output(print.data.frame(ucov_data)))
    } else {
      print.data.frame(ucov_data)
    }

    if (trim && n_ucov > 20) {
      log_output(paste0("[set argument trim = FALSE -- omitted ", n_ucov - 20, " rows]"))
    }

    if (!is.null(object$summary)){

    log_output("\n--------------------------------------------\n")
    log_output("vachette transformations")
    log_output("\n--------------------------------------------\n")

    log_output("\nAsymptotes:")
    log_output(paste0(
      sprintf("\tAsymptote right:\t%s", object$summary$values$asymptote_right)
    ))
    log_output(paste0(
      sprintf("\tAsymptote left:\t\t%s", object$summary$values$asymptote_left)
    ))
    log_output(paste0(
      sprintf(
        "\tZero asymptote right:\t%s",
        object$summary$values$zero_asymptote_right
      )
    ))
    log_output(paste0(
      sprintf(
        "\tZero asymptote left:\t%s",
        object$summary$values$zero_asymptote_left
      )
    ))

    errors <- object$summary$errors

    if (length(errors) > 0) {
      log_output("\nErrors:")
      for (i in seq_along(errors)) {
        emsg <- errors[[i]]
        if (any(nzchar(emsg))) {
          log_output(paste0("\n\t---i.ucov = ", i, "---", "\n"))
          if (!is.null(log_file)) {
            log_output(capture.output(print.data.frame(object$tab.ucov[i,])))
          } else {
            print.data.frame(object$tab.ucov[i,])
          }
          log_output(emsg)
        }
      }
    } else {
      log_output("\nObservations:")
      log_output(paste0(sprintf(
        "\tN transformed:\t\t%s",
        nrow(object$obs.all) - nrow(object$obs.excluded)
      )))
      log_output(paste0(sprintf(
        "\tN excluded:\t\t%s", nrow(object$obs.excluded)
      )))
      log_output("\nObservations Excluded:\n")
      obs_excl <-
        object$obs.excluded %>% dplyr::select(-c("REP", "COV", "PRED", "exclude"))
      if (trim) {
        obs_excl <- head(obs_excl, 20)
      }
      if (!is.null(log_file)) {
        log_output(capture.output(print.data.frame(obs_excl)))
      } else {
        print.data.frame(obs_excl)
      }
      if (trim && n_obs_excl > 20) {
        log_output(paste0(
          "[set argument trim = FALSE -- omitted ",
          n_obs_excl - 20,
          " rows]"
        ))
      }
    }

    warnings <- object$summary$warnings

    if (length(warnings) > 0 && verbose) {
      log_output("\n\nWarnings:")
      for (i in seq_along(warnings)) {
        wmsg <- warnings[[i]]
        if (any(nzchar(wmsg))) {
          log_output(paste0("\n\t---i.ucov = ", i, "---", "\n"))
          if (!is.null(log_file)) {
            log_output(capture.output(print.data.frame(object$tab.ucov[i,])))
          } else {
            print.data.frame(object$tab.ucov[i,])
          }
          log_output(wmsg)
        }
      }
    }
   }
    # output error type
    invisible(object)
  }
