vachette_env <- new.env()

vachette_env$summary_out <- list(errors = list(),
                    warnings = list(),
                    values = list()
)

vachette_env$warning_out <- NULL

log_output <- function(message) {

  if (exists("vachette_log_file", envir = vachette_env)){
    log_file <- get("vachette_log_file", envir = vachette:::vachette_env)
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
#' @param ... Additional args.
#' @return Returns \code{x} invisibly.
#' @export
summary.vachette_data <- function(x, ...) {
  stopifnot(inherits(x, "vachette_data"))

  cat(sprintf("Model name:\t\t%s", x$model.name), "\n")
  cat(sprintf("Covariate reference:\t%s", paste0(paste0(names(x$covariates), "=", x$covariates), collapse = " , ")),"\n")
  cat(sprintf("Unique covariate combinations:\t\t%s", x$summary$values$n_iter), "\n")

  errors <- x$summary$errors

  if(length(errors) > 0){
    cat("\nErrors:\n")
    for (i in seq_along(errors)){
      emsg <- errors[[i]]
      if (any(nzchar(emsg))) {
        cat("\n\t---i.ucov = ", i, "---", "\n")
        cat(emsg)
      }
    }
    return(invisible(x))
  }

  cat("\n--------------------------------------------\n")
  cat("vachette transformations")
  cat("\n--------------------------------------------\n")

  cat("\nAsymptotes:\n")
  cat(sprintf("\tAsymptote right:\t%s", x$summary$values$asymptote_right), "\n")
  cat(sprintf("\tAsymptote left:\t\t%s", x$summary$values$asymptote_left), "\n")
  cat(sprintf("\tZero asymptote right:\t%s", x$summary$values$zero_asymptote_right), "\n")
  cat(sprintf("\tZero asymptote left:\t%s", x$summary$values$zero_asymptote_left), "\n")

  cat("\nObservations:\n")
  cat(sprintf("\tN transformed:\t\t%s", nrow(x$obs.all) - nrow(x$obs.excluded)), "\n")
  cat(sprintf("\tN excluded:\t\t%s", nrow(x$obs.excluded)), "\n")


  cat("\nObservations Excluded:\n")
  print.data.frame(x$obs.excluded %>% dplyr::select(-c("REP", "COV", "PRED", "exclude")))

  warnings <- x$summary$warnings

  if(length(warnings) > 0){
    cat("\nWarnings:\n")
    for (i in seq_along(warnings)){
      wmsg <- warnings[[i]]
      if (any(nzchar(wmsg))) {
        cat("\n\t---i.ucov = ", i, "---", "\n")
        cat(wmsg)
      }
    }
  }
  # output error type
  invisible(x)
}
