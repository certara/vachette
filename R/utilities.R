vachette_env <- new.env()

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
