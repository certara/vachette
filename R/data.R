#' Oral Absorption Observed Data
#'
#' This dataset contains a single covariate \code{WT} and a single dose.
#' The dose event row is recorded by the \code{EVID} column.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{\code{ID}}{Subject ID}
#'   \item{\code{time}}{Time used as x-variable}
#'   \item{\code{EVID}}{Dose event column}
#'   \item{\code{DV}}{Dependent variable}
#'   \item{\code{PRED}}{Population predicted value}
#'   \item{\code{WT}}{Weight covariate}
#' }
#' @examples
#' \dontrun{
#' data(oral_absorption_obs)
#' head(oral_absorption_obs)
#' }
"oral_absorption_obs"

#' Oral Absorption Typical Curves Data
#'
#' This dataset contains a single covariate \code{WT} and a single dose.
#' The dose event row is recorded by the \code{EVID} column.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{\code{ID}}{Subject ID}
#'   \item{\code{time}}{Time used as x-variable}
#'   \item{\code{EVID}}{Dose event column}
#'   \item{\code{PRED}}{Population predicted value}
#'   \item{\code{WT}}{Weight covariate}
#' }
#' @examples
#' \dontrun{
#' data(oral_absorption_typ)
#' head(oral_absorption_typ)
#' }
"oral_absorption_typ"

#' Oral Absorption Simulated Data
#'
#' Simulated dataset containing 100 replicates. Time points are
#' simulated given observation time points in  \code{\link{oral_absorption_obs}}.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{\code{REP}}{Replicate number}
#'   \item{\code{ID}}{Subject ID}
#'   \item{\code{time}}{Time used as x-variable}
#'   \item{\code{EVID}}{Dose event column}
#'   \item{\code{DV}}{Dependent variable}
#'   \item{\code{PRED}}{Population predicted value}
#'   \item{\code{WT}}{Weight covariate}
#' }
#' @examples
#' \dontrun{
#' data(oral_absorption_sim)
#' head(oral_absorption_sim)
#' }
"oral_absorption_sim"
