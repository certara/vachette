
#' Plot theme for vachette
#'
#' A ggplot2 theme for vachette plot output
#'
#' @return \code{theme, gg}
#' @export
#'
#' @examples
#' theme_vachette()
#'
theme_vachette <- function() {
  ggplot2::theme_bw()+
  ggplot2::theme(text = ggplot2::element_text(size=16),
        axis.text.x = ggplot2::element_text(size=14),
        axis.text.y = ggplot2::element_text(size=14),
        plot.title = ggplot2::element_text(size=14),
        plot.subtitle=ggplot2::element_text(size=12))

}
