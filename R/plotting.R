
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


# For same looks and feel of plots:
render <- theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.title = element_text(size=14),
        plot.subtitle=element_text(size=12))

# Multiple plots on one page:
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
