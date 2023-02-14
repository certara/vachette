
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


#' p.scaled.typical.curves.landmarks
#'
#' @param vachette_data Object of class\code{vachette_data}
#'
#' @export
p.scaled.typical.curves.landmarks <- function(vachette_data) {

  curves.all <- vachette_data$curves.all
  lm.all <- vachette_data$lm.all
  model.name <- vachette_data$model.name
  xstop <- vachette_data$xstop

  stopifnot(!is.null(curves.all))

  curves.all %>%
  ggplot(aes(x=x,y=y,group=ucov,
             col = factor(seg)))+
  geom_line(lwd=1)+
  geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
  # Add landmark positions
  geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)+
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title=paste0(model.name," - Typical curve segments"),
       subtitle = "Dashed line: Reference typical curve",
      # caption=script,
       col="Segment")+
  render
}

#' p.scaled.typical.full.curves.landmarks
#'
#' @param vachette_data Object of class\code{vachette_data}
#'
#' @export
p.scaled.typical.full.curves.landmarks <- function(vachette_data) {
  curves.all <- vachette_data$curves.all
  lm.all <- vachette_data$lm.all
  model.name <- vachette_data$model.name

  curves.all %>%
  ggplot(aes(x=x,y=y,group=ucov,
             col = factor(seg)))+
  geom_line(lwd=1)+
  geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
  # Add landmark positions
  geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)+
  labs(title=paste0(model.name," - Typical curve segments"),
       subtitle = "Dashed line: Reference typical curve",
      # caption=script,
       col="Segment")+
  render
}

#' p.scaling.factor
#'
#' @param vachette_data Object of class\code{vachette_data}
#'
#' @export
p.scaling.factor <- function(vachette_data) {
  curves.all <- vachette_data$curves.all
  obs.all <- vachette_data$obs.all
  model.name <- vachette_data$model.name

 curves.all %>%
  mutate(mycurve = ifelse(ref=='Yes','Reference','Query')) %>%
  ggplot(aes(x=x,y=x.scaling,col=factor(seg)))+
  geom_line(lwd=1)+
  facet_wrap(~paste(mycurve," covariate =",COV))+
  coord_cartesian(xlim=c(NA,max(obs.all$x,obs.all$x.scaled)),
                  ylim=c(0,max(curves.all$x.scaling[curves.all$x<=max(obs.all$x,obs.all$x.scaled)])))+
  labs(title=paste0(model.name," - x-scaling factors"),
       subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error",if(vachette_data$PROP_TR) "Proportional Error"),
       #caption=script,
       col="Segment") +
  render
}
