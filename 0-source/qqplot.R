qqplot <- function(p,
                   filename = NULL,
                   width = 4,
                   height = 4){
  require(ggtext)
  df <- data.frame(o = -log10(sort(p,decreasing=FALSE)),
                   e = -log10(stats::ppoints(length(p))))
  
  plot <- df |> 
    ggplot(aes(x = e,
               y = o)) +
    geom_point(size = 0.5) +
    geom_abline(intercept = 0, slope = 1,
                colour = "grey40") +
    labs(x = "Expected -log<sub>10</sub>(p)", 
         y = "Observed -log<sub>10</sub>(p)") + 
    theme_minimal() +
    theme( 
      legend.position = "top",
      axis.title.y = element_markdown(),
      axis.title.x = element_markdown()
    )
  
  if(!is.null(filename)){
    ggsave(plot = plot, 
           filename = filename, width = width, height = height)
  } else {
    return(plot)
  }
  }