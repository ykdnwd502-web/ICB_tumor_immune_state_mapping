############################################################
## 00_plot_helpers.R
##
## Shared plotting utilities
##
## Loaded before all figure scripts
############################################################


if(!exists("safe_ggsave")){
  
  safe_ggsave <- function(
    filename,
    plot,
    width,
    height,
    dpi = 300
  ){
    
    dir.create(
      dirname(filename),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      units = "in"
    )
    
    message(
      "Saved: ",
      filename
    )
    
  }
  
}