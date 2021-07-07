# https://rpubs.com/Koundy/71792

# add new fonts
#install.packages('showtext', dependencies = TRUE)
#suppressPackageStartupMessages(library(showtext))
#showtext::showtext_auto()

# https://fonts.google.com/
#font_add_google( "Roboto Condensed", "RobotoCondensed-Regular")
#font_add_google("IBM Plex Serif", "IBM Plex Serif")
# Check the current search path for fonts
#font_paths()  

# syntax: font_add(family = "<family_name>", regular = "/path/to/font/file")
# font_add("Comic Sans MS", "comic.ttf")
#font_families()

library(extrafont)
# import fonts - only once
# font_import()
# load fonts - every session
loadfonts(device = "win", quiet = TRUE)

# find the name of a font you need for the family parameter of element_text
#library(tidyverse)
#fonttable() %>% filter(if_any(everything(), ~str_detect(.,"Comic")))

#fonttable() %>% dplyr::filter(str_detect(fontfile,"IBM")) %>% dplyr::pull(FamilyName)
#fonttable() %>% dplyr::filter(str_detect(fontfile,"Roboto")) %>% dplyr::pull(FamilyName)

# Automatically use showtext for new devices
#showtext_auto()
# Turn off if no longer needed
#showtext_auto(FALSE)

library(hrbrthemes)
library(ggthemes)

# Using ggthemr package https://github.com/Mikata-Project/ggthemr
# set ggthemr theme
#ggthemr("<theme name>") 
# plot your existing figure with the new theme
#plt
# to remove all ggthemr effects later:
#ggthemr_reset()


# colors
#scales::show_col(c("#AD6F3B", "indianred3","#673770", "#8569D5", "#5E738F","coral4",
#"lightskyblue4", "firebrick4", "rosybrown1",
#"#CBD588", "orange","#DA5724", "#508578","#652926", "#C84248", "#D1A33D", "#8A7C64", "#599861","#5F7FC7","#CD9BCD"))
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


#---------------------

theme_Publication_1 <- function(base_size = 14,
                                strip_text_size = 11,
                                strip_text_margin = 5,
                                subtitle_size = 13,
                                subtitle_margin = 10,
                                plot_title_size = 16,
                                plot_title_margin = 10,
                                title_axis_size = rel(1), # relative to the default
                                axis_text_size = 11,
                                ...) {
  
  ret <- ggplot2::theme_minimal(base_family = "helvetica",
                                base_size = base_size, ...)
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    # face = 'bold',
    margin = ggplot2::margin(b = strip_text_margin),
    family = "helvetica"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "helvetica"
  )
  ret$plot.title <- ggplot2::element_text(
    face = "bold",
    hjust = 0.5, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "helvetica"
  )
  ret$axis.title <- ggplot2::element_text(
    face = "bold",
    size = title_axis_size,
    family = "helvetica"
  )
  ret$axis.text <- ggplot2::element_text(
    face = "bold", 
    size = axis_text_size,
    family = "helvetica" 
  )
  ret
  
}

theme_Publication <- function(base_size=14, base_family="helvetica",...) {
  library(grid)
  library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.border = element_rect(fill = "transparent",colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(face = "bold",angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
            axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
          #  plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold", size = 11)
    ))
 
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# using the theme_article function form egg package

theme_Publication_2 <- function(...){

  egg::theme_article() +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 11), 
          axis.text.y = element_text(face = "bold", size = 11),
          strip.text = element_text(face="bold", size = 11),
          strip.background = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.margin = unit(0, "cm"),
          legend.title = element_text(face="italic", size = 12),
          legend.text = element_text(size = 11),
          text = element_text(family = "Garamond", color = "grey20"),
          legend.background = element_blank(),
          plot.margin=unit(c(10,5,5,5),"mm"))
 
  
}



#' Minimal ggplot2 theme using the Roboto Condensed and Roboto Bold fonts
#'
#' @param base_size base font size
#' @param strip_text_size,strip_text_margin plot strip text size and margin
#' @param subtitle_size,subtitle_margin plot subtitle size and margin
#' @param plot_title_size,plot_title_margin plot title size and margin
#' @param ... Other arguments passed to \code{theme_minimal}
#'
#' @details The Roboto Condensed and Roboto Bold fonts are both Google fonts;
#' they can be found at \url{https://fonts.google.com/specimen/Roboto+Condensed}
#' and \url{https://fonts.google.com/specimen/Roboto}. These fonts must be
#' installed locally on your computer for this theme to work.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     labs(title = "A Lovely Plot",
#'          subtitle = "What can the subtitle tell us?") +
#'     theme_roboto()
#'}
#'
#' @export

# for windows
theme_roboto <- function(base_size = 11,
                         strip_text_size = 12,
                         strip_text_margin = 5,
                         subtitle_size = 13,
                         subtitle_margin = 10,
                         plot_title_size = 16,
                         plot_title_margin = 10,
                         ...) {
  # Automatically use showtext for new devices
  
  ret <- ggplot2::theme_minimal(base_family = "Roboto",
                                base_size = base_size, ...)
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    margin = ggplot2::margin(b = strip_text_margin),
    family = "Roboto",
    face = 'bold'
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "Roboto",
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "Roboto",
    face = 'bold'
  ) 
  ret
  # Turn off if no longer needed
 
}




#' Minimal ggplot2 theme using the IBM Plex Sans fonts
#'
#' @param base_size base font size
#' @param strip_text_size,strip_text_margin plot strip text size and margin
#' @param subtitle_size,subtitle_margin plot subtitle size and margin
#' @param plot_title_size,plot_title_margin plot title size and margin
#' @param ... Other arguments passed to \code{theme_minimal}
#'
#' @details The IBM Plex fonts are open source and can be found at
#' \url{https://ibm.github.io/type/}. These fonts must be installed locally on
#' your computer for this theme to work.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     labs(title = "A Lovely Plot",
#'          subtitle = "What can the subtitle tell us?") +
#'     theme_plex()
#'
#' ggplot(diamonds, aes(carat, price, color = clarity)) +
#'     geom_point(alpha = 0.7) +
#'     facet_wrap(~cut) +
#'     labs(title = "A Lovely Plot",
#'          subtitle = "What can the subtitle tell us?") +
#'          theme_plex()
#'
#'}
#'
#' @export
theme_plex <- function(base_size = 11,
                       strip_text_size = 12,
                       strip_text_margin = 5,
                       subtitle_size = 13,
                       subtitle_margin = 10,
                       plot_title_size = 16,
                       plot_title_margin = 10,
                       ...) {

  ret <- ggplot2::theme_minimal(base_family = "IBM Plex Sans",
                                base_size = base_size, ...)
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    margin = ggplot2::margin(b = strip_text_margin),
    family = "IBM Plex Sans"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "IBM Plex Sans Sans"
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "IBM Plex Sans Sans"
  ) 
  ret
 
}


theme_Publication_3 <- function(base_size = 12,
                       strip_text_size = 12,
                       strip_text_margin = 5,
                       subtitle_size = 13,
                       subtitle_margin = 10,
                       plot_title_size = 16,
                       plot_title_margin = 10,
                       ...) {

  ret <- ggplot2::theme_minimal(base_family = "Comic Sans MS",
                                base_size = base_size, ...) +
    theme(panel.border = element_rect(fill = "NA",colour = "grey80")
    )
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
   # face = 'bold',
    margin = ggplot2::margin(b = strip_text_margin),
    family = "Comic Sans MS"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "Comic Sans MS"
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "Comic Sans MS"
  ) 
  ret

}



