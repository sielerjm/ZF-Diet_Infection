# Plot Settings



# Color Palette -----------------------------------------------------------
# Source: RcolorBrewer: https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/


## Sequential (Ordinal)--------------------------------------------------------------
#   Sequential palettes (first list of colors), which are suited to ordered data 
#     that progress from low to high (gradient).

pal.BuPu <- brewer.pal(9, "BuPu") 
pal.Greys <- brewer.pal(9, "Greys")


## Qualitative (categorical,  non-ordinal) ---------------------------------
#   Qualitative palettes (second list of colors), which are best suited to represent nominal or categorical data. 
#     They do not imply magnitude differences between groups. 

pal.Set1 <- brewer.pal(9, "Set1")  # Not colorblind friendly
pal.Set2 <- brewer.pal(8, "Set2")
pal.Dark2 <- brewer.pal(8, "Dark2")
pal.Paired <- brewer.pal(12, "Paired")


## Diverging (Numeric, spectrum) ------------------------------------------------------
#   Diverging palettes (third list of colors), which put equal emphasis on mid-range 
#     critical values and extremes at both ends of the data range.

pal.RdYlGn <- brewer.pal(11, "RdYlGn")
pal.Spectral <- brewer.pal(9, "Spectral")
pal.BrBg <- brewer.pal(11, "BrBG")


# GGplot2 Theme -----------------------------------------------------------
#     Standardizing ggplot settings
#     (adapted from Keaton Stagaman)

my_theme <- theme_update(
  legend.position = "top",
  legend.box = "vertical",
  legend.box.just = "center",
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 12),
  legend.key = element_rect(fill = "white"),
  legend.key.size = unit(1, "line"), # legend symbol size
  legend.spacing.y = unit(-.25, 'cm'),
  
  strip.text = element_text(size = 12),
  
  plot.caption = element_text(hjust = 0, size = 12),
  
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 16),
  
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  panel.background = element_rect(fill = pal.Greys[1]),
  panel.grid.major = element_line(color = pal.Greys[3]),
  panel.grid.minor = element_line(color = pal.Greys[3]),
) 



################################################################################
################################################################################