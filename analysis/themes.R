# Themes that can be used for plotting for lung project




themeAnno_Job <- ggplot2::theme(
  legend.position = 'right',
  axis.ticks = ggplot2::element_blank(),
  axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
  axis.text.x = ggplot2::element_blank(),
  text = ggplot2::element_text(size = 8, family='Helvetica', face = 'bold'),
  panel.background = ggplot2::element_rect(fill = NA),
  panel.grid.major = ggplot2::element_line(NULL),
  panel.grid.major.x = ggplot2::element_line(NULL),
  panel.grid.major.y = ggplot2::element_line(NULL),
  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
  legend.text = ggtext::element_markdown(),
  legend.title = ggtext::element_markdown(),
  legend.key = ggplot2::element_blank()
)

themeTrack_Job <- ggplot2::theme(
  legend.position = 'right',
  legend.direction = 'horizontal',
  text = ggplot2::element_text(size = 7, family = 'Helvetica', face = 'bold'),
  axis.title.y = ggtext::element_textbox_simple(size = 6, orientation = 'left-rotated', width = NULL, halign = .5),
  axis.text.x = ggplot2::element_blank(),
  axis.ticks.x = ggplot2::element_blank(),
  panel.grid.major.x = ggplot2::element_blank(),
  panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
  panel.grid.minor.y = ggplot2::element_blank(),
  axis.title.x = ggplot2::element_blank(),
  panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
  panel.border = ggplot2::element_rect(fill = NA, colour = NA),
  strip.background = ggplot2::element_rect(colour = 'grey20', fill = 'white'),
  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
  legend.text = ggtext::element_markdown(),
  legend.title = ggtext::element_markdown(),
  legend.key = ggplot2::element_blank()
)

theme_Job <- ggplot2::theme(
  legend.position = 'bottom',
  legend.direction = 'horizontal',
  text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
  axis.text.x = ggtext::element_markdown(),
  axis.text.y = ggtext::element_markdown(),
  axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
  axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
  strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
  panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
  panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
  panel.grid.minor.y = ggplot2::element_blank(),
  panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
  panel.border = ggplot2::element_rect(fill = NA, colour = NA),
  strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
  legend.text = ggtext::element_markdown(),
  legend.title = ggtext::element_markdown(),
  legend.key = ggplot2::element_blank()
)

# Half-Half Boxplot(s)
geomJob_HalfHalfBox <- list(
  gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, color = 'black'),
  gghalves::geom_half_point_panel(side = 'r', color = 'black', shape = 21, transformation = ggbeeswarm::position_quasirandom(width = .1, groupOnX = T)),
  ggplot2::stat_summary(fun = median, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, hjust = .5, vjust = -15, angle = 90)
)

geomJob_HalfHalfBox_Small <- list(
  gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, color = 'black'),
  gghalves::geom_half_point(side = 'r', color = 'black', alpha = .66, size = .75, shape = 21, transformation = ggbeeswarm::position_quasirandom(width = .075)),
  ggplot2::stat_summary(fun = min, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, hjust = 1.25, angle = 90)
)