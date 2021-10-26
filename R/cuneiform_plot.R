#' Creates cuneiform plots of result table from longdat_disc() or
#' longdat_cont()
#' @param result_table The result table from longdat_disc() or
#' longdat_cont() output, or any table that has the same format.
#' @param x_axis_order The plotting order of the x axis.
#' It should be a character vector
#' (eg. c("Effect_1_2", "Effect_2_3", "Effect_1_3")).
#' @param pos_color The color for a positive effect size.
#'  It should be a hex color code (eg. "#b3e6ff") or the colors recognized
#'  by R (see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf).
#'  The default is "red".
#' @param neg_color The color for a negative effect size.
#'  It should be a hex color code (eg. "#b3e6ff") or the colors recognized
#'  by R (see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf).
#'  The default is "blue".
#' @param panel_width The width of the effect size panel on the left
#' relative to the confounding status panel on the right (width set to 1).
#' It should be a numerical vector. The default is 4.
#' @param title The name of the plot title. The default is
#' "LongDat result cuneiform plot".
#' @param title_size The size of the plot title. The default is 20.
#' @param confound_text_size The size of the text in the confounding status
#' panel. The default is 4.
#' @param x_label_size The size of the x label. The default is 10.
#' @param y_label_size The size of the y label. The default is 10.
#' @param legend_title_size The size of the legend title. The default is 12.
#' @param legend_text_size The size of the legend text The default is 10.
#' @export
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import patchwork
#' @importFrom rlang .data
#' @importFrom magrittr '%>%'
#' @name cuneiform_plot
#' @return a 'ggplot' object
#' @details
#' This function creates a cuneiform plot which displays the
#' result of longdat_disc() or longdat_cont(). It plots the effect sizes
#'  within each time interval for each feature, and also shows the confounding
#'  status. Only the features with non-NS signals will be included in the plot.
#'  The output is a ggplot object in patchwork structure. For further
#'  customization of the plot, please refer to the vignette.
#' @examples
#'\dontrun{
#' # Get the path of example dataset
#' system.file("Fasting_disc.txt", package = "LongDat")
#'
#' # Paste the directory to the input below
#' test_disc <- longdat_disc(input = "your_path_to/Fasting_disc.txt",
#' data_type = "count", test_var = "Time_point", variable_col = 7,
#' fac_var = c(1:3))
#'
#' test_plot <- cuneiform_plot(result_table = test_disc[[1]])
#'}
#utils::globalVariables(c(""))

cuneiform_plot <- function(result_table,
                           x_axis_order = NULL,
                           pos_color = "red",
                           neg_color = "blue",
                           panel_width = 4,
                           title = "LongDat result cuneiform plot",
                           title_size = 20,
                           confound_text_size = 4,
                           x_label_size = 10,
                           y_label_size = 10,
                           legend_title_size = 12,
                           legend_text_size = 10
                           ) {
  if (missing(result_table)) {
    stop('Error! Necessary argument "result_table" is missing.')}

  result_table <- result_table
  sig_result <- result_table %>%
    dplyr::filter(Signal != "NS")

  if (nrow(sig_result) == 0) {
    stop('All results are insignificant! There is nothing to plot.')}

  # Select the required columns
  sig_wide <- sig_result %>%
    dplyr::select(c(Signal, stringr::str_which(string = colnames(sig_result),
                                               pattern = "Effect"))) %>%
    tibble::rownames_to_column("Variables")

  # Divide them into 2 tables: Effect and EffectSize
  Effect_wide <- sig_wide %>%
    dplyr::select(stringr::str_which(string = colnames(sig_wide),
                                     pattern = "EffectSize",
                                     negate = T))
  EffectSize_wide <- sig_wide %>%
    dplyr::select(c(Variables, Signal,
                    stringr::str_which(string = colnames(sig_wide),
                                       pattern = "EffectSize")))

  # Transform into long format
  Effect_long <- Effect_wide %>%
    tidyr::pivot_longer(stringr::str_which(string = colnames(Effect_wide),
                                           pattern = "Effect"),
                        names_to = "Effect_name", values_to = "Effect")
  EffectSize_long <- EffectSize_wide %>%
    tidyr::pivot_longer(stringr::str_which(string = colnames(EffectSize_wide),
                                           pattern = "EffectSize"),
                        names_to = "EffectSize_name", values_to = "EffectSize")

  # Merge the two long tables
  All_long <- cbind(Effect_long, EffectSize_long[ , c(3:4)])

  # Define plotting parameters
  All_long$Alpha <- ifelse(All_long$Effect != "NS",
                           yes = "Significant", no = "Insignificant")
  All_long$Shape <- ifelse(All_long$EffectSize > 0,
                           yes = "24",
                           no = ifelse(All_long$EffectSize < 0,
                                       yes = "25", no = "1"))

  # Reorder plotting sequence
  if (is.null(x_axis_order)) {
    All_long$Effect_name <-  All_long$Effect_name
  } else {
    All_long$Effect_name <- factor(All_long$Effect_name,
                                   levels = x_axis_order)
  }

  # Plotting
  g1 <- ggplot2::ggplot(All_long, aes(x = Effect_name, y = Variables)) +
    geom_point(aes(shape = Shape, fill = EffectSize,
                   alpha = Alpha), size = 3.5) +
    scale_shape_manual(values = c(1, 24, 25)) +
    scale_fill_gradient2(midpoint = 0, low = neg_color, mid = "white",
                         high = pos_color, n.breaks = 8,
                         limits = c(-1, 1) * max(abs(All_long$EffectSize))) +
    scale_alpha_manual(breaks = c("Insignificant", "Significant"),
                       values=c(0.4, 1)) +
    ggtitle(title) +
    labs(fill = "Effect size", alpha = "Significance") +
    theme_light() +
    theme(title = element_text(size = title_size),
          axis.text.y = element_text(size = y_label_size),
          axis.title.x=element_blank(),
          axis.text.x=element_text(size = x_label_size),
          axis.title.y=element_blank(),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size)) +
    guides(shape = "none") # Remove shape legend

  g2 <- ggplot2::ggplot(Effect_wide, aes(x = "Confounding status",
                                         y = Variables)) +
    geom_text(aes(label = Signal),
              size = confound_text_size, color = "gray30") +
    theme_light() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(size = x_label_size),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.x = element_blank())

  final_plot <- (g1|g2) + patchwork::plot_layout(guides = "collect",
                                 widths = c(panel_width, 1))
  return(final_plot)
  print("Finished plotting successfully!")
}
