#'Remove the dependent variables that are below the threshold of
#'sparsity when the data type is count data in longdat_cont()
#' @param values Internal function argument.
#' @param data Internal function argument.
#' @param nonzero_count_cutoff1 Internal function argument.
#' @param nonzero_count_cutoff2 Internal function argument.
#' @param theta_cutoff Internal function argument.
#' @param Ps_null_model Internal function argument.
#' @param prevalence Internal function argument.
#' @param absolute_sparsity Internal function argument.
#' @param mean_abundance Internal function argument.
#' @param p_poho Internal function argument.
#' @param assoc Internal function argument.
#' @importFrom rlang .data
#' @importFrom stats as.formula confint cor.test kruskal.test
#'             na.omit p.adjust wilcox.test
#' @importFrom magrittr '%>%'
#' @import tibble
#' @import dplyr
#' @name rm_sparse_cont
utils::globalVariables(c("NB_theta"))

rm_sparse_cont <- function(values, data, nonzero_count_cutoff1,
                           nonzero_count_cutoff2, theta_cutoff,
                           Ps_null_model, prevalence,
                           absolute_sparsity, mean_abundance,
                           p_poho, assoc) {

  absolute_sparsity <- c()
  for (i in 1:ncol(values)) {
    absolute_sparsity[i] <- sum(values[ , i] == 0)
  }
  non_zero_count <- nrow(data) - absolute_sparsity
  Ps_null_model <- Ps_null_model %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(non_zero_count = non_zero_count) %>%
    tibble::column_to_rownames()

  bac_exclude_1 <-
    subset(Ps_null_model, non_zero_count <= nonzero_count_cutoff1 &
             NB_theta >= theta_cutoff)
  bac_exclude_2 <-
    subset(Ps_null_model, non_zero_count <= nonzero_count_cutoff2)
  bac_exclude <- unique(c(rownames(bac_exclude_1), rownames(bac_exclude_2)))
  bac_include <-
    rownames(Ps_null_model)[!rownames(Ps_null_model) %in% bac_exclude]

  prevalence <- prevalence[match(bac_include, table = rownames(Ps_null_model))]
  absolute_sparsity <- absolute_sparsity[match(bac_include,
                                               table = rownames(Ps_null_model))]
  mean_abundance <- mean_abundance[match(bac_include,
                                         table = rownames(Ps_null_model))]

  Ps_null_model <- Ps_null_model %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include) %>%
    tibble::column_to_rownames()

  p_poho <- p_poho %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include) %>%
    tibble::column_to_rownames()

  assoc <- assoc %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include) %>%
    tibble::column_to_rownames()

  variables <- bac_include
  return(list(prevalence, absolute_sparsity, mean_abundance,
              Ps_null_model, p_poho, assoc, variables, bac_include))
}
