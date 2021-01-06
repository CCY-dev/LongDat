#'Remove the dependent variables that are below the threshold of sparsity when the data type is count data in longdat_disc()

rm_sparse_disc <- function(values, data, nonzero_count_cutoff1, nonzero_count_cutoff2, theta_cutoff, Ps_null_model,
                      prevalence, absolute_sparsity, mean_abundance, Ps_poho_fdr, delta) {
  absolute_sparsity <- c()
  for (i in 1:ncol(values)) {
    absolute_sparsity[i] <- sum(values[ , i] == 0)
  }
  non_zero_count <- nrow(data) - absolute_sparsity
  Ps_null_model <- Ps_null_model %>%
    rownames_to_column() %>%
    mutate(non_zero_count = non_zero_count) %>%
    column_to_rownames()

  bac_exclude_1 <- subset(Ps_null_model, non_zero_count <= nonzero_count_cutoff1 & NB_theta >= theta_cutoff)
  bac_exclude_2 <- subset(Ps_null_model, non_zero_count <= nonzero_count_cutoff2)
  bac_exclude <- unique(c(rownames(bac_exclude_1), rownames(bac_exclude_2)))
  bac_include <- rownames(Ps_null_model)[!rownames(Ps_null_model) %in% bac_exclude]

  prevalence <- prevalence[match(bac_include, table = rownames(Ps_null_model))]
  absolute_sparsity <- absolute_sparsity[match(bac_include, table = rownames(Ps_null_model))]
  mean_abundance <- mean_abundance[match(bac_include, table = rownames(Ps_null_model))]

  Ps_null_model <- Ps_null_model %>%
    rownames_to_column() %>%
    dplyr::filter(rowname %in% bac_include) %>%
    column_to_rownames()

  Ps_poho_fdr <- Ps_poho_fdr %>%
    rownames_to_column() %>%
    dplyr::filter(rowname %in% bac_include) %>%
    column_to_rownames()

  delta <- delta %>%
    rownames_to_column() %>%
    dplyr::filter(rowname %in% bac_include) %>%
    column_to_rownames()

  variables <- bac_include

  return(list(prevalence, absolute_sparsity, mean_abundance,
              Ps_null_model, Ps_poho_fdr, delta, variables, bac_include))
}
