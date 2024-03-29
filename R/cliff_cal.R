#' Effect size (Cliff's delta) calculation in longdat_disc() pipeline
#' @param melt_data Internal function argument.
#' @param Ps_poho_fdr Internal function argument.
#' @param variables Internal function argument.
#' @param test_var Internal function argument.
#' @param data Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data
#' @importFrom stats as.formula confint cor.test kruskal.test na.omit
#'             p.adjust wilcox.test
#' @import dplyr
#' @importFrom magrittr '%>%'
#' @name cliff_cal
utils::globalVariables(c("variable", "Individual"))

cliff_cal <- function(melt_data, Ps_poho_fdr, variables, test_var,
                      data, verbose) {
  case_pairs <- combn(sort(unique(melt_data[ , test_var])), m = 2)
  delta <- data.frame(matrix(nrow = length(row.names(Ps_poho_fdr)),
                             ncol = ncol(case_pairs)))
  case_pairs_name <- c()
  for (i in seq_len(length(row.names(Ps_poho_fdr)))) {
    # loop through all variables
    if (verbose == TRUE) {print(i)}
    bVariable <- variables[i]
    subdata_pre <- subset(melt_data, variable == bVariable)
    counts <- subdata_pre %>% dplyr::count(.data$Individual)
    # Exclude the ones not having data points at ALL timepoints
    exclude <-
      counts$Individual[which(counts$n != length(unique(data[ , test_var])))]
    if (length(exclude) > 0) {
      subdata2 <- subset(subdata_pre, !Individual %in% exclude)
    } else {
      subdata2 <- subdata_pre
    }
    for (k in seq_len(ncol(case_pairs))) { # loop through each case pair
      sub3 <- subdata2[subdata2[ , test_var] == case_pairs[1,k], ] %>%
        dplyr::arrange(Individual)
      sub4 <- subdata2[subdata2[ , test_var] == case_pairs[2,k], ] %>%
        dplyr::arrange(Individual)
      #d <- as.numeric(dmes(x = sub3$value, y = sub4$value)$dw)
      #delta[i, k] <- d
      ### Orddom is outdated, so calculate Cliff's delta as below:
      #Note that the order is treatment - control
      delta[i, k] <- mean(sign(sub4$value-sub3$value), na.rm = TRUE)
      name <- paste(case_pairs[1,k], sep = "_", case_pairs[2,k])
      case_pairs_name <- c(case_pairs_name, name)
    }
  }
  row.names(delta) <- row.names(Ps_poho_fdr)
  colnames(delta) <- paste("effect_size", sep = "_", unique(case_pairs_name))
  return(list(case_pairs, delta))
}
