#' Effect size (Cliff's delta) calculation in longdat_disc() pipeline
#' @param melt_data Internal function argument.
#' @param Ps_poho_fdr Internal function argument.
#' @param variables Internal function argument.
#' @param test_var Internal function argument.
#' @param data Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data

cliff_cal <- function(melt_data, Ps_poho_fdr, variables, test_var, data, verbose) {
  case_pairs <- combn(sort(unique(melt_data[ , test_var])), m = 2)
  delta <- data.frame(matrix(nrow = length(row.names(Ps_poho_fdr)), ncol = ncol(case_pairs)))
  case_pairs_name <- c()
  library(orddom)
  for (i in 1:length(row.names(Ps_poho_fdr))) { # loop through all variables
    if (verbose == T) {print(i)}
    bVariable = variables[i]
    subdata_pre <- subset(melt_data, variable == bVariable)
    counts <- subdata_pre %>% dplyr::count(.data$Individual)
    # Exclude the ones not having data points at ALL timepoints
    exclude <- counts$Individual[which(counts$n != length(unique(data[ , test_var])))]
    if (length(exclude) > 0) {
      subdata2 <- subset(subdata_pre, !Individual %in% exclude)
    } else {
      subdata2 <- subdata_pre
    }
    for (k in 1:ncol(case_pairs)) { # loop through each case pair
      sub3 <- subdata2[subdata2[ , test_var] == case_pairs[1,k], ]
      sub4 <- subdata2[subdata2[ , test_var] == case_pairs[2,k], ]
      d <- as.numeric(orddom::dmes(x = sub3$value, y = sub4$value)$dw)
      delta[i, k] <- d
      name <- paste(case_pairs[1,k], sep = "_", case_pairs[2,k])
      case_pairs_name <- c(case_pairs_name, name)
    }
  }
  row.names(delta) <- row.names(Ps_poho_fdr)
  colnames(delta) <- paste("effect_size", sep = "_", unique(case_pairs_name))
  return(list(case_pairs, delta))
}
