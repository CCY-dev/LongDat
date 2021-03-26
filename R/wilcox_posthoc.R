#' Wilcoxon post-hoc test
#' @param result_neg_ctrl Internal function argument.
#' @param model_q Internal function argument.
#' @param melt_data Internal function argument.
#' @param test_var Internal function argument.
#' @param variables Internal function argument.
#' @param data Internal function argument.
#' @param N Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data

wilcox_posthoc <- function(result_neg_ctrl, model_q, melt_data, test_var, variables, data, N, verbose) {
  #Count false positives
  false_pos <- result_neg_ctrl %>%
    filter(.data$Final_signal == "False_positive" & .data$Signal_of_CI_signs == "Good")
  false_pos_count <- nrow(false_pos)

  #Do it if there are false positive in the randomized result
  if (false_pos_count > 0) {
    case_pairs <- combn(x = sort(unique(melt_data[ , test_var])), m = 2)
    p_wilcox <- data.frame(matrix(nrow = length(variables), ncol = ncol(case_pairs)))
    case_pairs_name <- c()
    for (i in 1:N) { # loop through all variables
      if (verbose == T) {print(i)}
      bVariable = variables[i]
      subdata_pre <- subset(melt_data, variable == bVariable)
      counts <- subdata_pre %>% dplyr::count(.data$Individual)
      # Exclude the samples that don't have value at all time points
      exclude <- counts$Individual[which(counts$n != length(unique(data[ , test_var])))]
      if (length(exclude) > 0) {
        subdata2 <- subset(subdata_pre, !Individual %in% exclude)
      } else {
        subdata2 <- subdata_pre
      }
      for (k in 1:ncol(case_pairs)) { # loop through each case pair
        sub3 <- subdata2[subdata2[ , test_var] == case_pairs[1,k], ]
        sub4 <- subdata2[subdata2[ , test_var] == case_pairs[2,k], ]
        # Here use "paired wilcoxon test because it's longitudinal data
        suppressWarnings(
          p_w <- stats::wilcox.test(sub3$value, sub4$value, paired = T)$p.value)
        p_wilcox[i, k] <- p_w
        name <- paste(case_pairs[1,k], sep = "_", case_pairs[2,k])
        case_pairs_name <- c(case_pairs_name, name)
      }
    }
    row.names(p_wilcox) <- variables
    colnames(p_wilcox) <- unique(paste("p_", case_pairs_name))
  } else {
    p_wilcox <- data.frame(matrix(NA, nrow = 2, ncol = 2))
  }
  return(list(false_pos_count, p_wilcox))
}
