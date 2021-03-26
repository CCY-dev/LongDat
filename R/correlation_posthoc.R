#' Post-hoc test based on correlation test for longdat_cont().
#' @param N Internal function argument.
#' @param variables Internal function argument.
#' @param melt_data Internal function argument.
#' @param test_var Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data

correlation_posthoc <- function(variables, verbose, melt_data, test_var, N) {
  # Here uses Spearman's correlation
  p_poho <- as.data.frame(matrix(nrow = length(variables), ncol = 1))
  assoc <- as.data.frame(matrix(nrow = length(variables), ncol = 1))

  for (i in 1:N) { # loop through all variables
    if (verbose == T) {print(i)}
    bVariable = variables[i]
    subdata <- subset(melt_data, variable == bVariable)
    # Here set the "test_var" to numeric
    subdata[ , test_var] <- as.numeric(subdata[ , test_var])
    c <- stats::cor.test(subdata[ , test_var], subdata$value, method = "spearman")
    p_c <- c$p.value
    a_c <- c$estimate
    p_poho[i, 1] <- p_c
    assoc[i, 1] <- a_c
  }
  row.names(p_poho) <- variables
  row.names(assoc) <- variables
  colnames(p_poho) <- "p_post-hoc"
  colnames(assoc) <- "association"
  return(list(p_poho, assoc))
}
