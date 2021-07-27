#' Randomized negative control for count data in longdat_cont()
#' @param test_var Internal function argument.
#' @param variable_col Internal function argument.
#' @param fac_var Internal function argument.
#' @param not_used Internal function argument.
#' @param factors Internal function argument.
#' @param data Internal function argument.
#' @param N Internal function argument.
#' @param data_type Internal function argument.
#' @param variables Internal function argument.
#' @param adjustMethod Internal function argument.
#' @param model_q Internal function argument.
#' @param posthoc_q Internal function argument.
#' @param theta_cutoff Internal function argument.
#' @param nonzero_count_cutoff1 Internal function argument.
#' @param nonzero_count_cutoff2 Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data
#' @importFrom stats as.formula confint cor.test kruskal.test
#'             na.omit p.adjust wilcox.test
#' @importFrom car Anova
#' @importFrom magrittr '%>%'
#' @import reshape2
#' @import tidyr
#' @import dplyr
#' @import glmmTMB
#' @import tibble
#' @import dplyr
#' @name random_neg_ctrl_cont
utils::globalVariables(c("non_zero_count", "NB_theta"))

random_neg_ctrl_cont <- function(test_var, variable_col, fac_var, not_used,
                                 factors, data, N, data_type, variables,
                                 adjustMethod, model_q, posthoc_q,
                                 theta_cutoff, nonzero_count_cutoff1,
                                 nonzero_count_cutoff2, verbose) {
  ######### Randomize the raw data first
  # Shuffle the rows randomly
  value_for_random <- data[ , variable_col:ncol(data)]
  value_random <- value_for_random[sample(nrow(value_for_random)), ]
  data_randomized <- as.data.frame(cbind(data[ , 1:(variable_col-1)],
                                         value_random))

  # Remove all the symbols from variables
  colnames(data_randomized)[variable_col:ncol(data_randomized)] <-
    fix_name_fun(colnames(data_randomized)[variable_col:ncol(data_randomized)])

  # Melt randomized data
  predictor_names <- (colnames(data_randomized))[1: (variable_col - 1)]
  melt_data_random <- reshape2::melt (data_randomized, id = predictor_names)
  # Omit the rows whose value column equals to NA
  melt_data_random <- melt_data_random %>% tidyr::drop_na(.data$value)
  # Remove all dots in the bacteria name or it will cause problem
  melt_data_random$variable <- gsub(".", "_", melt_data_random$variable,
                                    fixed = TRUE)

  # Make sure that all the columns are in the right class
  # Columns mentioned in fac_var,
  # and the second last column in melt data are factors
  # Columns not in fac_var, and the last column in
  # melt data are numerical numbers
  "%notin%" <- Negate("%in%")
  num_var <- c(which(1:(ncol(melt_data_random)-2) %notin% fac_var),
               ncol(melt_data_random))
  fac_var <- c(fac_var, ncol(melt_data_random)-1)
  for (i in fac_var) {
    melt_data_random[ ,i] <- as.factor(melt_data_random[ ,i])
  }
  for (i in num_var) {
    melt_data_random[ ,i] <- as.numeric(as.character(melt_data_random[ ,i]))
  }

  # Remove the not-used columns in melt_data
  if (!is.null(not_used)) {
    melt_data_random <- melt_data_random %>% dplyr::select(-c(not_used))
  }
  # Change the first column name of melt_data to "Individual"
  colnames(melt_data_random)[1] <- "Individual"

  ######### Then do the model test
  Ps_neg_ctrl <- as.data.frame(matrix(data = NA, nrow = N, ncol = 2))
  rownames(Ps_neg_ctrl) <- gsub(".", "_", variables, fixed = TRUE)
  colnames(Ps_neg_ctrl) <- c("Neg_ctrl_model_p", "Signal_of_CI_signs")
  Theta_random <- c()

  suppressWarnings(
    for (i in 1:N) { # loop through all variables
      aVariable = variables[i]
      if (verbose == T) {print(i)}
      subdata_random <- subset(melt_data_random, variable == aVariable)
      tryCatch({
        # Negative binomial
        fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
        m3 <- glmmTMB::glmmTMB(formula = fmla2, data = subdata_random,
                               family = nbinom2, na.action = na.omit, REML = F)

        # Extract dispersion theta out of model
        Theta_random[i] <- glmmTMB::sigma(m3)

        # Wald Chisq test
        Ps_neg_ctrl[i, 1] <-
          car::Anova(m3, type=c("II"),  test.statistic=c("Chisq"))$"Pr(>Chisq)"

        # Calculate confidence interval for test_var
        # Followed by the determination of the signs. For "- -" or "+ +",
        # it means that the doesn't span 0.
        # But "- +" means that the CI spans 0.
        # Here I sum the signs over the row, sum != 0 means it's OK.
        remove(ci)
        ci <- as.data.frame(confint(m3))
        ci <- ci[str_detect(row.names(ci), test_var), 1:2]
        if (nrow(ci) > 1) {
          if (any(apply(sign(ci), 1, sum, na.rm = T) != 0)) {
            Ps_neg_ctrl[i, 2] <- "Good"
          } else {
            Ps_neg_ctrl[i, 2] <- "Bad"
          }
        } else if (nrow(ci) == 1) {
          if (sign(ci)[1] == sign(ci)[2]) {
            Ps_neg_ctrl[i, 2] <- "Good"
          } else {
            Ps_neg_ctrl[i, 2] <- "Bad"
          }
        }
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })
  Ps_neg_ctrl <- Ps_neg_ctrl %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(NB_theta = Theta_random) %>%
    tibble::column_to_rownames()

  ######### Post-hoc test and effect size
  # Here uses Spearman's correlation
  p_poho_random <- as.data.frame(matrix(nrow = length(variables), ncol = 1))
  assoc_random <- as.data.frame(matrix(nrow = length(variables), ncol = 1))

  for (i in 1:N) { # loop through all variables
    if (verbose == T) {print(i)}
    cVariable = variables[i]
    subdata_random <- subset(melt_data_random, variable == cVariable)
    # Here set the "test_var" to numeric
    subdata_random[ , test_var] <- as.numeric(subdata_random[ , test_var])
    e <- cor.test(subdata_random[ , test_var], subdata_random$value,
                  method = "spearman")
    p_c_random <- e$p.value
    a_c_random <- e$estimate
    p_poho_random[i, 1] <- p_c_random
    assoc_random[i, 1] <- a_c_random
  }
  row.names(p_poho_random) <- variables
  row.names(assoc_random) <- variables
  colnames(p_poho_random) <- "p_post-hoc"
  colnames(assoc_random) <- "association"

  #### Remove high theta and low prevalence ones from randomized result
  absolute_sparsity_random <- c()
  for (i in 1:ncol(value_random)) {
    absolute_sparsity_random[i] <- sum(value_random[ , i] == 0, na.rm = T)
  }

  non_zero_count_randomized <- nrow(data_randomized) - absolute_sparsity_random
  Ps_neg_ctrl <- Ps_neg_ctrl %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(non_zero_count = non_zero_count_randomized) %>%
    tibble::column_to_rownames()

  bac_exclude_1_random <-
    subset(Ps_neg_ctrl, non_zero_count <= nonzero_count_cutoff1 &
             NB_theta >= theta_cutoff)
  bac_exclude_2_random <-
    subset(Ps_neg_ctrl, non_zero_count <= nonzero_count_cutoff2)
  bac_exclude_random <- unique(c(rownames(bac_exclude_1_random),
                                 rownames(bac_exclude_2_random)))
  bac_include_random <-
    rownames(Ps_neg_ctrl)[!rownames(Ps_neg_ctrl) %in% bac_exclude_random]

  Ps_neg_ctrl_filterd <- Ps_neg_ctrl %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include_random) %>%
    tibble::column_to_rownames()

  p_poho_random_filtered <- p_poho_random %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include_random) %>%
    tibble::column_to_rownames()

  assoc_random_filtered <- assoc_random %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include_random) %>%
    tibble::column_to_rownames()

  ####### FDR correction
  adjust_fun <- function(x) p.adjust(p = x, method = adjustMethod)
  Ps_neg_ctrl_fdr <- adjust_fun(Ps_neg_ctrl_filterd$Neg_ctrl_model_p)
  Ps_neg_ctrl_filterd <- Ps_neg_ctrl_filterd %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(P_fdr = Ps_neg_ctrl_fdr) %>%
    tibble::column_to_rownames()

  ####### Write randomized control table
    signal_neg_ctrl_tbl <-
      data.frame(matrix(nrow = length(row.names(Ps_neg_ctrl_filterd)),
                        ncol = 1, data = NA, byrow = F))
    colnames(signal_neg_ctrl_tbl) <- "Signal"
    rownames(signal_neg_ctrl_tbl) <- rownames(Ps_neg_ctrl_filterd)
    for (i in 1:nrow(Ps_neg_ctrl_filterd)) {
        signal_neg_ctrl_tbl[i, 1] <-
          ifelse(Ps_neg_ctrl_filterd$P_fdr[i] < model_q &
                   p_poho_random_filtered$`p_post-hoc`[i] < posthoc_q &
                   !is.na(Ps_neg_ctrl_filterd$P_fdr[i]) &
                   !is.na(p_poho_random_filtered$`p_post-hoc`[i]),
                 yes = "False_positive", no = "Negative")
      }
    result_neg_ctrl <- cbind(Ps_neg_ctrl_filterd[ ,c(2, 5)], signal_neg_ctrl_tbl,
                             p_poho_random_filtered$`p_post-hoc`,
                             assoc_random_filtered)
    colnames(result_neg_ctrl) <- c("Signal_of_CI_signs", "Model_q", "Signal",
                                   "Posthoc_q", "Effect_size")

    # Change the rownames to avoid confusion for the users
    rownames(result_neg_ctrl) <- paste0("Randomized_feature_",
                                        1:nrow(result_neg_ctrl))

    result_neg_ctrl_sig <- result_neg_ctrl %>%
      dplyr::filter(.data$Signal == "False_positive" &
                      .data$Signal_of_CI_signs == "Good") %>%
      dplyr::select(-1)
    false_pos_count <- nrow(result_neg_ctrl_sig)

    if (nrow(result_neg_ctrl_sig) > 0) {
    #write.table(x = result_neg_ctrl_sig, file =
      # paste0("_randomized_control.txt"), sep = "\t",
    #            row.names = T, col.names = NA, quote = F)
    }
    return(list(result_neg_ctrl, false_pos_count))
  }


