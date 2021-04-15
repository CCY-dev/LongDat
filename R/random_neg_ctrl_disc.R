#' Randomized negative control for count data in longdat_disc()
#' @param test_var Internal function argument.
#' @param variable_col Internal function argument.
#' @param fac_var Internal function argument.
#' @param not_used Internal function argument.
#' @param factors Internal function argument.
#' @param data Internal function argument.
#' @param N Internal function argument.
#' @param data_type Internal function argument.
#' @param variables Internal function argument.
#' @param case_pairs Internal function argument.
#' @param adjustMethod Internal function argument.
#' @param model_q Internal function argument.
#' @param posthoc_q Internal function argument.
#' @param theta_cutoff Internal function argument.
#' @param nonzero_count_cutoff1 Internal function argument.
#' @param nonzero_count_cutoff2 Internal function argument.
#' @param output_tag Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data

random_neg_ctrl_disc <- function(test_var, variable_col, fac_var, not_used, factors, data, N, data_type, variables, case_pairs,
                            adjustMethod, model_q, posthoc_q, theta_cutoff, nonzero_count_cutoff1, nonzero_count_cutoff2,
                            output_tag, verbose) {
  ######### Randomize the raw data first
  # Shuffle the rows randomly
  value_for_random <- data[ , variable_col:ncol(data)]
  value_random <- value_for_random[sample(nrow(value_for_random)), ]
  data_randomized <- as.data.frame(cbind(data[ , 1:(variable_col-1)], value_random))

  # Remove all the symbols from variables
  colnames(data_randomized)[variable_col:ncol(data_randomized)] <- fix_name_fun(colnames(data_randomized)[variable_col:ncol(data_randomized)])

  # Melt randomized data
  predictor_names <- (colnames(data_randomized))[1: (variable_col - 1)]
  melt_data_random <- reshape2::melt (data_randomized, id = predictor_names)
  # Omit the rows whose value column equals to NA
  melt_data_random <- melt_data_random %>% tidyr::drop_na(.data$value)
  # Remove all dots in the bacteria name or it will cause problem
  melt_data_random$variable <- gsub(".", "_", melt_data_random$variable, fixed = TRUE)

  # Make sure that all the columns are in the right class
  # Columns mentioned in fac_var, and the second last column in melt data are factors
  # Columns not in fac_var, and the last column in melt data are numerical numbers
  "%notin%" <- Negate("%in%")
  num_var <- c(which(1:(ncol(melt_data_random)-2) %notin% fac_var), ncol(melt_data_random))
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
  p_poho_neg_crtl <- data.frame(matrix(nrow = length(variables), ncol = ncol(case_pairs)))

  suppressWarnings(
    for (i in 1:N) { # loop through all variables
      aVariable = variables[i]
      if (verbose == T) {print(i)}
      subdata_random <- subset(melt_data_random, melt_data_random$variable == as.character(aVariable))
      tryCatch({
        # Negative binomial
        fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
        m3 <- glmmTMB(formula = fmla2, data = subdata_random, family = nbinom2, na.action = na.omit, REML = F)

        # Extract dispersion theta out of model
        Theta_random[i] <- sigma(m3)

        # Wald Chisq test
        Ps_neg_ctrl[i, 1] <- car::Anova(m3, type=c("II"),  test.statistic=c("Chisq"))$"Pr(>Chisq)"

        # Post-hoc test
        m_means_neg_ctrl <- emmeans(object = m3, specs = test_var)
        e <- as.data.frame(pairs(m_means_neg_ctrl, adjust = "fdr", infer = c(F, T), reverse = F))
        for (m in 1:ncol(case_pairs)) {
          p_poho_neg_crtl[i, m] <- e$p.value[m]
        }

        # Calculate confidence interval for test_var
        # Followed by the determination of the signs. For "- -" or "+ +", it means that the doesn't span 0.
        # But "- +" means that the CI spans 0. Here I sum the signs over the row, sum != 0 means it's OK.
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
    rownames_to_column() %>%
    mutate(NB_theta = Theta_random) %>%
    column_to_rownames()
  row.names(p_poho_neg_crtl) <- variables
  case_pairs_name <- c()
  for (i in 1:ncol(case_pairs)) {
    cpn <- paste0(case_pairs[1, i], "_", case_pairs[2, i])
    case_pairs_name <- c(case_pairs_name, cpn)
  }
  colnames(p_poho_neg_crtl) <-  paste0("p_", case_pairs_name)

  ####### Effect size
  case_pairs_random <- combn(x = sort(unique(melt_data_random[ , test_var])), m = 2)
  delta_random <- data.frame(matrix(nrow = length(row.names(Ps_neg_ctrl)), ncol = ncol(case_pairs_random)))
  case_pairs_random_name <- c()
  library(orddom)
  for (i in 1:length(row.names(Ps_neg_ctrl))) { # loop through all variables
    if (verbose == T) {print(i)}
    cVariable = variables[i]
    subdata_pre_random <- subset(melt_data_random, variable == as.character(cVariable))
    counts_random <- subdata_pre_random %>% dplyr::count(.data$Individual)
    # Exclude the ones not having data points at ALL timepoints
    exclude_random <- counts_random$Individual[which(counts_random$n != length(unique(data_randomized[ , test_var])))]
    if (length(exclude_random) > 0) {
      subdata2_random <- subset(subdata_pre_random, !Individual %in% exclude_random)
    } else {
      subdata2_random <- subdata_pre_random
    }
    for (k in 1:ncol(case_pairs_random)) { # loop through each case pair
      sub5 <- subdata2_random[subdata2_random[ , test_var] == case_pairs_random[1,k], ]
      sub6 <- subdata2_random[subdata2_random[ , test_var] == case_pairs_random[2,k], ]
      d_random <- as.numeric(orddom::dmes(x = sub5$value, y = sub6$value)$dw)
      delta_random[i, k] <- d_random
      name_random <- paste(case_pairs_random[1,k], sep = "_", case_pairs_random[2,k])
      case_pairs_random_name <- c(case_pairs_random_name, name_random)
    }
  }
  row.names(delta_random) <- row.names(Ps_neg_ctrl)
  colnames(delta_random) <- paste("effect_size", sep = "_", unique(case_pairs_random_name))

  #### Remove high theta and low prevalence ones from randomized result
  absolute_sparsity_random <- c()
  for (i in 1:ncol(value_random)) {
    absolute_sparsity_random[i] <- sum(value_random[ , i] == 0, na.rm = T)
  }

  non_zero_count_randomized <- nrow(data_randomized) - absolute_sparsity_random
  Ps_neg_ctrl <- Ps_neg_ctrl %>%
    rownames_to_column() %>%
    mutate(non_zero_count = non_zero_count_randomized) %>%
    column_to_rownames()

  bac_exclude_1_random <- subset(Ps_neg_ctrl, non_zero_count <= nonzero_count_cutoff1 & NB_theta >= theta_cutoff)
  bac_exclude_2_random <- subset(Ps_neg_ctrl, non_zero_count <= nonzero_count_cutoff2)
  bac_exclude_random <- unique(c(rownames(bac_exclude_1_random), rownames(bac_exclude_2_random)))
  bac_include_random <- rownames(Ps_neg_ctrl)[!rownames(Ps_neg_ctrl) %in% bac_exclude_random]

  Ps_neg_ctrl_filterd <- Ps_neg_ctrl %>%
    rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include_random) %>%
    column_to_rownames()

  delta_random_filtered <- delta_random %>%
    rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include_random) %>%
    column_to_rownames()

  p_poho_neg_crtl_filtered <- p_poho_neg_crtl %>%
    rownames_to_column() %>%
    dplyr::filter(.data$rowname %in% bac_include_random) %>%
    column_to_rownames()

  ####### FDR correction
  adjust_fun <- function(x) p.adjust(p = x, method = adjustMethod)
  Ps_neg_ctrl_fdr <- adjust_fun(Ps_neg_ctrl_filterd$Neg_ctrl_model_p)
  Ps_neg_ctrl_filterd <- Ps_neg_ctrl_filterd %>%
    rownames_to_column() %>%
    mutate(P_fdr = Ps_neg_ctrl_fdr) %>%
    column_to_rownames()

  ####### Write randomized control table
    signal_neg_ctrl_tbl <- data.frame(matrix(nrow = length(row.names(Ps_neg_ctrl_filterd)), ncol = ncol(case_pairs_random),
                                             data = NA, byrow = F))
    colnames(signal_neg_ctrl_tbl) <- paste0("Signal_", case_pairs_name)
    rownames(signal_neg_ctrl_tbl) <- rownames(Ps_neg_ctrl_filterd)
    for (i in 1:nrow(p_poho_neg_crtl_filtered)) {
      for (j in 1:ncol(case_pairs_random)) {
        signal_neg_ctrl_tbl[i, j] <- ifelse(Ps_neg_ctrl_filterd$P_fdr[i] < model_q & p_poho_neg_crtl_filtered[i, j] < posthoc_q & !is.na(Ps_neg_ctrl_filterd$P_fdr[i]) & !is.na(p_poho_neg_crtl_filtered[i, j]),
                                     yes = "False_positive", no = "Negative")
      }}
    Summary_signal <- c()
    for (i in 1:nrow(signal_neg_ctrl_tbl)) {
      Summary_signal[i] <- ifelse(sum(str_detect(signal_neg_ctrl_tbl[i, ], "False_positive")) > 0,
                                  yes = "False_positive", no = "Negative")
    }
    result_neg_ctrl <- cbind(Ps_neg_ctrl_filterd[ ,c(2, 5)], Summary_signal, signal_neg_ctrl_tbl,
                             p_poho_neg_crtl_filtered, delta_random_filtered)
    colnames(result_neg_ctrl) <- c("Signal_of_CI_signs", "Model_q", "Final_signal", paste0("Signal_", case_pairs_name),
                                   paste0("Posthoc_q_", case_pairs_name), paste0("Effect_size_", case_pairs_name))
    # Change the rownames to avoid confusion for the users
    rownames(result_neg_ctrl) <- paste0("Randomized_feature_", 1:nrow(result_neg_ctrl))

     result_neg_ctrl_sig <- result_neg_ctrl %>%
      filter(.data$Final_signal == "False_positive" & .data$Signal_of_CI_signs == "Good") %>%
      dplyr::select(-1)
     if (nrow(result_neg_ctrl_sig) > 0) {
       write.table(x = result_neg_ctrl_sig, file = paste0(output_tag, "_randomized_control.txt"), sep = "\t",
                   row.names = T, col.names = NA, quote = F)
     }
     return(result_neg_ctrl)
  }


