#'Generate result table as output in longdat_disc()
#' @param variable_col Internal function argument.
#' @param N Internal function argument.
#' @param Ps_conf_inv_model_unlist Internal function argument.
#' @param variables Internal function argument.
#' @param sel_fac Internal function argument.
#' @param Ps_conf_model_unlist Internal function argument.
#' @param model_q Internal function argument.
#' @param posthoc_q Internal function argument.
#' @param Ps_null_model_fdr Internal function argument.
#' @param Ps_null_model Internal function argument.
#' @param delta Internal function argument.
#' @param prevalence Internal function argument.
#' @param mean_abundance Internal function argument.
#' @param not_used Internal function argument.
#' @param Ps_effectsize Internal function argument.
#' @param data_type Internal function argument.
#' @param false_pos_count Internal function argument.
#' @param case_pairs Internal function argument.
#' @param case_pairs_name Internal function argument.
#' @param p_wilcox_final Internal function argument.
#' @param Ps_poho_fdr Internal function argument.
#' @import dplyr
#' @import stringr
#' @importFrom rlang .data
#' @importFrom stats as.formula confint cor.test kruskal.test
#'             na.omit p.adjust wilcox.test
#' @importFrom magrittr '%>%'
#' @name final_result_summarize_disc

final_result_summarize_disc <- function(variable_col, N,
                                        Ps_conf_inv_model_unlist, variables,
                                        sel_fac, Ps_conf_model_unlist,
                                   model_q, posthoc_q, Ps_null_model_fdr,
                                   Ps_null_model, delta, case_pairs,
                                   prevalence, mean_abundance, Ps_poho_fdr,
                                   not_used, Ps_effectsize, case_pairs_name,
                                   data_type, false_pos_count, p_wilcox_final){
  if (variable_col-1-2-length(not_used) > 0) {
    # There are potential confounders in raw input data
    # Generate potential confounding factors as output
    # Find the max number of how many confounders each bacteria has
    num <- c()
    for (i in 1:N) {
      num[i] <- sum(!is.na(Ps_conf_inv_model_unlist[i, ]))
    }

    prep_conf <- data.frame(matrix(NA, nrow = length(variables), ncol = 2))
    rownames(prep_conf) <-variables
    colnames(prep_conf) <- c("sel_fac_length", "Signal")
    for (i in seq_len(length(variables))) {
      prep_conf[i, 1] <- length(sel_fac[[i]])
      if (Ps_null_model_fdr[i, 1] >= model_q | is.na(Ps_null_model_fdr[i, 1]) |
          sum((Ps_poho_fdr[i, ]) < posthoc_q, na.rm = TRUE) == 0) {
        # Null time model q >= model_q is NS
        prep_conf[i, 2] <- "NS"
      } else {
        prep_conf[i, 2] <- "not_NS"
      }
    }

    confs <- prep_conf %>%
      rownames_to_column("Bacteria") %>%
      dplyr::filter(.data$sel_fac_length > 0 & .data$Signal == "not_NS")

    confs_num <- match(confs$Bacteria, variables)

    confound <- data.frame(matrix(NA, nrow = length(confs_num),
                                  ncol = 3*max(num + 1)))
    rownames(confound) <- confs$Bacteria
    colnames(confound) <- paste(rep(c("Confounder", "Confounding_type",
                                      "Effect_size"),
                                    time = max(num + 1)), sep = "",
                                rep(1:max(num + 1), each = 3))
    for (i in confs_num) {
      # loop through only the ones that do have confounding effect
      # (Strictly/Ambiguously/Completely confounded/deconfounded)
      for (j in seq_len(length(sel_fac[[i]]))) {
        # loop through different confounders
        tryCatch({
          c_name <- sel_fac[[i]][j]
          c_effectsize <- Ps_effectsize[i, c_name]
          c_type <- if (is.null(Ps_conf_model_unlist[i, sel_fac[[i]][j]])) {
            # if Ps_conf_model_unlist is null
            print("Strictly_deconfounded")
          } else { #Ps_conf_model_unlist isn't  null
            if (Ps_conf_model_unlist[i, sel_fac[[i]][j]] < posthoc_q &
                !is.na(Ps_conf_model_unlist[i, sel_fac[[i]][j]])) {
              # Confounding model p < posthoc_q
              print("Strictly_deconfounded")
            } else {# Confounding model p >= posthoc_q
              if (Ps_conf_inv_model_unlist[i, sel_fac[[i]][j]] < posthoc_q &
                  !is.na(Ps_conf_inv_model_unlist[i, sel_fac[[i]][j]])) {
                # Inverse confounding model p < posthoc_q
                print("Confounded")
              } else {# Inverse onfounding model p >= posthoc_q
                print("Ambiguously_deconfounded")
              }
            }}
          confound[match(i, confs_num), 3*j-2] <- c_name
          confound[match(i, confs_num), 3*j-1] <- c_type
          confound[match(i, confs_num), 3*j] <- c_effectsize
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
    confound <- confound %>%
      tibble::rownames_to_column("Feature")
    #write.table(x = confound, file = paste0(output_tag, "_confounders.txt"),
    # sep = "\t",
    #            row.names = TRUE, col.names = NA, quote = FALSE)

    # Create the result table
    result_table <- data.frame(matrix(nrow = length(variables),
                                      ncol = ncol(case_pairs)))
    row.names(result_table) <- variables
    colnames(result_table) <- paste("effect", sep = "_",
                                    unique(case_pairs_name))
    # The final determination of this feature's signal
    final_sig <- as.data.frame(matrix(nrow = N, ncol = 1, data = NA))
    row.names(final_sig) <- variables
    colnames(final_sig) <- "Final_sig"
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] >= model_q | is.na(Ps_null_model_fdr[i, 1])){
        # Null time model q >= model_q is NS
        final_sig[i, 1] <- "NS"
      } else { # Null time model q < model_q
        if (sum((Ps_poho_fdr[i, ]) < posthoc_q, na.rm = TRUE) == 0){
          # All post-hoc q >= posthoc_q
          final_sig[i, 1] <- "NS"
        } else {# At least one post-hoc q < posthoc_q
          if (length(sel_fac[[i]]) == 0) {
            # No sel_fac, meaning no confounding effect
            if (Ps_null_model[i, 2] == "Good" & !is.na(Ps_null_model[i, 2])) {
              # Proper confidence interval: OK and no confounder
              final_sig[i, 1] <- "OK_nc"
            } else {# Improper confidence interval: OK but doubtful
              final_sig[i, 1] <- "OK_d"
            }
          } else {# There are sel_fac, so decide the signal
                  # based on confound table
            subconfound <- as.data.frame(
              confound[str_which(string = confound$Feature,
                                 pattern = as.character(variables[i])), ])
            subconfound_columns <-
              subconfound[ , str_which(string = colnames(confound),
                                       pattern = "Confounding_type")]
            if (sum(stringr::str_detect(string = subconfound_columns,
                                        pattern = "Confounded"),
                    na.rm = TRUE) > 0) {
              # There is "confounding" signals: Confounded
              final_sig[i, 1] <- "C"
            } else if(sum(stringr::str_detect(string = subconfound_columns,
                                              pattern = "Confounded"),
                          na.rm = TRUE) == 0){
              # There isn't "confounding" signals
              if (sum(
                stringr::str_detect(string = subconfound_columns,
                                    pattern = "Ambiguously_deconfounded"),
                na.rm = TRUE) > 0) {
                # If there is "ambiguously deconfounded" signals:
                # Ambiguously deconfounded
                final_sig[i, 1] <- "AD"
              } else { # There isn't "ambiguously deconfounded" signals:
                # OK and strictly deconfounded
                final_sig[i, 1] <- "OK_sd"
              }
            }
          }
        }
      }
    }

    # The effect type determination
    effect_type <- as.data.frame(matrix(nrow = N, ncol = ncol(case_pairs),
                                        data = NA))
    row.names(effect_type) <- variables
    colnames(effect_type) <- unique(case_pairs_name)
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] < model_q & !is.na(Ps_null_model_fdr[i, 1])){
        # Null time model q < model_q and isn't NA
        for (j in seq_len(ncol(Ps_poho_fdr))) {
          if (Ps_poho_fdr[i, j] < posthoc_q & !is.na(Ps_poho_fdr[i, j])){
            # Post-hoc q < posthoc_q and isn't NA
            if (delta[i, j] > 0 & !is.na(delta[i, j])) {
              # Delta > 0 and isn't NA
              effect_type[i, j] <- "Enriched"
            } else if (delta[i, j] < 0 & !is.na(delta[i, j])) {
              # Delta < 0 and isn't NA
              effect_type[i, j] <- "Decreased"
            }
          } else {
            effect_type[i, j] <- "NS"
          }
        }
      } else {
        effect_type[i, ] <- "NS"
      }
    }

    # Bind the columns together
    result_table <- cbind(prevalence, mean_abundance, final_sig, effect_type,
                          delta, Ps_null_model_fdr, Ps_poho_fdr)
    colnames(result_table) <- c("Prevalence_percentage",
                                "Mean_abundance", "Signal",
                                paste("Effect_", sep = "",
                                      unique(case_pairs_name)),
                                paste("EffectSize_", sep = "",
                                      unique(case_pairs_name)),
                                "Null_time_model_q",
                                paste("Post-hoc_q_", sep = "",
                                      unique(case_pairs_name)))
    if (data_type == "count") {
      if (false_pos_count > 0) {
        result_table <- cbind(result_table, p_wilcox_final)
      }
    }
  } else if (variable_col-1-2-length(not_used) == 0) {
    # No potential confounders in raw input data
    result_table <- data.frame(matrix(nrow = length(variables),
                                      ncol = ncol(case_pairs)))
    row.names(result_table) <- variables

    # The final determination of this feature's signal
    final_sig <- as.data.frame(matrix(nrow = N, ncol = 1, data = NA))
    row.names(final_sig) <- variables
    colnames(final_sig) <- "Final_sig"
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] >= model_q| is.na(Ps_null_model_fdr[i, 1])) {
        # Null time model q >= model_q is NS
        final_sig[i, 1] <- "NS"
      } else { # Null time model q < model_q
        if (sum((Ps_poho_fdr[i, ]) < posthoc_q, na.rm = TRUE) == 0) {
          # All post-hoc q >= posthoc_q
          final_sig[i, 1] <- "NS"
        } else {# At least one post-hoc q < posthoc_q
          if (Ps_null_model[i, 2] == "Good" & !is.na(Ps_null_model[i, 2])) {
            # Proper confidence interval: OK and no confounder
            final_sig[i, 1] <- "OK_nc"
          } else {# Improper confidence interval: OK but doubtful
            final_sig[i, 1] <- "OK_d"
          }
        }
      }
    }

    # The effect type determination
    effect_type <- as.data.frame(matrix(nrow = N, ncol = ncol(case_pairs),
                                        data = NA))
    row.names(effect_type) <- variables
    colnames(effect_type) <- unique(case_pairs_name)
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] < model_q & !is.na(Ps_null_model_fdr[i, 1])){
        # Null time model q < model_q and isn't NA
        for (j in seq_len(ncol(Ps_poho_fdr))) {
          if (Ps_poho_fdr[i, j] < posthoc_q & !is.na(Ps_poho_fdr[i, j])) {
            # Post-hoc q < posthoc_q and isn't NA
            if (delta[i, j] > 0 & !is.na(delta[i, j])) {
              # Delta > 0 and isn't NA
              effect_type[i, j] <- "Enriched"
            } else if (delta[i, j] < 0 & !is.na(delta[i, j])) {
              # Delta < 0 and isn't NA
              effect_type[i, j] <- "Decreased"
            }
          } else {
            effect_type[i, j] <- "NS"
          }
        }
      } else {
        effect_type[i, ] <- "NS"
      }
    }

    # Bind the columns together
    result_table <- cbind(prevalence, mean_abundance, final_sig, effect_type,
                          delta, Ps_null_model_fdr, Ps_poho_fdr)
    colnames(result_table) <-
      c("Prevalence_percentage", "Mean_abundance",
        "Signal", paste("Effect_", sep = "",
                        unique(case_pairs_name)),
        paste("EffectSize_", sep = "",
              unique(case_pairs_name)),
        "Null_time_model_q",
        paste("Post-hoc_q_", sep = "", unique(case_pairs_name)))
    if (data_type == "count") {
      if (false_pos_count > 0) {
        result_table <- cbind(result_table, p_wilcox_final)
      }
    }
  }
  result_table <- result_table %>%
    rownames_to_column("Feature")
  #write.table(x = result_table, file = paste0(output_tag,
  # "_result_table.txt"), sep = "\t",
  #            row.names = TRUE, col.names = NA, quote = FALSE)
  if (variable_col-1-2-length(not_used) > 0) {
    return(list(confound, result_table))
  } else if (variable_col-1-2-length(not_used) == 0) {
    return(result_table)
  }
}
