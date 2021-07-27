#'Generate result table as output in longdat_cont()
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
#' @param assoc Internal function argument.
#' @param prevalence Internal function argument.
#' @param mean_abundance Internal function argument.
#' @param p_poho Internal function argument.
#' @param not_used Internal function argument.
#' @param Ps_effectsize Internal function argument.
#' @param data_type Internal function argument.
#' @param false_pos_count Internal function argument.
#' @importFrom stats as.formula confint cor.test kruskal.test na.omit p.adjust wilcox.test
#' @importFrom magrittr '%>%'
#' @importFrom rlang .data
#' @import tibble
#' @import dplyr
#' @name final_result_summarize_cont

final_result_summarize_cont <- function(variable_col, N, Ps_conf_inv_model_unlist, variables, sel_fac, Ps_conf_model_unlist,
                                        model_q, posthoc_q, Ps_null_model_fdr, Ps_null_model, assoc, prevalence,
                                        mean_abundance, p_poho, not_used, Ps_effectsize, data_type, false_pos_count) {
  if (variable_col-1-2-length(not_used) > 0) {# There are potential confounders in raw input data
    # Generate potential confounding factors as output
    # Find the max number of how many confounders each bacteria has
    num <- c()
    for (i in 1:N) {
      num[i] <- sum(!is.na(Ps_conf_inv_model_unlist[i, ]))
    }

    prep_conf <- data.frame(matrix(NA, nrow = length(variables), ncol = 2))
    rownames(prep_conf) <-variables
    colnames(prep_conf) <- c("sel_fac_length", "Signal")
    for (i in 1:length(variables)) {
      prep_conf[i, 1] <- length(sel_fac[[i]])
      if (Ps_null_model_fdr[i, 1] >= model_q | is.na(Ps_null_model_fdr[i, 1]) | (p_poho[i, 1]) >= posthoc_q | is.na(p_poho[i, 1])) {# Null time model q >= model_q is NS
        prep_conf[i, 2] <- "NS"
      } else {
        prep_conf[i, 2] <- "not_NS"
      }
    }

    confs <- prep_conf %>%
      tibble::rownames_to_column("Bacteria") %>%
      dplyr::filter(.data$sel_fac_length > 0 & .data$Signal == "not_NS")

    confs_num <- match(confs$Bacteria, variables)

    confound <- data.frame(matrix(NA, nrow = length(confs_num), ncol = 3*max(num + 1)))
    rownames(confound) <- confs$Bacteria
    colnames(confound) <- paste(rep(c("Confounder", "Confounding_type", "Effect_size"), time = max(num + 1)), sep = "", rep(1:max(num + 1), each = 3))
    for (i in confs_num) { # loop through only the ones that do have confounding effect (Strictly/Ambiguously/Completely confounded/deconfounded)
      for (j in 1:length(sel_fac[[i]])) {# loop through different confounders
        tryCatch({
          c_name <- sel_fac[[i]][j]
          c_effectsize <- Ps_effectsize[i, c_name]
          c_type <- if (is.null(Ps_conf_model_unlist[i, sel_fac[[i]][j]])) { # if Ps_conf_model_unlist is null
            print("Strictly deconfounded")
          } else { #Ps_conf_model_unlist isn't  null
            if (Ps_conf_model_unlist[i, sel_fac[[i]][j]] < posthoc_q & !is.na(Ps_conf_model_unlist[i, sel_fac[[i]][j]])) {# Confounding model p < posthoc_q
              print("Strictly deconfounded")
            } else {# Confounding model p >= posthoc_q
              if (Ps_conf_inv_model_unlist[i, sel_fac[[i]][j]] < posthoc_q & !is.na(Ps_conf_inv_model_unlist[i, sel_fac[[i]][j]])) {# Inverse confounding model p < posthoc_q
                print("Confounded")
              } else {# Inverse onfounding model p >= posthoc_q
                print("Ambiguously deconfounded")
              }
            }}
          confound[match(i, confs_num), 3*j-2] <- c_name
          confound[match(i, confs_num), 3*j-1] <- c_type
          confound[match(i, confs_num), 3*j] <- c_effectsize
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
    #write.table(x = confound, file = paste0(output_tag, "_confounders.txt"), sep = "\t",
    #            row.names = T, col.names = NA, quote = F)

    # Create the result table
    result_table <- data.frame(matrix(nrow = length(variables), ncol = 1))
    row.names(result_table) <- variables
    colnames(result_table) <- "Correlation"

    # The final determination of this feature's signal
    final_sig <- as.data.frame(matrix(nrow = N, ncol = 1, data = NA))
    row.names(final_sig) <- variables
    colnames(final_sig) <- "Final_sig"
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] >= model_q | is.na(Ps_null_model_fdr[i, 1])) {# Null time model q >= model_q is NS
        final_sig[i, 1] = "NS"
      } else { # Null time model q < model_q
        if ((p_poho[i, 1]) >= posthoc_q | is.na(p_poho[i, 1])) {# Post-hoc q >= posthoc_q
          final_sig[i, 1] = "NS"
        } else {# Post-hoc q < posthoc_q
          if (length(sel_fac[[i]]) == 0) {# No sel_fac, meaning no confounding effect
            if (Ps_null_model[i, 2] == "Good" & !is.na(Ps_null_model[i, 2])) {# Proper confidence interval
              final_sig[i, 1] = "OK"
            } else {# Improper confidence interval
              final_sig[i, 1] = "OK but doubtful"
            }
          } else {# There are sel_fac, so decide the signal based on confound table
            subconfound <- as.data.frame(confound[str_which(string = rownames(confound), pattern = as.character(variables[i])), ])
            subconfound_columns <- subconfound[ , str_which(string = colnames(confound), pattern = "Confounding_type")]
            if (sum(str_detect(string = subconfound_columns, pattern = "Confounded"), na.rm = T) > 0) { # There is "confounding" signals
              final_sig[i, 1] = "Confounded"
            } else if(sum(str_detect(string = subconfound_columns, pattern = "Confounded"), na.rm = T) == 0){ # There isn't "confounding" signals
              if (sum(str_detect(string = subconfound_columns, pattern = "Ambiguously deconfounded"), na.rm = T) > 0) { # If there is "ambiguously deconfounded" signals
                final_sig[i, 1] = "Ambiguously deconfounded"
              } else { # There isn't "ambiguously deconfounded" signals
                final_sig[i, 1] = "OK and strictly deconfounded"
              }
            }
          }
        }
      }
    }

    # The effect type determination
    effect_type <- as.data.frame(matrix(nrow = N, ncol = 1, data = NA))
    row.names(effect_type) <- variables
    colnames(effect_type) <- "Effect"
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] < model_q & !is.na(Ps_null_model_fdr[i, 1])) { # Null time model q < model_q and isn't NA
        if (p_poho[i, 1] < posthoc_q & !is.na(p_poho[i, 1])) {# Post-hoc q < posthoc_q and isn't NA
          if (assoc[i, 1] > 0 & !is.na(assoc[i, 1])) {# Assoc > 0 and isn't NA
            effect_type[i, 1] = "Enriched"
          } else if (assoc[i, 1] < 0 & !is.na(assoc[i, 1])) {# Assoc < 0 and isn't NA
            effect_type[i, 1] = "Decreased"
          }
        } else {
          effect_type[i, 1] = "NS"
        }
      } else {
        effect_type[i, 1] = "NS"
      }
    }

    # Bind the columns together
    result_table <- cbind(prevalence, mean_abundance, final_sig, effect_type, assoc, Ps_null_model_fdr, p_poho)
    colnames(result_table) <- c("Prevalence_percentage", "Mean_abundance", "Signal", "Effect", "Effect_size", "Null_time_model_q", "Post-hoc_q")

  } else if (variable_col-1-2-length(not_used) == 0) {
    result_table <- data.frame(matrix(nrow = length(variables), ncol = 1))
    row.names(result_table) <- variables
    colnames(result_table) <- "Correlation"

    # The final determination of this feature's signal
    final_sig <- as.data.frame(matrix(nrow = N, ncol = 1, data = NA))
    row.names(final_sig) <- variables
    colnames(final_sig) <- "Final_sig"
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] >= model_q | is.na(Ps_null_model_fdr[i, 1])) {# Null time model q >= model_q is NS
        final_sig[i, 1] = "NS"
      } else { # Null time model q < model_q
        if (p_poho[i, 1] >= posthoc_q | is.na(p_poho[i, 1])) {# Post-hoc q >= posthoc_q
          final_sig[i, 1] = "NS"
        } else {# Post-hoc q < posthoc_q
          if (Ps_null_model[i, 2] == "Good" & !is.na(Ps_null_model[i, 2])) {# Proper confidence interval
            final_sig[i, 1] = "OK"
          } else {# Improper confidence interval
            final_sig[i, 1] = "OK but doubtful"
          }
        }
      }
    }
    # The effect type determination
    effect_type <- as.data.frame(matrix(nrow = N, ncol = 1, data = NA))
    row.names(effect_type) <- variables
    colnames(effect_type) <- "Effect"
    for (i in 1:N) {
      if (Ps_null_model_fdr[i, 1] < model_q & !is.na(Ps_null_model_fdr[i, 1])) { # Null time model q < model_q and isn't NA
        if (p_poho[i, 1] < posthoc_q & !is.na(p_poho[i, 1])) {# Post-hoc q < posthoc_q and isn't NA
          if (assoc[i, 1] > 0 & !is.na(assoc[i, 1])) {# Assoc > 0 and isn't NA
            effect_type[i, 1] = "Enriched"
          } else if (assoc[i, 1] < 0 & !is.na(assoc[i, 1])) {# Assoc < 0 and isn't NA
            effect_type[i, 1] = "Decreased"
          }
        } else {
          effect_type[i, 1] = "NS"
        }
      } else {
        effect_type[i, 1] = "NS"
      }
    }
    # Bind the columns together
    result_table <- cbind(prevalence, mean_abundance, final_sig, effect_type, assoc, Ps_null_model_fdr, p_poho)
    colnames(result_table) <- c("Prevalence_percentage", "Mean_abundance", "Signal", "Effect", "Effect_size", "Null_time_model_q", "Post-hoc_q")
  }
  #write.table(x = result_table, file = paste0(output_tag, "_result_table.txt"), sep = "\t",
  #           row.names = T, col.names = NA, quote = F)
  if (variable_col-1-2-length(not_used) > 0) {
    return(list(confound, result_table))
  } else if (variable_col-1-2-length(not_used) == 0) {
    return(result_table)
  }
}
