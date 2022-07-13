#' Longitudinal analysis with time as discrete variable
#' @description
#' longdat_disc calculates the p values, effect sizes and discover covariate
#'  effects of time variables from longitudinal data.
#' @param input A data frame with the first column as "Individual"
#' and all the columns of dependent variables (features, e.g. bacteria)
#'         at the end of the table. The time variable here should be
#'         discrete, if time is continuous, please apply longdat_cont()
#'         instead. Please avoid using characters that don't belong to
#'         ASCII printable characters for potential covariates names
#'         (covariates are any column apart from individual, test_var and
#'         dependent variables).
#' @param data_type The data type of the dependent variables (features).
#' Can either be "proportion", "measurement", "count", "binary", "ordinal" or
#' "others". Proportion (or ratio) data range from 0 to 1. Measurement data are
#'        continuous and can be measured at finer and finer scale (e.g. weight).
#'        Count data consist of discrete non-negative integers resulted from
#'        counting. Binary data are the data of sorting things into
#'        one of two mutually exclusive categories. Ordinal data consist of
#'        ranks. Any data that doesn't belong to the previous categories
#'        should be classified as "others".
#' @param test_var The name of the independent variable you are testing for,
#' should be a string (e.g. "Time")
#'        identical to its column name and make sure there is no space in it.
#' @param variable_col The column number of the position where the dependent
#' variable columns (features, e.g. bacteria) start in the table.
#' @param fac_var The column numbers of the position where the columns that
#' aren't numerical  (e.g. characters, categorical numbers, ordinal numbers).
#' This should be a numerical vector (e.g. c(1, 2, 5:7)).
#' @param not_used The column position of the columns not are irrelevant and
#' can be ignored when in the analysis. This should be a numerical vector, and
#' the default is NULL.
#' @param adjustMethod Multiple testing p value correction. Choices are the
#' ones in p.adjust(), including 'holm', 'hochberg', 'hommel', 'bonferroni',
#' 'BH', 'BY' and 'fdr'. The default is 'fdr'.
#' @param model_q The threshold for significance of model test after multiple
#' testing correction. The default is 0.1.
#' @param posthoc_q The threshold for significance of post-hoc test of the
#' model after multiple testing correction. The default is 0.05.
#' @param theta_cutoff Required when the data type is set as "count".
#'  Variable with theta value from negative binomial regression
#'        larger than or equal to the cutoff will be filtered out if it
#'        also doesn't meet the non-zero count threshold.
#'        Users can use the function "theta_plot()" to help with
#'        specifying the value for theta_cutoff. The default is 2^20.
#' @param nonzero_count_cutoff1 Required when the data type is set as "count".
#'  Variable with non-zero counts lower than or equal to this value
#'        will be filtered out if it doesn't meet the theta threshold either.
#'        Users can use the function "theta_plot()" to help with
#'        specifying the value for nonzero_count_cutoff1. The default is 9.
#' @param nonzero_count_cutoff2 Required when the data type is set as "count".
#'  Variable with non-zero counts lower than or equal to this value
#'        will be filtered out. Users can use the function "theta_plot()" to
#'        help with specifying the value for nonzero_count_cutoff2.
#'        The default is 5.
#' @param verbose A boolean vector indicating whether to print detailed
#' message. The default is TRUE.
#' @name longdat_disc
#' @export
#' @import utils
#' @import graphics
#' @import dplyr
#' @import tibble
#' @importFrom magrittr '%>%'
#' @importFrom rlang .data
#' @importFrom stats as.formula confint cor.test kruskal.test na.omit p.adjust
#'             wilcox.test
#' @details
#' The brief workflow of longdat_disc() is as below:
#'
#' When there's no potential covariates in the input data
#' (covariates are anything apart from individual, test_var and
#' dependent variables):
#' First, the model test tests the significance of test_var on dependent
#' variables. Different generalized linear mixed effect models are implemented
#' for different types of dependent variable. Negative binomial mixed model for
#'  "count", linear mixed model (dependent variables normalized first) for
#'  "measurement",
#' beta mixed model for "proportion", binary logistic mixed model for
#' "binary", and proportional odds logistic mixed model for "ordinal".
#' Then, post-hoc test ('emmeans') on
#' the model is done. When the data type is "count" mode, a control model
#'  test will be run on randomized data (the rows are shuffled). If there
#'   are false positive signals in this control model
#' test, then additional Wilcoxon post-hoc test will be done because it
#' is more conservative.
#'
#' When there are potential covariates in the input data:
#' After the model test and post-hoc test described above, a covariate
#'  model test will be added to the work flow. The potential covariates will
#'   be added to the model
#' one by one and test for its significance on each dependent variable. The
#' rest are the same as the description above.
#'
#' Also, when your data type is count data, please use set.seed() before
#' running longdat_disc() so that you can get reproducible randomized
#' negative check.
#'
#' @return
#' longdat_disc() returns a list which contains a "Result_table",
#' and if there are covariates in the input data frame,
#' there will be another table called "Covariate_table". For count mode,
#'  if there is any false positive
#' in the randomized control result, then another table named
#' "Randomized_control_table" will also be
#' generated. The detailed description is as below.
#'
#' Result_table
#'
#' 1. The first column: The dependent variables in the input data.
#' This can be used as row name when being imported into R.
#'
#' 2. Prevalence_percentage: The percentage of each dependent variable
#' present across individuals and time points.
#'
#' 3. Mean_abundance: The mean value of each dependent variable across
#' individuals and time points.
#'
#' 4. Signal: The final decision of the significance of the test_var
#'  (independent variable) on each dependent variable.
#'    NS: This represents "Non-significant",
#'    which means that there’s no effect of time.
#'
#'    OK_nc: This represents "OK and no covariate".
#'    There’s an effect of time and there’s no potential covariate.
#'
#'    OK_d: This represents "OK but doubtful".
#'    There’s an effect of time and there’s no
#'    potential covariate, however the confidence interval of the test_var
#'    estimate in the model test covers zero, and thus it is doubtful of
#'    this signal.
#'
#'    OK_nrc: This represents "OK and not reducible to covariate".
#'    There are potential covariates, however
#'    there’s an effect of time and it is independent of those of covariates.
#'
#'    EC: This represents "Entangled with covariate".
#'    There are potential covariates, and it isn’t
#'    possible to conclude whether the effect is resulted from time or
#'    covariates.
#'
#'    RC: This represents "Effect reducible to covariate".
#'    There’s an effect of time, but it can be reduced to the
#'    covariate effects.
#'
#' 5. 'Effect_a_b': The "a" and "b" here are the names of the time points.
#'  These columns describe the value of each dependent
#'                variable decreases/increases/NS(non-significant) at time
#'                 point b comparing with time point a. The number of
#'                Effect columns depends on how many combinations of
#'                 time points in the input data.
#'
#' 6. 'EffectSize_a_b': The "a" and "b" here are the names of the time points.
#'  These columns describe the effect size (Cliff's delta) of each dependent
#'                    variable between time point b and a. The number of
#'                    'EffectSize' columns depends on how many combinations of
#'                    time points in the input data.
#'
#'  7. 'Null_time_model_q': This column shows the multiple-comparison-adjusted
#'  p values (Wald test) of the significance of test_var in the models.
#'
#'  8. 'Post-hoc_q_a_b': The "a" and "b" here are the names of the time points.
#'   These are the multiple-comparison-adjusted p values from the post-hoc
#'                     test of the model. The number of Post-hoc_q columns
#'                     depends on how many combinations of time points
#'                      in the input data.
#'
#'  9. 'Wilcox_p_a_b': The "a" and "b" here are the names of the time points.
#'   These columns only appear when data type is "count" and there exist
#'   false positives in the model test on randomized data. Wilcoxon test
#'   are more conservative than the default post-hoc test ('emmeans'),
#'   and thus it is a good reference for getting a more conservative result
#'   of the significant outcomes.
#'
#' Covariate_table
#'
#' The first column contains the dependent variables in the input data.
#' This can be used as row name when being imported into R.
#' Then every 3 columns are a group. Covariate column shows the covariate's
#'  name; Covariate_type column shows how effect is affected by covariate;
#'  Effect_size column shows the effect
#' size of dependent variable value between different values of
#' covariate. Due to the different number of covariates for each dependent
#' variable, there may be NAs in the table and they can
#' simply be ignored. If the covariate table is totally empty, this means
#' that there are no covariates detected.
#'
#' Randomized_control_table (for user's reference)
#'
#' We assume that there shouldn't be positive results in the randomized control
#' test, because all the rows in the original dataset are
#' shuffled randomly. Therefore, any signal that showed significance here will
#'  be regarded as false positive. And if there's false
#' positive in this randomized control result, longdat_disc() will warn the user
#' at the end of the run. This Randomized_control table is only
#' generated when there is false positive in the randomized control test.
#'  It is intended to be a reference for users to see the effect size of
#' false positive features.
#'
#'  1. "Model_q": It shows the multiple-comparison-adjusted p values
#'   (Wald test)of the significance of test_var in the negative-binomial
#'    models in the randomized dataset. Only the features with Model_q lower
#'     than the defined model_q (default = 0.1) will be listed in this table.
#'
#'  2. Final_signal: It show the overall signal being either false positive
#'   or negative. "False positive" indicates that test_var is significant,
#'   while "Negative" indicates non-significance.
#'
#'  3. 'Signal_a_b': The "a" and "b" here are the names of the time points.
#'   These columns describe if test_var is significant on each dependent
#'   variable between each time point based on the post-hoc test p values
#'    (listed right to Signal_a_b). "False positive" indicates that test_var
#'     is significant, while "Negative" indicates non-significance. The number
#'     of Signal_a_b columns depends on how many combinations of time points
#'      in the input data.
#'
#'  4. 'Posthoc_q_a_b': The "a" and "b" here are the names of the time points.
#'  These columns describe the multiple-comparison-adjusted p values from the
#'  post-hoc test of the model between time point b and a in the randomized
#'  control dataset. The number of 'Posthoc_q_a_b' columns depends on how
#'  many combinations of time points in the input data.
#'
#'
#'  5. 'Effect_size_a_b': The "a" and "b" here are the names of the time points.
#'  These columns describe the effect size (Cliff's delta) of each dependent
#'   variable between time point b and a in the randomized control dataset.
#'   The number of Effect_size_a_b columns depends on how many combinations
#'    of time points in the input data.
#'
#' @examples
#' test_disc <- longdat_disc(input = LongDat_disc_master_table,
#' data_type = "count", test_var = "Time_point",
#' variable_col = 7, fac_var = c(1:3))


longdat_disc <- function(input, data_type, test_var, variable_col, fac_var,
                         not_used = NULL, adjustMethod = "fdr", model_q = 0.1,
                         posthoc_q = 0.05, theta_cutoff = 2^20,
                         nonzero_count_cutoff1 = 9,
                         nonzero_count_cutoff2 = 5, verbose = TRUE) {
  if (missing(input)) {
    stop('Error! Necessary argument "input" is missing.')
  }
  if (missing(data_type)) {
    stop('Error! Necessary argument "data_type" is missing.')
  }
  if (missing(test_var)) {
    stop('Error! Necessary argument "test_var" is missing.')
  }
  if (missing(variable_col)) {
    stop('Error! Necessary argument "variable_col" is missing.')
  }
  if (missing(fac_var)) {
    stop('Error! Necessary argument "fac_var" is missing.')
  }

  ############## Data preprocessing #################
  if (verbose == TRUE) {print("Start data preprocessing.")}
  preprocess_lists <- data_preprocess(input, test_var, variable_col,
                                      fac_var, not_used)
  mean_abundance <- preprocess_lists[[1]]
  prevalence <- preprocess_lists[[2]]
  variables_original <- preprocess_lists[[3]]
  melt_data <- as.data.frame((preprocess_lists[[4]]))
  variables <- preprocess_lists[[5]]
  factor_columns <- preprocess_lists[[6]]
  factors <- preprocess_lists[[7]]
  data <- preprocess_lists[[8]]
  values <- preprocess_lists[[9]]
  N <- length (variables)
  if (verbose == TRUE) {print("Finish data preprocessing.")}

  ########## Calculate the p values for every factor
  #          (used for selecting factors later)
  if (variable_col-1-2-length(not_used) > 0) {
    if (verbose == TRUE) {print("Start selecting potential covariates.")}
    factor_p_lists <- suppressWarnings(factor_p_cal(melt_data, variables,
                                                    factor_columns, factors,
                                                    data, N, verbose))
    Ps <- as.data.frame(factor_p_lists[[1]])
    Ps_effectsize <- as.data.frame(factor_p_lists[[2]])
    sel_fac <- factor_p_lists[[3]]
    if (verbose == TRUE) {print("Finished selecting potential covariates.")}
  }

  ################## Null Model Test and Post-Hoc Test #################
  if (verbose == TRUE) {print("Start null model test and post-hoc test.")}
  NuModel_lists <- NuModelTest_disc(N, data_type, test_var, melt_data,
                                    variables, verbose)
  Ps_null_model <- as.data.frame(NuModel_lists[[1]])
  Ps_poho_fdr <- as.data.frame(NuModel_lists[[2]])
  case_pairs_name <- NuModel_lists[[3]]

  if (data_type %in% c("measurement", "others")) {
    Normalize_method <- NuModel_lists[[4]]
  }
  if (verbose == TRUE) {print("Finish null model test and post-hoc test.")}

  ################## Covariate model test ###############
  if (variable_col-1-2-length(not_used) > 0) {
    if (verbose == TRUE) {print("Start covariate model test.")}
    ConModel_lists <- ConModelTest_disc(N, variables, melt_data, sel_fac,
                                        data_type, test_var, verbose)
    Ps_conf_model <- ConModel_lists[[1]]
    Ps_inv_conf_model <-ConModel_lists[[2]]
    if (verbose == TRUE) {print("Finish covariate model test.")}
  }

  ####### Unlist the Ps_conf_model and Ps_inv_conf_model ########
  suppressWarnings(
    if (variable_col-1-2-length(not_used) > 0) {
      if (verbose == TRUE) {print(
        "Start unlisting tables from covariate model result.")}
      Ps_conf_model_unlist <- unlist_table(Ps_conf_model, N, variables)
      Ps_conf_inv_model_unlist <- unlist_table(Ps_inv_conf_model, N, variables)
      if (verbose == TRUE) {print(
        "Finish unlisting tables from covariate model result.")}
      })

  ############## Calculate cliff's delta for all variables #################
  if (verbose == TRUE) {print("Start calculating effect size.")}
  cliff_lists <- cliff_cal(melt_data, Ps_poho_fdr, variables, test_var,
                           data, verbose)
  case_pairs <- cliff_lists[[1]]
  delta <- cliff_lists[[2]]
  if (verbose == TRUE) {print("Finish calculating effect size.")}

  ######### Randomized negative control model test #########
  if (data_type == "count") {
    if (verbose == TRUE) {print(
      "Start randomized negative control model test.")}
    result_neg_ctrl <- random_neg_ctrl_disc(test_var, variable_col, fac_var,
                                            not_used, factors, data, N,
                                            data_type, variables, case_pairs,
                                            adjustMethod, model_q, posthoc_q,
                                            theta_cutoff, nonzero_count_cutoff1,
                                            nonzero_count_cutoff2, verbose)
    if (verbose == TRUE) {print(
      "Finish randomized negative control model test.")}
  }

  ################### Wilcoxon post-hoc test (conditional) ###################
  suppressWarnings(
    if (data_type == "count") {
      if (verbose == TRUE) {print("Start Wilcoxon post-hoc test.")}
      wilcox_poho_lists <- wilcox_posthoc(result_neg_ctrl, model_q, melt_data,
                                          test_var, variables, data, N, verbose)
      false_pos_count <- wilcox_poho_lists[[1]]
      p_wilcox <- wilcox_poho_lists[[2]]
      if (verbose == TRUE) {print("Finish Wilcoxon post-hoc test.")}
    }
  )

  ##### Reset the variable names and covariate names to oringinal ####
  rownames(Ps_null_model) <- variables_original
  rownames(Ps_poho_fdr) <- variables_original
  rownames(delta) <- variables_original
  variables <- variables_original

  if (data_type == "count" & variable_col-1-2-length(not_used) > 0) {
    names(sel_fac) <- variables_original
    rownames(Ps_conf_model_unlist) <- variables_original
    rownames(Ps_conf_inv_model_unlist) <- variables_original
  }

  if (data_type == "count") {
    if (false_pos_count > 0) {
      rownames(p_wilcox) <- variables_original
    }}

  if (data_type %in% c("measurement", "others")) {
    Normalize_method$Feature <- variables_original
  }

  ############## Remove the excluded one when data_type = count ##########
  if (data_type == "count") {
    if (verbose == TRUE) {print(
      "Start removing the dependent variables to be exlcuded.")}
    rm_sparse_lists <- rm_sparse_disc(values, data, nonzero_count_cutoff1,
                                      nonzero_count_cutoff2, theta_cutoff,
                                      Ps_null_model, prevalence,
                                      absolute_sparsity, mean_abundance,
                                      Ps_poho_fdr, delta)
    prevalence <- rm_sparse_lists[[1]]
    absolute_sparsity <- rm_sparse_lists[[2]]
    mean_abundance <- rm_sparse_lists[[3]]
    Ps_null_model <-  rm_sparse_lists[[4]]
    Ps_poho_fdr <- rm_sparse_lists[[5]]
    delta <- rm_sparse_lists[[6]]
    variables <- rm_sparse_lists[[7]]
    bac_include <-  rm_sparse_lists[[8]]
    N <- length(variables)
    if (variable_col-1-2-length(not_used) > 0) {
      sel_fac <-  sel_fac[match(bac_include, table = names(sel_fac))]
      Ps_conf_model_unlist <- Ps_conf_model_unlist %>%
        tibble::rownames_to_column() %>%
        dplyr::filter(.data$rowname %in% bac_include) %>%
        tibble::column_to_rownames()
      Ps_conf_inv_model_unlist <- Ps_conf_inv_model_unlist %>%
        tibble::rownames_to_column() %>%
        dplyr::filter(.data$rowname %in% bac_include) %>%
        tibble::column_to_rownames()
    }
    if (false_pos_count > 0) {
      p_wilcox <- p_wilcox %>%
        tibble::rownames_to_column() %>%
        dplyr::filter(.data$rowname %in% bac_include) %>%
        tibble::column_to_rownames()
    }
    if (verbose == TRUE) {print(
      "Finish removing the dependent variables to be exlcuded.")}
  }

  ################## FDR correction on Ps_null_model (model test) ############
  # Do FDR correction (correct for the number of features)
  if (verbose == TRUE) {print(
    "Start multiple test correction on null model test p values.")}
  adjust_fun <- function(x) p.adjust(p = x, method = adjustMethod)
  suppressWarnings(
  Ps_null_model_fdr <- apply(X = Ps_null_model, MARGIN = 2,
                             FUN = adjust_fun))
  if (class(Ps_null_model_fdr)[1] != "numeric") {
    Ps_null_model_fdr <- as.data.frame(Ps_null_model_fdr[ , 1])
  } else {
    Ps_null_model_fdr <- as.data.frame(Ps_null_model_fdr[1])
    rownames(Ps_null_model_fdr) <- rownames(Ps_null_model)
  }
  #Reorder the rows of p_lm_fdr to make it follow the order of variables
  Ps_null_model_fdr <-
    as.data.frame(Ps_null_model_fdr[match(variables,
                                          rownames(Ps_null_model_fdr)), ])
  rownames(Ps_null_model_fdr) <- variables
  colnames(Ps_null_model_fdr) <- c("Ps_null_model_fdr")
  if (verbose == TRUE) {print(
    "Finish multiple test correction on null model test p values.")}

  ################## FDR correction on p_wilcoxon (conditional) ##############
  suppressWarnings(
    if (data_type == "count") {
      if (false_pos_count > 0 ) {
      if (verbose == TRUE) {print(
        "Start multiple test correction on wilcoxon test p values.")}
      if (ncol(p_wilcox) > 1) {
        p_wilcox_t <- as.data.frame(t(p_wilcox))
        p_wilcox_fdr <- apply(X = p_wilcox_t, MARGIN = 2, FUN = adjust_fun)
        p_wilcox_final <- as.data.frame(t(p_wilcox_fdr))
        new_wilcox_name <- paste0("Wilcox_", colnames(p_wilcox_final))
        colnames(p_wilcox_final) <- new_wilcox_name
      } else if (ncol(p_wilcox) == 1){ #when ncol(p_wilcox) = 1
        p_wilcox_t <- as.data.frame(t(p_wilcox))
        p_wilcox_fdr <- apply(X = p_wilcox_t, MARGIN = 2, FUN = adjust_fun)
        p_wilcox_final <- as.data.frame(p_wilcox_fdr)
        new_wilcox_name <- paste0("Wilcox_", colnames(p_wilcox_final))
        colnames(p_wilcox_final) <- new_wilcox_name
      }
      if (verbose == TRUE) {print(
        "Finish multiple test correction on wilcoxon test p values.")}
      }
    }
  )

  ############## Generate result table as output #################
  if (verbose == TRUE) {print("Start generating result tables.")}
  final_result <- final_result_summarize_disc(variable_col, N,
                                              Ps_conf_inv_model_unlist,
                                              variables, sel_fac,
                                              Ps_conf_model_unlist,
                                              model_q, posthoc_q,
                                              Ps_null_model_fdr,
                                              Ps_null_model, delta,
                                              case_pairs, prevalence,
                                              mean_abundance,
                                              Ps_poho_fdr, not_used,
                                              Ps_effectsize, case_pairs_name,
                                              data_type,
                                              false_pos_count, p_wilcox_final)
  if (variable_col-1-2-length(not_used) > 0) {
    Confounder_table <- final_result[[1]]
    Result_table <- final_result[[2]]
  } else if (variable_col-1-2-length(not_used) == 0) {
    Result_table <- final_result
  }
  print("Finished successfully!")
  if (data_type == "count") {
    if (false_pos_count > 0) {
      print(paste0("Attention! Since there are false ",
            "positives in the randomized control test, ",
            "it's recommended to check the Wilcoxon post-hoc ",
            "p values of significant signals in the ",
            "Result_table to get a more conservative result. ",
            "See the documentation for more details."))
    }
  }

  if (data_type == "count") {
    if (nrow(result_neg_ctrl) > 0 & variable_col-1-2-length(not_used) > 0) {
      return(list(Result_table = Result_table,
                  Covariate_table = Confounder_table,
                  Randomized_control_table = result_neg_ctrl))
    } else if (nrow(result_neg_ctrl) > 0 &
               variable_col-1-2-length(not_used) == 0) {
      return(list(Result_table = Result_table,
                  Randomized_control_table = result_neg_ctrl))
    } else if (nrow(result_neg_ctrl) == 0 &
               variable_col-1-2-length(not_used) > 0) {
      return(list(Result_table = Result_table,
                  Covariate_table = Confounder_table))
    } else if (nrow(result_neg_ctrl) == 0 &
               variable_col-1-2-length(not_used) == 0) {
      return(list(Result_table = Result_table))}
  } else if(data_type %in% c("measurement", "others")) {
    if (variable_col-1-2-length(not_used) > 0) {
      return(list(Result_table = Result_table,
                  Covariate_table = Confounder_table,
                  Normalize_method = Normalize_method))
    } else if (variable_col-1-2-length(not_used) == 0) {
      return(Result_table = Result_table,
             Normalize_method = Normalize_method)
    }
  } else if(data_type %in% c("proportion", "binary", "ordinal")) {
    if (variable_col-1-2-length(not_used) > 0) {
      return(list(Result_table = Result_table,
                  Covariate_table = Confounder_table))
    } else if (variable_col-1-2-length(not_used) == 0) {
      return(Result_table = Result_table)
    }
  }
}


