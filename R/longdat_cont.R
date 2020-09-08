#' Longitudinal analysis with time as continuous variable
#' @description
#' Longdat_cont calculates the p values, effect sizes and discover confounding effects of time variables from longitudinal data.
#' @param input A character vector. This is the path to a txt file with the first column as "Individual", and all the dependent variables (ex: bacteria)
#'         should be at the end of the table. The time variable here should be continuous, if time is discrete, please apply longdat_disc() instead.
#'         Please avoid symbols in the names of potential confounders (confounders are any column apart from individual, test_var and dependent variables).
#' @param data_type The data type of the dependent variable. Can either be "proportion", "measurement", "count", "binary", "ordinal" or "others".
#'        Proportion (or ratio) data range from 0 to 1. Measurement data are continuous and can be measured at finer and finer scale.
#'        Count data consist of discrete non-negative integers resulted from counting. Binary data are the data of sorting things into
#'        one of two mutually exclusive categories. Ordinal data consist of ranks. Any data that doesn't belong to the previous categories should be classified as "others".
#' @param test_var The name of the independent variable you are testing for, should be a character vector (ex: c("Time"))
#'        identical to its column name and make sure there is no space in it.
#' @param variable_col The column number of the position where the dependent variable columns (ex: bacteria) start in the table
#' @param fac_var The column numbers of the position where the columns that aren't
#'         numerical  (e.g. characters, categorical numbers, ordinal numbers), should be a numerical vector (ex: c(1, 2, 5:7))
#' @param not_used The column position of the columns not are irrevelant and can be ignored when in the analysis.
#'        This should be a number vector, and the default is NULL.
#' @param adjustMethod Multiple testing p value correction. Choices are the ones in p.adjust(), including
#'         "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr". The default is "fdr".
#' @param model_q The threshold for significance of model test after multiple testing correction.
#'                The default is 0.1.
#' @param posthoc_q The threshold for significance of post-hoc test after multiple testing correction.
#'                The default is 0.05.
#' @param output_tag The name tag for the output files. This should be a character vector.
#' @param theta_cutoff Required when the data type is set as "count". Variable with theta value from negative binomial regression
#'        larger than or equal to the cutoff will be filtered out if it also doesn't meet the non-zero count threshold. The default is 2^20.
#' @param nonzero_count_cutoff1 Required when the data type is set as "count". Variable with non-zero counts lower than or equal to this value
#'        will be filtered out if it also doesn't meet the theta threshold. The default is 9.
#' @param nonzero_count_cutoff2 Required when the data type is set as "count". Variable with non-zero counts lower than or equal to this value
#'        will be filtered out. The default is 5.
#' @export
#' @import lme4
#' @import tidyverse
#' @import reshape2
#' @import orddom
#' @import lmtest
#' @import glmmTMB
#' @import bestNormalize
#' @import MASS
#' @details
#' The brief workflow of longdat_cont() is as below:
#'
#' When there's no potential confounder in the input data (confounders are anything apart from individual, test_var and dependent variables):
#' First, the model test tests the significance of test_var on dependent variables. Different generalized linear mixed effect models are implemented
#' for different types of dependent variable. Negative binomial mixed model for "count", linear mixed model (dependent variables normalized first) for "measurement",
#' beta mixed model for "proportion", Binary logistic mixed model for "binary", and proportional odds logistic mixed model for "ordinal". Then, post-hoc test
#' (Spearman's correlation test) on the model is done.
#'
#' When there are potential confounders:
#' After the model test and post-hoc test described above, a confounding model test will be added to the work flow. The potential confounders will be added to the model
#' one by one and test for its significance on each dependent variable. The rest are the same as the description above.
#' @return
#' The result table will be in your output directory. If there are confounders in the input,
#' there will be another table called "confounders". The detailed description is as below.
#'
#' Result table
#'
#' 1. The first column: The dependent variables in the input data. This can be used as row name when being imported into R.
#'
#' 2. Prevalence_percentage: The percentage of each dependent variable present across individuals and time points
#'
#' 3. Mean_abundance: The mean value of each dependent variable across individuals and time points
#'
#' 4. Signal: The final decision of the significance of the test_var (independent variable) on each dependent variable.
#'    NS: Non-significant, meaning there’s no effect of time.
#'
#'    OK: There’s an effect of time and there’s no potential confounder.
#'
#'    OK but doubtful: There’s an effect of time and there’s no potential confounder, however the confidence interval of the test_var
#'                             estimate in the model test covers zero, and thus it is doubtful of this signal.
#'
#'    OK and strictly deconfounded: There are potential confounders, however there’s an effect of time and
#'                                          it is independent of those of confounders.
#'
#'    Ambiguously deconfounded: There are potential confounders, and it isn’t possible to conclude
#'                                      whether the effect is resulted from time or confounders.
#'
#'    Confounded: There’s an effect of time, but it can be reduced to the confounder effects.
#'
#' 5. Effect: This column contains the value of each dependent variable decreases/increases/NS(non-significant) along the time. A positive correlation between
#'             with time dependent variable value yields "increase", while a negative correlation yields "decrease". NS means no significant correlation.
#'
#' 6. EffectSize: This column reports the correlation coeffecient (Spearman's rho) between each dependent variable value and time.
#'
#' 7. Null_time_model_q: This column shows the multiple-comparison-adjusted p values (Wald test) of the significance of test_var in the models.
#'
#' 8. Post-hoc_q: These are the multiple-comparison-adjusted p values from the post-hoc test (Spearman's correlation test) of the model.

#'
#' Confounder table
#'
#' The first column contains the dependent variables in the input data. This can be used as row name when being imported into R.
#' Then every 3 columns are a group. Confounder column shows the confounder's name; Confounding_type column shows how test_var is
#' confounded with this confounder; Effect_size column shows the effect size of dependent variable value between different values of
#' confounders. Due to the different number of confounders for each dependent variable, there may be NAs in the table and they can
#' simply be ignored. If the confounder table is totally empty, this means that there are no confounders detected.
#' @examples
#' # Get the path of example dataset
#' system.file("Fasting_cont.txt", package = "longdat")
#' # Paste the directory to the input below
#' longdat_cont(input = "your_path_to/Fasting_cont.txt", data_type = "count",
#'              test_var = "Day", variable_col = 7, fac_var = c(1, 3),
#'              output_tag = "longdat_cont_example")

longdat_cont <- function(input, data_type, test_var, variable_col, fac_var, not_used = NULL,
                         output_tag, adjustMethod = "fdr", model_q = 0.1,
                         posthoc_q = 0.05, theta_cutoff = 2^20, nonzero_count_cutoff1 = 9, nonzero_count_cutoff2 = 5) {
  if (missing(input)) {
    stop('Error! Necessary argument "input" missing.')
  }
  if (missing(data_type)) {
    stop('Error! Necessary argument "data_type" is missing.')
  }
  if (missing(test_var)) {
    stop('Error! Necessary argument "test_var" missing.')
  }
  if (missing(variable_col)) {
    stop('Error! Necessary argument "variable_col" missing.')
  }
  if (missing(fac_var)) {
    stop('Error! Necessary argument "fac_var" missing.')
  }
  if (missing(output_tag)) {
    stop('Error! Necessary argument "output_tag" missing.')
  }

  library(lme4)
  library(tidyverse)
  library(reshape2)
  library(orddom)
  library(lmtest)
  library(glmmTMB)
  library(bestNormalize)
  library(MASS)

  data <- read.table (file = input, header = T, sep = "\t", check.names = F, stringsAsFactors = F)

  # Remove the features (bacteria) whose column sum is 0
  values <- as.data.frame(data[ , variable_col:ncol(data)])
  values <- as.data.frame(apply(values, 2, as.numeric))
  values <- values[ , colSums(values,  na.rm = T) > 0]
  data <- as.data.frame(cbind(data[ , 1:(variable_col-1)], values))
  mean_abundance <- round(colSums(values, na.rm = T)/nrow(data), 3)
  prevalence <- c()
  for (i in 1:ncol(values)) {
    prevalence[i] <- round(sum(values[ , i] > 0, na.rm = T)/nrow(data) * 100, 3)
  }

  # Here value column is defined by the user
  predictor_names <- (colnames(data))[1: (variable_col - 1)]
  melt_data <- reshape2::melt (data, id = predictor_names)
  # Omit the rows whose value column equals to NA
  melt_data <- melt_data %>% tidyr::drop_na(value)
  # Remove all dots in the bacteria name or it will cause problem
  melt_data$variable <- gsub(".", "_", melt_data$variable, fixed = TRUE)

  # Remove all symbols from the colnames of melt_data
  n1 <- gsub(x = colnames(melt_data), pattern = "-", replacement = "_", fixed = T)
  n2 <- gsub(x = n1, pattern = "#", replacement = "_", fixed = T)
  n3 <- gsub(x = n2, pattern = "%", replacement = "_", fixed = T)
  n4 <- gsub(x = n3, pattern = " ", replacement = "_",fixed = T)
  n5 <- gsub(x = n4, pattern = "(", replacement = "_", fixed = T)
  n6 <- gsub(x = n5, pattern = ")", replacement = "", fixed = T)
  n7 <- gsub(x = n6, pattern = "%", replacement = "percent", fixed = T)
  n8 <- gsub(x = n7, pattern = ".", replacement = "_", fixed = T)
  n9 <- gsub(x = n8, pattern = "&", replacement = "and", fixed = T)
  n10 <- gsub(x = n9, pattern = "$", replacement = "_", fixed = T)
  n11 <- gsub(x = n10, pattern = "@", replacement = "at", fixed = T)
  n12 <- gsub(x = n11, pattern = "!", replacement = "_", fixed = T)
  n13 <- gsub(x = n12, pattern = "+", replacement = "add", fixed = T)
  n14 <- gsub(x = n13, pattern = "/", replacement = "_", fixed = T)
  n15 <- gsub(x = n14, pattern = "?", replacement = "_", fixed = T)
  n16 <- gsub(x = n15, pattern = ",", replacement = "_", fixed = T)
  n17 <- gsub(x = n16, pattern = ">", replacement = "_", fixed = T)
  n18 <- gsub(x = n17, pattern = "<", replacement = "_", fixed = T)
  n19 <- gsub(x = n18, pattern = "*", replacement = "_", fixed = T)
  colnames(melt_data) <- n19

  # Make sure that all the columns are in the right class
  # Columns mentioned in fac_var, and the second last column in melt data are factors
  # Columns not in fac_var, and the last column in melt data are numerical numbers
  "%notin%" <- Negate("%in%")
  num_var <- c(which(1:(ncol(melt_data)-2) %notin% fac_var), ncol(melt_data))
  fac_var <- c(fac_var, ncol(melt_data)-1)
  for (i in fac_var) {
    melt_data[ ,i] <- as.factor(melt_data[ ,i])
  }
  for (i in num_var) {
    melt_data[ ,i] <- as.numeric(as.character(melt_data[ ,i]))
  }

  # Remove the not-used columns in melt_data
  if (!is.null(not_used)) {
    melt_data <- melt_data %>% dplyr::select(-c(not_used))
  }

  # Change the first column name of melt_data to "Individual"
  colnames(melt_data)[1] <- "Individual"

  # Variables are all the bacteria taxanomies
  variables <- unique (melt_data$variable)

  # Extract all the column names from melt_data except for the last two columns
  # which are variables and values, and also exclude "Individual" and test_var
  factors <- colnames(melt_data)[-c(ncol(melt_data), ncol(melt_data)-1)]
  factors <- factors[-which(factors %in% c("Individual", test_var))]
  factor_columns <- match(factors, colnames(melt_data))
  N <- length (variables)

  ########## Calculate the p values for every factor (used for selecting factors later)
  suppressWarnings(
   if (length(factor_columns) != 0) {
    Ps <- matrix(data = NA, nrow = N, ncol = length(factors))
    Ps_effectsize <- matrix(data = NA, nrow = N, ncol = length(factors))
    rownames(Ps) <- variables
    rownames(Ps_effectsize) <- variables
    colnames(Ps) <- factors
    colnames(Ps_effectsize) <- factors
    col_NA <- c()

    for (i in 1:N) {# loop through all variables
      aVariable = variables[i]
      print(i)
      print(aVariable)
      subdata <- subset(melt_data, variable == aVariable)

      for (j in 1:length(factor_columns)) {# loop through all factor columns

        # If the factor has only two kinds of values, run wilcoxon test
        if (length(unique(subdata[ , factor_columns[j]])) == 2) {
          # sub1 selects all the data that has the first kind of value
          sub1 <- subdata[which(subdata[, factor_columns[j]] == unique(subdata[ , factor_columns[j]])[1]),]
          # sub2 selects all the data that has the second kind of value
          sub2 <- subdata[which(subdata[, factor_columns[j]] == unique(subdata[ , factor_columns[j]])[2]),]
          p <- wilcox.test(sub1$value, sub2$value, paired = F)$p.value
          Ps[i,j] <- p
          d <- as.numeric(orddom::dmes(x = sub1$value, y = sub2$value)$dc)
          Ps_effectsize[i, j] <- d

          #If the factor has more than two kinds of values and is continuous (it's class is "numeric"), run spearman's correlation test
        } else if (length(unique(subdata[ , factor_columns[j]])) > 2 & class(subdata[ , factor_columns[j]]) == "numeric") {
          p <- cor.test(subdata[ , factor_columns[j]], subdata$value, method = "spearman")$p.value
          Ps[i,j] <- p
          d <- cor.test(subdata[ , factor_columns[j]], subdata$value, method = "spearman")$estimate
          Ps_effectsize[i, j] <- d

          # If the factor has more than two kinds of values and it's class is "factor" (not numeric), run kruskal-wallis test
        } else if (length(unique(subdata[ , factor_columns[j]])) > 2 & class(subdata[ , factor_columns[j]]) == "factor") {
          fmla <- as.formula(paste("value ~ ", colnames(subdata)[factor_columns[j]], sep = ""))
          p <- as.list(kruskal.test(fmla , data = subdata))$p.value
          Ps[i,j] <- p
          d <- rstatix::kruskal_effsize(data = subdata, formula = fmla)$effsize
          Ps_effectsize[i, j] <- d

          # If the factor has only one value
        } else if (length(unique(subdata[ , factor_columns[j]])) < 2) {
          Ps[i,j] <- "NA"
          Ps_effectsize[i, j] <- "NA"

          # Record the column numbers who are "NA"
          col_NA <- c(col_NA, j)
        }
      }
    }

    Ps_ori <- Ps # Save the original Ps matrix to Ps_ori

    # Exclude columns containing NA, and then turn the class of matrix into numeric
    if (length(col_NA) > 0) {
      Ps <- Ps[ , -unique(col_NA)]
    }
    mode(Ps) <- "numeric"

    if (length(col_NA) > 0) {
      Ps_effectsize <- Ps_effectsize[ , -unique(col_NA)]
    }
    mode(Ps_effectsize) <- "numeric"

    # Extract the selected factors whose p value < 0.05
    sel_fac_ini <- list() # Make an empty list first
    for (i in 1:N) {  # loop through all variables
      facs <- c()
      for (j in 1:(ncol(Ps))) {
        if (Ps[i, j] < 0.05 & !is.na(Ps[i, j])) {
          fac <- colnames(Ps)[j]
          facs <- c(facs, fac)
        }
      }
      if (!is.null(facs)) {
        sel_fac_ini[[i]] <- facs
      } else {
        sel_fac_ini[[i]] <- NA
      }
    }
    names(sel_fac_ini) <- variables
    sel_fac_ini <- lapply(sel_fac_ini, function(x) x[!is.na(x)])

    # Here remove the factors that are fixed to the individual(eg. age, sex...)
    sel_fac <- list()
    for (i in 1:N) { # loop through all variables
      aVariable = variables[i]
      print(i)
      subdata <- subset(melt_data, variable == aVariable)
      facs <- c()
      for (k in 1:length(sel_fac_ini[[i]])) {# Loop through all sel_fac
        fac <- c()
        for (m in 1:length(unique(subdata$Individual))) { # Loop through all idividuals
          subdata_individual <- subset(subdata, Individual == (unique(subdata$Individual)[m]))
          if (class(sel_fac_ini[[i]]) != "logical") {
            if (length(unique(subdata_individual[ , (sel_fac_ini[[i]][k])])) == 1) {
              f <- 0 }
            else {
              f <- 1 }
          } else if (class(sel_fac_ini[[i]]) == "logical"){
            f <- 0
          }
          fac <- c(fac, f)
        }
        if (sum(fac) == 0) { # If all individuals have the same value on this factor
          ff <- NA
        } else { # If there's any individual having different value on this factor
          ff <- sel_fac_ini[[i]][k]}
        facs <- c(facs, ff)
      }
      sel_fac[[i]] <- facs
    }
    names(sel_fac) <- variables
    sel_fac <- lapply(sel_fac, function(x) x[!is.na(x)])
    print("Finished selecting factors.")
   }
  )

  ################## Null Model Test #################
  Ps_null_model <- as.data.frame(matrix(data = NA, nrow = N, ncol = 2))
  if (data_type == "count") {
    Theta <- c()
  }
  suppressWarnings(
    for (i in 1:N) { # loop through all variables
      aVariable = variables[i]
      print(i)
      subdata <- subset(melt_data, variable == aVariable)

      tryCatch({
        if (data_type %in% c("measurement", "others")) {
          subdata <- subdata %>% mutate(value_norm = bestNormalize(value)$x.t)
        }
        if (data_type == "count") {
          # Negative binomial
          fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
          m2 <- glmmTMB(formula = fmla2, data = subdata, family = nbinom2, REML = F)
          # Extract dispersion theta out of model
          Theta[i] <- sigma(m2)
        } else if (data_type == "proportion") {
          fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
          m2 <- glmmTMB(fmla2, data = subdata, family = beta_family(), REML = F)
        } else if (data_type %in% c("measurement", "others")) {
          fmla2 <- as.formula(paste("value_norm ~ (1|Individual) +", test_var))
          m2 <- lme4::lmer(data = subdata, fmla2, REML = F)
        } else if (data_type == "binary") {
          fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
          m2 <- glmmTMB(fmla2, data = subdata, family = "binomial", REML = F)
        } else if (data_type == "ordinal") {
          fmla2 <- as.formula(paste("as.factor(value) ~ (1| Individual) +", test_var))
          m2 <- polr(fmla2, data = subdata, method = "logistic")
        }

        # Wald Chisq test
        Ps_null_model[i, 1] <- car::Anova(m2, type=c("II"),  test.statistic=c("Chisq"))$"Pr(>Chisq)"

        # Calculate confidence interval for test_var
        # Followed by the determination of the signs. For "- -" or "+ +", it means that the doesn't span 0.
        # But "- +" means that the CI spans 0. Here I sum the signs over the row, sum != 0 means it's OK.
        remove(ci)
        ci <- as.data.frame(confint(m2))
        ci <- ci[str_detect(row.names(ci), test_var), 1:2]
        if (nrow(ci) > 1) {
          if (any(apply(sign(ci), 1, sum, na.rm = T) != 0)) {
            Ps_null_model[i, 2] <- "Good"
          } else {
            Ps_null_model[i, 2] <- "Bad"
          }
        } else if (nrow(ci) == 1) {
          if (sign(ci)[1] == sign(ci)[2]) {
            Ps_null_model[i, 2] <- "Good"
          } else {
            Ps_null_model[i, 2] <- "Bad"
          }
        }
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })
  rownames(Ps_null_model) <- gsub(".", "_", variables, fixed = TRUE)
  colnames(Ps_null_model) <- c("Null_model_p", "Signal_of_CI_signs")
  if (data_type == "count") {
    Ps_null_model <- Ps_null_model %>%
      rownames_to_column() %>%
      mutate(NB_theta = Theta) %>%
      column_to_rownames()
  }

  ################## Confounding model test ###############
  if (variable_col-1-2-length(not_used) > 0) {
    Ps_conf_model <- list()
    Ps_inv_conf_model <- list()
    suppressWarnings(
      for (i in 1:N) { # loop through all variables
        aVariable = variables[i]
        print(i)
        subdata <- subset(melt_data, variable == aVariable)
        tryCatch({
        if (data_type %in% c("measurement", "others")) {
          subdata <- subdata %>% mutate(value_norm = bestNormalize(value)$x.t)
        }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        ps_lm <- c()
        ps_lm_inv <- c()
        if (length(sel_fac[[i]]) >= 1) {
          #tryCatch skips error in for loop
          tryCatch({
            for (k in 1:length(sel_fac[[i]])) { # Sel_fac exists
              if (data_type == "count") {
                fmla2 <- as.formula(paste("value ~ (1| Individual) +" , sel_fac[[i]][k], "+", test_var))
                m2 <- glmmTMB(formula = fmla2, data = subdata, family = nbinom2, REML = F)
              } else if (data_type == "proportion") {
                fmla2 <- as.formula(paste("value ~ (1| Individual) +" , sel_fac[[i]][k], "+", test_var))
                m2 <- glmmTMB(fmla2, data = subdata, family = beta_family(), REML = F)
              } else if (data_type %in% c("measurement", "others")) {
                fmla2 <- as.formula(paste("value_norm ~ (1| Individual) +" , sel_fac[[i]][k], "+", test_var))
                m2 <- lme4::lmer(data = subdata, fmla2, REML = F)
              } else if (data_type == "binary") {
                fmla2 <- as.formula(paste("value ~ (1| Individual) +" , sel_fac[[i]][k], "+", test_var))
                m2 <- glmmTMB(fmla2, data = subdata, family = "binomial", REML = F)
              } else if (data_type == "ordinal") {
                fmla2 <- as.formula(paste("as.factor(value) ~ (1| Individual) +" , sel_fac[[i]][k], "+", test_var))
                m2 <- polr(fmla2, data = subdata, method = "logistic")
              }

              # Wald Chisq test
              p_lm <- car::Anova(m2, type=c("II"), test.statistic=c("Chisq"))$`Pr(>Chisq)`[2]
              names(p_lm) <- paste(sel_fac[[i]][k])
              ps_lm <- c(ps_lm, p_lm)

              p_lm_inv <- car::Anova(m2, type=c("II"), test.statistic=c("Chisq"))$`Pr(>Chisq)`[1]
              names(p_lm_inv) <- paste(sel_fac[[i]][k])
              ps_lm_inv <- c(ps_lm_inv, p_lm_inv)
            }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        } else if (length(sel_fac[[i]]) == 0) { # No sel_fac existing
          ps_lm <- "No selected factors for confounding model test"
          ps_lm_inv <- "No selected factors for inverse confounding model test"
        }
        Ps_conf_model[[i]] <- ps_lm
        Ps_inv_conf_model[[i]] <- ps_lm_inv
      }
    )
    names(Ps_conf_model) <- gsub(".", "_", variables, fixed = TRUE)
    names(Ps_inv_conf_model) <- gsub(".", "_", variables, fixed = TRUE)
    }

  ############## Unlist the Ps_conf_model #################
  suppressWarnings(
    # Spread the Ps_conf_model list
    if (variable_col-1-2-length(not_used) > 0) {
      unlist = unlist(Ps_conf_model, use.names=FALSE) # output only the numbers
      mode(unlist) <- "numeric" # it's numeric now!
      unlist_noNA <-unlist[!is.na(unlist)] # Exclude NA in unlist
      unlist_withname = unlist(Ps_conf_model) # numbers with names
      unlist_name <- names(unlist_withname) # outputs the name of every element
      # if there's NA, then remove names of NA value
      if (TRUE %in% is.na(unlist)) {
        unlist_name_noNA <- unlist_name[-which(is.na(unlist))]
      } else {
        unlist_name_noNA <- unlist_name
      }

      # Create a empty matrix first
      ps_lm_m <- data.frame(matrix(nrow = length(unlist), ncol = 2))
      # Specify column names and fill in the values
      colnames(ps_lm_m) <- c("name", "p_lm")
      # Fill in the values
      ps_lm_m[ , 1] <- unique(unlist_name)
      ps_lm_m[ , 2] <- unlist

      # Split the names of each p_lm, ex: Actinobacteria.PPI becomes "Actinobacteria" "PPI"
      library(stringr)
      out <- stringr::str_split_fixed(ps_lm_m$name, pattern = stringr::fixed("."), n = 2)
      # Add the 2 new columns with split names to the matrix
      ps_lm_m2 <- cbind(out, ps_lm_m)
      # Remove the original unsplit name column, and then name the columns
      ps_lm_m2 $name = NULL
      colnames(ps_lm_m2) <- c("feature", "factor", "p_lm")
      # Reshaping the matrix so that can do fdr correction in the next step
      cast.ps_lm_m2 <- reshape2::dcast(data = ps_lm_m2, feature ~ factor)
      if (sum(is.na(cast.ps_lm_m2$Var.2)) == N) {
        cast.ps_lm_m2 <- as.data.frame(cast.ps_lm_m2[ , -2])
      }
      rownames(cast.ps_lm_m2) <- cast.ps_lm_m2$feature
      cast.ps_lm_m2$feature <- NULL
      Ps_conf_model_unlist <- cast.ps_lm_m2
      #Reorder the rows to make it follow the order of variables
      Ps_conf_model_unlist <- as.data.frame(Ps_conf_model_unlist[match(variables,rownames(Ps_conf_model_unlist)), ])
    }
  )
  ############## Unlist the Ps_inv_conf_model #################
  suppressWarnings(
    if (variable_col-1-2-length(not_used) > 0) {
      unlist_inv = unlist(Ps_inv_conf_model, use.names=FALSE) # output only the numbers
      mode(unlist_inv) <- "numeric" # it's numeric now!
      unlist_inv_noNA <-unlist[!is.na(unlist_inv)] # Exclude NA in unlist
      unlist_inv_withname = unlist(Ps_inv_conf_model) # numbers with names
      unlist_inv_name <- names(unlist_inv_withname) # outputs the name of every element
      # if there's NA, then remove names of NA value
      if (TRUE %in% is.na(unlist_inv)) {
        unlist_inv_name_noNA <- unlist_inv_name[-which(is.na(unlist_inv))]
      } else {
        unlist_inv_name_noNA <- unlist_inv_name
      }

      # Create a empty matrix first
      ps_inv_lm_m <- data.frame(matrix(nrow = length(unlist_inv), ncol = 2))
      # Specify column names and fill in the values
      colnames(ps_inv_lm_m) <- c("name", "p_lm_inv")
      # Fill in the values
      ps_inv_lm_m[ , 1] <- unique(unlist_inv_name)
      ps_inv_lm_m[ , 2] <- unlist_inv

      # Split the names of each p_lm, ex: Actinobacteria.PPI becomes "Actinobacteria" "PPI"
      library(stringr)
      out_inv <- stringr::str_split_fixed(ps_inv_lm_m$name, pattern = stringr::fixed("."), n = 2)
      # Add the 2 new columns with split names to the matrix
      ps_inv_lm_m2 <- cbind(out_inv, ps_inv_lm_m)
      # Remove the original unsplit name column, and then name the columns
      ps_inv_lm_m2 $name = NULL
      colnames(ps_inv_lm_m2) <- c("feature", "factor", "p_inv_lm")
      # Reshaping the matrix so that can do fdr correction in the next step
      cast.ps_inv_lm_m2 <- reshape2::dcast(data = ps_inv_lm_m2, feature ~ factor)
      if (sum(is.na(cast.ps_inv_lm_m2$Var.2)) == N) {
        cast.ps_inv_lm_m2 <- as.data.frame(cast.ps_inv_lm_m2[ , -2])
      }
      rownames(cast.ps_inv_lm_m2) <- cast.ps_inv_lm_m2$feature
      cast.ps_inv_lm_m2$feature <- NULL

      Ps_conf_inv_model_unlist <- cast.ps_inv_lm_m2
      #Reorder the rows to make it follow the order of variables
      Ps_conf_inv_model_unlist <- as.data.frame(Ps_conf_inv_model_unlist[match(variables,rownames(Ps_conf_inv_model_unlist)), ])
    }
  )

  ############## Post-hoc test (p value and association) #################
  # Here uses Spearman's correlation
  p_poho <- as.data.frame(matrix(nrow = length(variables), ncol = 1))
  assoc <- as.data.frame(matrix(nrow = length(variables), ncol = 1))

  for (i in 1:N) { # loop through all variables
    print(i)
    bVariable = variables[i]
    subdata <- subset(melt_data, variable == bVariable)
    # Here set the "test_var" to numeric
    subdata[ , test_var] <- as.numeric(subdata[ , test_var])
    c <- cor.test(subdata[ , test_var], subdata$value, method = "spearman")
    p_c <- c$p.value
    a_c <- c$estimate
    p_poho[i, 1] <- p_c
    assoc[i, 1] <- a_c
  }
  row.names(p_poho) <- variables
  row.names(assoc) <- variables
  colnames(p_poho) <- "p_post-hoc"
  colnames(assoc) <- "association"

  print("Finished post-hoc correlation test.")

  ######################### Remove the excluded ones #########################
  if (data_type == "count") {
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

    p_poho <- p_poho %>%
      rownames_to_column() %>%
      dplyr::filter(rowname %in% bac_include) %>%
      column_to_rownames()

    assoc <- assoc %>%
      rownames_to_column() %>%
      dplyr::filter(rowname %in% bac_include) %>%
      column_to_rownames()

    variables <- bac_include
    N <- length(variables)
  }

  if (data_type == "count" & variable_col-1-2-length(not_used) > 0) {
    sel_fac <-  sel_fac[match(bac_include, table = names(sel_fac))]
    Ps_conf_model_unlist <- Ps_conf_model_unlist %>%
      rownames_to_column() %>%
      dplyr::filter(rowname %in% bac_include) %>%
      column_to_rownames()

    Ps_conf_inv_model_unlist <- Ps_conf_inv_model_unlist %>%
      rownames_to_column() %>%
      dplyr::filter(rowname %in% bac_include) %>%
      column_to_rownames()
  }
  ################## FDR correction on Ps_null_model (model test) ##############
  # Do FDR correction (correct for the number of features)
  adjust_fun <- function(x) p.adjust(p = x, method = adjustMethod)
  suppressWarnings(
    Ps_null_model_fdr <- apply(X = Ps_null_model, MARGIN = 2, FUN = adjust_fun))
  Ps_null_model_fdr <- as.data.frame(Ps_null_model_fdr[ , 1])
  #Reorder the rows of p_lm_fdr to make it follow the order of variables
  Ps_null_model_fdr <- as.data.frame(Ps_null_model_fdr[match(variables,rownames(Ps_null_model_fdr)), ])
  rownames(Ps_null_model_fdr) <- gsub(".", "_", variables, fixed = TRUE)
  colnames(Ps_null_model_fdr) <- c("Ps_null_model_fdr")
  ############## Generate result table as output #################
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
      rownames_to_column("Bacteria") %>%
      filter(sel_fac_length > 0 & Signal == "not_NS")

    confs_num <- match(confs$Bacteria, variables)

    confound <- data.frame(matrix(NA, nrow = length(confs_num), ncol = 3*max(num + 1)))
    rownames(confound) <- rownames(confs)
    colnames(confound) <- paste(rep(c("Confounder", "Confounding_type", "Effect_size"), time = max(num + 1)), sep = "", rep(1:max(num + 1), each = 3))
    for (i in confs_num) { # loop through only the ones that do have confounding effect (Strictly/Ambiguously/Completely confounded/deconfounded)
      for (j in 1:length(sel_fac[[i]])) {# loop through different confounders
        tryCatch({
          c_name <- sel_fac[[i]][j]
          c_effectsize <- Ps_effectsize[i, c_name]
          c_type <- if (is.null(Ps_conf_model_unlist[i, sel_fac[[i]][j]])) { # if Ps_conf_model_unlist is null
            print("Strictly deconfounded")
          } else { #Ps_conf_model_unlist isn't  null
            if (Ps_conf_model_unlist[i, sel_fac[[i]][j]] < 0.05 & !is.na(Ps_conf_model_unlist[i, sel_fac[[i]][j]])) {# Confounding model p < 0.05
              print("Strictly deconfounded")
            } else {# Confounding model p >= 0.05
              if (Ps_conf_inv_model_unlist[i, sel_fac[[i]][j]] < 0.05 & !is.na(Ps_conf_inv_model_unlist[i, sel_fac[[i]][j]])) {# Inverse onfounding model p < 0.05
                print("Confounded")
              } else {# Inverse onfounding model p >= 0.05
                print("Ambiguously deconfounded")
              }
            }}
          confound[match(i, confs_num), 3*j-2] <- c_name
          confound[match(i, confs_num), 3*j-1] <- c_type
          confound[match(i, confs_num), 3*j] <- c_effectsize
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
    write.table(x = confound, file = paste0(output_tag, "_confounders.txt"), sep = "\t",
                row.names = T, col.names = NA, quote = F)

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
          } else {# There are sel_fac
            subconfound <- as.data.frame(confound[str_which(string = rownames(confound), pattern = as.character(variables[i])), ])
            subconfound_columns <- subconfound[ , str_which(string = colnames(confound), pattern = "Confounding_type")]
            if (sum(str_detect(string = subconfound_columns, pattern = "Confounded"), na.rm = T) > 0) { # There is "confounding" signals
              final_sig[i, 1] = "Confounded"
            } else { # There isn't "confounding" signals
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
        if ((p_poho[i, 1]) < posthoc_q & !is.na(p_poho[i, 1])) {# Post-hoc q >= posthoc_q
          final_sig[i, 1] = "NS"
        } else {# At least one post-hoc q < posthoc_q
          if (Ps_null_model[i, 2] == "Good") {# Proper confidence interval
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
  write.table(x = result_table, file = paste0(output_tag, "_result_table.txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)
  print("Finished! The results are now in your directory.")
}


