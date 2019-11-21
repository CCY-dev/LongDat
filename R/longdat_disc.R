#' Calculate the p values and effect size from longitudinal data
#' @param input a table with the first column as "Individual", and all the features
#'         should be at the end of the table. Also, place the columns that aren't
#'         numerical (e.g. categorical, ordinal) in the front of the table, right after "Individual"
#' @param test_var the name of the variable you are testing for, should be a character vector
#'        identical to its column name and make sure there is no space in it.
#' @param value_column the number of the position where the features columns start in the table
#' @param factor_column the number of the position where the columns that aren't
#'         numerical  (e.g. categorical, ordinal) ends in the table
#' @param non_factors a character vector containing the factors that are not considered as potential confounders to the factor tested for,
#'         for example, "Individual" and "SampleID" are random factors instead of confounders, so
#'         they should be listed as non_factors. Also make sure to include test_var, because it isn't a potential confounder in your test.
#' @param adjustMethod Multiple testing p value correction. Choices are the ones in p.adjust(), including
#'         "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr"). The default is "fdr".
#' @param model_p The threshold for significance of model test after multiple testing correction.
#'                The default is 0.1.
#' @param posthoc_p The threshold for significance of post-hoc test after multiple testing correction.
#'                The default is 0.05.
#' @export
#' @examples longdat_disc(table, 24, 5, c("Case", "Individual", "SampleID"))
#' @import lme4
#' @import orddom
#' @import reshape2
#' @import lmtest
#' @import tidyr
#' @import dplyr
#' @import stringr

longdat_disc <- function(input, test_var, value_column, factor_column, non_factors,
                         adjustMethod = "fdr", model_p = 0.1, posthoc_p = 0.05) {
  if (missing(input)) {
    stop('Error! Necessary argument "input" missing.')
  }
  if (missing(test_var)) {
    stop('Error! Necessary argument "test_var" missing.')
  }
  if (missing(value_column)) {
    stop('Error! Necessary argument "value_column" missing.')
  }
  if (missing(factor_column)) {
    stop('Error! Necessary argument "factor_column" missing.')
  }
  if (missing(non_factors)) {
    stop('Error! Necessary argument "non_factors" missing.')
  }

  data <- read.table (file = input, header = T, sep = "\t", check.names = F, row.names = NULL)
  # Here value column is defined by the user
  predictor_names <- (colnames(data))[1: (value_column - 1)]
  melt_data <- reshape::melt (data, id = predictor_names)
  # Omit the rows whose value column equals to NA
  melt_data <- melt_data %>% tidyr::drop_na(value)
  # Remove all dots in the bacteria name or it will cause problem
  melt_data$variable <- gsub(".", " ", melt_data$variable, fixed = TRUE)

  # Here factor column is defined the user
  # Make sure that all the columns are in the right class
  for (i in c((factor_column +1):((ncol(melt_data))-2))) {
    melt_data[ ,i] <- as.numeric(melt_data[ ,i])
  }
  for (i in c(1:factor_column)) {
    melt_data[ ,i] <- as.factor(melt_data[ ,i])
  }

  # Change the first column name of melt_data to "Individual"
  colnames(melt_data)[1] <- "Individual"

  # Variables are all the bacteria taxanomies
  variables <- unique (melt_data$variable)

  # Extract all the column names from melt_data except for the last two columns
  # which are variables and values
  factors <- colnames(melt_data)[-c(length(colnames(melt_data)),
                                    length(colnames(melt_data))-1)]

  ########## Calculate the p values for every factor (used for selecting factors later)
  N <- length (variables)
  Ps <- matrix(NA, N, ncol(melt_data)-2)
  rownames(Ps) <- variables
  colnames(Ps) <- factors
  col_NA <- c()

  for (i in 1:N) {   # loop through all variables
    aVariable = variables [i]
    print(i)
    print(aVariable)
    subdata <- subset(melt_data, variable == aVariable)

    for (j in 1:(ncol(subdata)-2)) {
      # If the factor has only two kinds of values
      if (length(unique(subdata[ , j])) == 2) {
        # sub1 selects all the data that has the first kind of value
        sub1 <- subdata[which(subdata[, j] == unique(subdata[ , j])[1]),]
        # sub2 selects all the data that has the second kind of value
        sub2 <- subdata[which(subdata[, j] == unique(subdata[ , j])[2]),]
        p <- wilcox.test(sub1$value, sub2$value, paired = F)$p.value
        Ps[i,j] <- p

        # If the factor has more than two kinds of values
      } else if (length(unique(subdata[ , j])) > 2) {
        fmla <- as.formula(paste("value ~ ", colnames(subdata)[j], sep = ""))
        p <- as.list(kruskal.test(fmla , data = subdata))$p.value
        Ps[i,j] <- p

        # If the factor has only one value
      } else if (length(unique(subdata[ , j])) < 2) {
        Ps[i,j] <- "NA"
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

  # Before extracting sel_fac, discard the factors defined by user "non_factors",
  # because we don't want them to be included in sel_fac
  Ps <- Ps[ , -match(non_factors, colnames(Ps))]

  # Extract the selected factors whose p value < 0.05
  sel_fac_ini <- c() # Make an empty list first
  for (i in 1:N) {  # loop through all variables
    facs <- c()
    for (j in 1:(ncol(Ps))) {
      if (Ps[i, j] < 0.05) {
        fac <- colnames(Ps)[j]
        facs <- c(facs, fac)
      }
    }
    sel_fac_ini[[i]] <- facs
  }
  names(sel_fac_ini) <- variables

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
        if (length(unique(subdata_individual[ , (sel_fac_ini[[i]][k])])) == 1) {
          f <- 0
        } else {
          f <- 1 }
        fac <- c(fac, f)
      }
      if (sum(fac) == 0) { # If all individuals have the same value on this factor
        ff <- NULL
      } else { # If there's any individual having different value on this factor
        ff <- sel_fac_ini[[i]][k]}
      facs <- c(facs, ff)
    }
    sel_fac[[i]] <- facs
  }
  names(sel_fac) <- variables
  print("Finished selecting factors.")

  ################## Model test ###############
  # Only focus on the comparison between "test_var" and other factors!!
  # The goal is to know if adding in "test_var" as a variable can enhance the fitness of
  # the model. Here, the "small model" is used, meaning that test_var is compared to other
  # factors separately. All q values <0.1 we can say "test_var" is sigificant for this bacteria
  Ps_model <- list() # Make a empty list first
  library(lme4)
  suppressWarnings(
    for (i in 1:N) {
      # loop through all variables
      aVariable = variables[i]
      print(i)
      subdata <- subset(melt_data, variable == aVariable)
      ps_lm <- c()
      if (length(sel_fac[[i]]) >= 1) {
        #tryCatch skips error in for loop
        tryCatch({
          for (k in 1:length(sel_fac[[i]])) {
            fmla1 <- as.formula(paste("rank(value) ~ (1| Individual) + (" , sel_fac[[i]][k], "| Individual)"))
            fmla2 <- as.formula(paste("rank(value) ~ (1| Individual) + (" , sel_fac[[i]][k], "| Individual)", "+", test_var))
            m1 <- lme4::glmer(data = subdata, fmla1, control=lmerControl(check.nobs.vs.nRE="ignore"))
            m2 <- lme4::glmer(data = subdata, fmla2, control=lmerControl(check.nobs.vs.nRE="ignore"))
            p <- lmtest::lrtest (m1, m2)
            # Here make sure that m2 likelihood is larger than m1, if m2 likelihood is smaller,
            # mark it as "1", meaning it's insignificant
            ifelse(test = p$LogLik[2] > p$LogLik[1],
                   yes = p_lm <- p$"Pr(>Chisq)"[2],
                   no = p_lm <- 1)
            names(p_lm) <- paste(sel_fac[[i]][k])
            ps_lm <- c(ps_lm, p_lm)
          }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      } else if (length(sel_fac[[i]]) == 0) {
        ps_lm <- "No significant factor for model test"
      }
      Ps_model[[i]] <- ps_lm
    }
  )
  names(Ps_model) <- gsub(".", " ", variables, fixed = TRUE)

  ##############Prepare to do FDR correction on Ps_model#################
  # Spread the Ps_model list
  unlist = unlist(Ps_model, use.names=FALSE) # output only the numbers
  mode(unlist) <- "numeric" # it's numeric now!
  unlist_noNA <-unlist[!is.na(unlist)] # Exclude NA in unlist
  unlist_withname = unlist(Ps_model) # numbers with names
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
  out <- stringr::str_split_fixed(ps_lm_m$name, pattern = fixed("."), n = 2)
  # Add the 2 new columns with split names to the matrix
  ps_lm_m2 <- cbind(out, ps_lm_m)
  # Remove the original unsplit name column, and then name the columns
  ps_lm_m2 $name = NULL
  colnames(ps_lm_m2) <- c("feature", "factor", "p_lm")
  # Reshaping the matrix so that can do fdr correction in the next step
  cast.ps_lm_m2 <- reshape2::dcast(data = ps_lm_m2, feature ~ factor)

  rownames(cast.ps_lm_m2) <- cast.ps_lm_m2$feature
  cast.ps_lm_m2 <- cast.ps_lm_m2[ , -1]
  cast.ps_lm_m2 <- cast.ps_lm_m2[ , -1]


  ################## FDR correction on Ps_model (model test) ###############
  # Do FDR correction (correct for the number of features)
  adjust_fun <- function(x) p.adjust(p = x, method = adjustMethod)
  p_lm_fdr <- apply(X = cast.ps_lm_m2, MARGIN = 2, FUN = adjust_fun)
  #Reorder the rows of p_lm_fdr to make it follow the order of variables
  p_lm_fdr <- p_lm_fdr[match(variables,rownames(p_lm_fdr)), ]

  Min <- function(x) min(x, na.rm = T)
  p_lm_fdr_min <- apply(p_lm_fdr, 1, FUN = Min)

  print("Finished model test.")

  ################# Inverse Model test #################
  Ps_inv_model <- list() # Make a empty list first
  library(lme4)
  suppressWarnings(
    for (i in 1:N) {
      # loop through all variables
      aVariable = variables[i]
      print(i)
      subdata <- subset(melt_data, variable == aVariable)
      ps_lm <- c()
      if (length(sel_fac[[i]]) >= 1) {
        tryCatch({
          for (k in 1:length(sel_fac[[i]])) {
            fmla1 <- as.formula(paste("rank(value) ~ ", "(1| Individual)", "+", test_var))
            fmla2 <- as.formula(paste("rank(value) ~ ", "(1| Individual) + (" , sel_fac[[i]][k], "| Individual)", "+", test_var))
            m1 <- lme4::glmer(data = subdata, fmla1, control=lmerControl(check.nobs.vs.nRE="ignore"))
            m2 <- lme4::glmer(data = subdata, fmla2, control=lmerControl(check.nobs.vs.nRE="ignore"))
            p <- lmtest::lrtest (m1, m2)
            # Here make sure that m2 likelihood is larger than m1, if m2 likelihood is smaller,
            # mark it as "1", meaning it's insignificant
            ifelse(test = p$LogLik[2] > p$LogLik[1],
                   yes = p_lm <- p$"Pr(>Chisq)"[2],
                   no = p_lm <- 1)
            names(p_lm) <- paste(sel_fac[[i]][k])
            ps_lm <- c(ps_lm, p_lm)
          }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      } else if (length(sel_fac[[i]]) == 0) {
        ps_lm <- "No significant factor for model test"
      }
      Ps_inv_model[[i]] <- ps_lm
    }
  )
  names(Ps_inv_model) <- gsub(".", " ", variables, fixed = TRUE)

  ############## Prepare to do FDR correction on Ps_inv_model #################
  # Spread the Ps_model list
  unlist_inv = unlist(Ps_inv_model, use.names=FALSE) # output only the numbers
  mode(unlist_inv) <- "numeric" # it's numeric now!
  unlist_inv_noNA <-unlist[!is.na(unlist_inv)] # Exclude NA in unlist
  unlist_inv_withname = unlist(Ps_inv_model) # numbers with names
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
  out_inv <- stringr::str_split_fixed(ps_inv_lm_m$name, pattern = fixed("."), n = 2)
  # Add the 2 new columns with split names to the matrix
  ps_inv_lm_m2 <- cbind(out_inv, ps_inv_lm_m)
  # Remove the original unsplit name column, and then name the columns
  ps_inv_lm_m2 $name = NULL
  colnames(ps_inv_lm_m2) <- c("feature", "factor", "p_inv_lm")
  # Reshaping the matrix so that can do fdr correction in the next step
  cast.ps_inv_lm_m2 <- reshape2::dcast(data = ps_inv_lm_m2, feature ~ factor)

  rownames(cast.ps_inv_lm_m2) <- cast.ps_inv_lm_m2$feature
  # Remove the first 2 columns that unwanted
  cast.ps_inv_lm_m2 <- cast.ps_inv_lm_m2[ , -1]
  cast.ps_inv_lm_m2 <- cast.ps_inv_lm_m2[ , -1]

  ################## FDR correction on inverse model ###############
  # Do FDR correction (correct for the number of features)
  adjust_fun <- function(x) p.adjust(p = x, method = adjustMethod)
  p_inv_lm_fdr <- apply(X = cast.ps_inv_lm_m2, MARGIN = 2, FUN = adjust_fun)
  #Reorder the rows of p_lm_fdr to make it follow the order of variables
  p_inv_lm_fdr <- p_inv_lm_fdr[match(variables,rownames(p_inv_lm_fdr)), ]

  print("Finished inverse model test.")

  ############## Post-hoc test for all variables #################
  case_pairs <- combn(x = sort(unique(melt_data[ , test_var])), m = 2)
  p_poho <- data.frame(matrix(nrow = length(variables), ncol = ncol(case_pairs)))
  case_pairs_name <- c()
  for (i in 1:N) { # loop through all variables
    print(i)
    bVariable = variables[i]
    subdata2 <- subset(melt_data, variable == bVariable)
    for (k in 1:ncol(case_pairs)) { # loop through each case pair
      sub3 <- subdata2[subdata2[ , test_var] == case_pairs[1,k], ]
      sub4 <- subdata2[subdata2[ , test_var] == case_pairs[2,k], ]
      # Here use "paired wilcoxon test because it's longitudinal data
      p_w <- wilcox.test(sub3$value, sub4$value, paired = T)$p.value
      p_poho[i, k] <- p_w
      name <- paste(case_pairs[1,k], sep = "_", case_pairs[2,k])
      case_pairs_name <- c(case_pairs_name, name)
    }
  }
  row.names(p_poho) <- variables
  colnames(p_poho) <- unique(paste("p_", case_pairs_name))
  library(tidyr)
  p_poho2t <- t(p_poho)

  ############## FDR correction on post-hoc test (on number of pairs) #################
  Ps_poho_fdrt <- apply(X = p_poho2t, MARGIN = 2, FUN = adjust_fun)

  if (ncol(case_pairs) == 1) {
    Ps_poho_fdrt <- as.data.frame(Ps_poho_fdrt)
    Ps_poho_fdr <- as.data.frame(Ps_poho_fdrt)
  } else {
    Ps_poho_fdr <- as.data.frame(t(Ps_poho_fdrt))
  }

  print("Finished post-hoc wilcoxon test.")

  ############## Calculate cliff's delta for all variables #################
  case_pairs <- combn(x = sort(unique(melt_data[ , test_var])), m = 2)
  delta <- data.frame(matrix(nrow = length(row.names(Ps_poho_fdr)), ncol = ncol(case_pairs)))
  case_pairs_name <- c()
  library(orddom)
  for (i in 1:length(row.names(Ps_poho_fdr))) { # loop through all variables
    print(i)
    bVariable = row.names(Ps_poho_fdr)[i]
    subdata2 <- subset(melt_data, variable == bVariable)
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

  print("Finished calculating effect size.")

  ############## Confirm confounding types #################
  # First determine the variables with post-hoc q(corrected p value) < 0.05,
  # then test whether inverse model q < 0.1 (confounded) or q > 0.1 (deconfounded)
  con <- c()
  Min <- function(x) min(x, na.rm = T)
  con <- apply(Ps_poho_fdr, 1, FUN = Min) < posthoc_p
  con_type <- c()
  for (i in 1:nrow(Ps_poho_fdr)) {
    if (con[i] == TRUE) { # post-hoc q is significant, meaning that other factors aren't confounding
      con_type[i] <- "Strictly deconfounded"
    } else if (con[i] == FALSE){
      if (min(p_inv_lm_fdr[i, ], na.rm = T) < model_p) { # post-hoc q isn't significant, while other factors are significant, so confounded
        con_type[i] <- "Confounded"
      }
      else if (min(p_inv_lm_fdr[i, ], na.rm = T) >= model_p){ # post-hoc and other factors are both insignificant
        con_type[i] <- "Ambiguously deconfounded"
      }
    }
  }

  ############## Generate potential confounding factors as output #################
  confound <- list() # Make an empty list first
  for (i in 1:nrow(p_inv_lm_fdr)) { # loop through all variables
    cc <- c()
    for (j in 1:ncol(p_inv_lm_fdr)) {
      if (!is.na(p_inv_lm_fdr[i, j]) && p_inv_lm_fdr[i, j] < model_p) {
        c <- colnames(p_inv_lm_fdr)[j]
      }
      else {
        c <- NA
      }
      cc <- c(cc, c)
    }
    confound[[i]] <- cc # listing out all the confounding factors for each variable
  }
  names(confound) <- variables

  # Remove "Pairlist of length 0" arrays from confound
  confound_noNA <- lapply(confound, function(x) x[!is.na(x)])
  confound_noNull <- Filter(length, confound_noNA)

  confound_noNull <- sapply(confound_noNull, "length<-", max(lengths(confound_noNull)))
  t_confound_noNull <- t(confound_noNull)

  write.table(x = t_confound_noNull, file = "confounders.txt", sep = "\t",
              row.names = T, col.names = NA, quote = F)

  ############## Generate result table as output #################
  result_table <- data.frame(matrix(nrow = length(variables), ncol = ncol(case_pairs)))
  row.names(result_table) <- variables
  colnames(result_table) <- paste("effect", sep = "_", unique(case_pairs_name))

  # Model q < 0.1 and Ps_poho_fdr < 0.05 and delta < 0 is depleted
  # Model q < 0.1 and Ps_poho_fdr < 0.05 and delta > 0 is enriched
  # Ps_poho_fdr > 0.05 is non-significant
  for (i in 1:nrow(Ps_poho_fdr)) {
    for (j in 1:ncol(Ps_poho_fdr)) {
      if (p_lm_fdr_min[i] < model_p | p_lm_fdr_min[i] == "Inf") {
        ifelse(test = Ps_poho_fdr[i, j] < posthoc_p, no = result_table[i, j] <- "NS",
               yes = ifelse(test = delta[i, j] < 0, yes = result_table[i, j] <- "Depleted",
                            no = result_table[i, j] <- "Enriched"))
      } else if (p_lm_fdr_min[i] > model_p) {
        result_table[i, j] <- "NS"
      }}}

  # Bind the columns together
  result_table <- cbind(result_table, Ps_poho_fdr, delta)

  # The final determination of this feature's signal
  sig_det <- function(x) {
    ifelse(test = sum(x == "NS") < ncol(case_pairs),
           yes = sig_det <- "OK", no =  sig_det <- "Insignificant")
  }
  result_table[ , "Signal"] <- apply(result_table, 1, FUN = sig_det)
  # Replace "NA" in signal column with "insignificant"
  result_table[is.na(result_table$Signal), "Signal"] <- "Insignificant"

  case_pairs <- combn(x = sort(unique(melt_data[ , test_var])), m = 2)
  case_pairs_name <- c()
  for (i in 1:ncol(case_pairs)) {
    pn <- paste(case_pairs[1, i], case_pairs[2, i], sep = "")
    case_pairs_name <- c(case_pairs_name, pn)
  }

  result_table <- cbind(result_table, con_type)
  # Change the signal of ambiguously deconfounded features from "Insignificant" to "Potential", meaning that
  # "test_var" is potentially significant for this feature
  result_table[result_table$con_type == "Ambiguously deconfounded", "Signal"] <- "Potential"

  colnames(result_table) <- c(paste("Effect_", sep = "", case_pairs_name), paste("Post-hoc_", sep = "", case_pairs_name),
                              paste("Effectsize_", sep = "", case_pairs_name), "Signal", "Confounding_type")

  write.table(x = result_table, file = "result_table.txt", sep = "\t",
              row.names = T, col.names = NA, quote = F)
  print("Result_table.txt is in your directory.")
  print("Finished! The result table and confounding table are now in your directory.")

}

