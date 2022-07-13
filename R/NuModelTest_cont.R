#' Null Model Test and post-hoc Test in longdat_cont() pipeline
#' @param N Internal function argument.
#' @param data_type Internal function argument.
#' @param test_var Internal function argument.
#' @param melt_data Internal function argument.
#' @param variables Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data
#' @importFrom  MASS polr
#' @importFrom car Anova
#' @import dplyr
#' @import glmmTMB
#' @import lme4
#' @import tibble
#' @import bestNormalize
#' @importFrom stats as.formula confint cor.test
#'             kruskal.test na.omit p.adjust wilcox.test
#' @importFrom magrittr '%>%'
#' @name NuModelTest_cont
utils::globalVariables(c("value"))

NuModelTest_cont <- function(N, data_type, test_var, melt_data,
                             variables, verbose) {
  Ps_null_model <- as.data.frame(matrix(data = NA, nrow = N, ncol = 2))
  if (data_type == "count") {
    Theta <- c()
  }
  if (data_type %in% c("measurement", "others")) {
    Normalize_method <- as.data.frame(matrix(data = NA, nrow = N, ncol = 2))
    colnames(Normalize_method) <- c("Feature", "Normalization_method")
    Normalize_method$Feature <- variables
  }
  suppressWarnings(
    for (i in 1:N) { # loop through all variables
      aVariable <- variables[i]
       if (verbose == TRUE) {print(i)}
      subdata <- subset(melt_data, variable == aVariable)

      tryCatch({
        if (data_type %in% c("measurement", "others")) {
          normalize_result <- bestNormalize::bestNormalize(subdata$value,
                                                           loo = TRUE)
          # Record the method of normalization
          # Find the one with the lowest Pearson P/df, lower = more normal
          Normalize_method$Normalization_method[i] <-
            names(normalize_result$norm_stats)[
              which.min(normalize_result$norm_stats)]

          subdata <- subdata %>%
            mutate(value_norm = normalize_result$x.t)
          remove(normalize_result)
        }
        if (data_type == "count") {
          # Negative binomial
          fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
          m2 <- glmmTMB::glmmTMB(formula = fmla2, data = subdata,
                                 family = nbinom2, na.action = na.omit,
                                 REML = FALSE)
          # Extract dispersion theta out of model
          Theta[i] <- glmmTMB::sigma(m2)
        } else if (data_type == "proportion") {
          fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
          m2 <- glmmTMB::glmmTMB(fmla2, data = subdata, family = beta_family(),
                                 na.action = na.omit, REML = FALSE)
        } else if (data_type %in% c("measurement", "others")) {
          fmla2 <- as.formula(paste("value_norm ~ (1|Individual) +", test_var))
          m2 <- lme4::lmer(data = subdata, fmla2, REML = FALSE)
        } else if (data_type == "binary") {
          fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
          m2 <- glmmTMB::glmmTMB(fmla2, data = subdata, family = "binomial",
                                 na.action = na.omit, REML = FALSE)
        } else if (data_type == "ordinal") {
          fmla2 <- as.formula(paste("as.factor(value) ~ (1| Individual) +",
                                    test_var))
          m2 <- MASS::polr(fmla2, data = subdata, method = "logistic")
        }

        # Wald Chisq test
        Ps_null_model[i, 1] <-
          car::Anova(m2, type=c("II"),  test.statistic=c("Chisq"))$"Pr(>Chisq)"

        # Calculate confidence interval for test_var
        # Followed by the determination of the signs. For "- -" or "+ +",
        # it means that the doesn't span 0.
        # But "- +" means that the CI spans 0.
        # Here I sum the signs over the row, sum != 0 means it's OK.
        remove(ci)
        ci <- as.data.frame(confint(m2))
        ci <- ci[str_detect(row.names(ci), test_var), 1:2]
        if (nrow(ci) > 1) {
          if (any(apply(sign(ci), 1, sum, na.rm = TRUE) != 0)) {
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
      tibble::rownames_to_column() %>%
      dplyr::mutate(NB_theta = Theta) %>%
      tibble::column_to_rownames()
  }
  if (data_type %in% c("measurement", "others")) {
    return(list(Ps_null_model, Normalize_method))
  } else  {
    return(Ps_null_model)
  }
}
