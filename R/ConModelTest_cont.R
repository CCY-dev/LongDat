#' Confounding model test in longdat_cont() pipeline
#' @param N Internal function argument.
#' @param variables Internal function argument.
#' @param melt_data Internal function argument.
#' @param sel_fac Internal function argument.
#' @param data_type Internal function argument.
#' @param test_var Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data
#' @importFrom stats as.formula confint cor.test kruskal.test na.omit p.adjust wilcox.test
#' @importFrom car Anova
#' @importFrom  MASS polr
#' @importFrom magrittr '%>%'
#' @import dplyr
#' @import bestNormalize
#' @import glmmTMB
#' @import lme4
#' @name ConModelTest_cont
utils::globalVariables(c("value"))

ConModelTest_cont <- function(N, variables, melt_data, sel_fac, data_type, test_var, verbose) {
  Ps_conf_model <- list()
  Ps_inv_conf_model <- list()
  suppressWarnings(
    for (i in 1:N) { # loop through all variables
      aVariable = variables[i]
      if (verbose == T) {print(i)}
      subdata <- subset(melt_data, variable == aVariable)
      colnames(subdata) <- fix_name_fun(colnames(subdata))
      tryCatch({
        if (data_type %in% c("measurement", "others")) {
          subdata <- subdata %>%
            dplyr::mutate(value_norm = bestNormalize::bestNormalize(value, loo = T)$x.t)
        }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      ps_lm <- c()
      ps_lm_inv <- c()
      if (length(sel_fac[[i]]) >= 1) {
        #tryCatch skips error in for loop
        tryCatch({
          for (k in 1:length(sel_fac[[i]])) { # Sel_fac exists
            if (data_type == "count") {
              fmla2 <- as.formula(paste("value ~ (1| Individual) +" , fix_name_fun(sel_fac[[i]][k]), "+", test_var))
              m2 <- glmmTMB::glmmTMB(formula = fmla2, data = subdata, family = nbinom2, na.action = na.omit, REML = F)
            } else if (data_type == "proportion") {
              fmla2 <- as.formula(paste("value ~ (1| Individual) +" , fix_name_fun(sel_fac[[i]][k]), "+", test_var))
              m2 <- glmmTMB::glmmTMB(fmla2, data = subdata, family = beta_family(), na.action = na.omit, REML = F)
            } else if (data_type %in% c("measurement", "others")) {
              fmla2 <- as.formula(paste("value_norm ~ (1| Individual) +" , fix_name_fun(sel_fac[[i]][k]), "+", test_var))
              m2 <- lme4::lmer(data = subdata, fmla2, REML = F)
            } else if (data_type == "binary") {
              fmla2 <- as.formula(paste("value ~ (1| Individual) +" , fix_name_fun(sel_fac[[i]][k]), "+", test_var))
              m2 <- glmmTMB::glmmTMB(fmla2, data = subdata, family = "binomial", na.action = na.omit, REML = F)
            } else if (data_type == "ordinal") {
              fmla2 <- as.formula(paste("as.factor(value) ~ (1| Individual) +" , fix_name_fun(sel_fac[[i]][k]), "+", test_var))
              m2 <- MASS::polr(fmla2, data = subdata, method = "logistic")
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
  names(Ps_conf_model) <- variables
  names(Ps_inv_conf_model) <- variables
  return(list(Ps_conf_model, Ps_inv_conf_model))
}
