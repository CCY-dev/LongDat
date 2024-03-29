#' Calculate the p values for every factor (used for selecting factors later)
#' @param melt_data Internal function argument.
#' @param variables Internal function argument.
#' @param factor_columns Internal function argument.
#' @param factors Internal function argument.
#' @param data Internal function argument.
#' @param N Internal function argument.
#' @param verbose Internal function argument.
#' @importFrom rlang .data
#' @importFrom stats as.formula confint cor.test kruskal.test
#'             na.omit p.adjust wilcox.test
#' @importFrom magrittr '%>%'
#' @importFrom effsize cliff.delta
#' @import tidyr
#' @import rstatix
#' @name factor_p_cal

factor_p_cal <- function(melt_data, variables, factor_columns,
                         factors, data, N, verbose) {
  if (length(factor_columns) != 0) {
    Ps <- matrix(data = NA, nrow = N, ncol = length(factors))
    Ps_effectsize <- matrix(data = NA, nrow = N, ncol = length(factors))
    rownames(Ps) <- variables
    rownames(Ps_effectsize) <- variables
    colnames(Ps) <- factors
    colnames(Ps_effectsize) <- factors
    col_NA <- c()

    for (i in 1:N) {# loop through all variables
      aVariable <- variables[i]
      if (verbose == TRUE) {print(i)}
      subdata <- subset(melt_data, variable == aVariable)
      colnames(subdata) <- fix_name_fun(colnames(subdata))

      for (j in seq_len(length(factor_columns))) {
        # loop through all factor columns
        # Remove NA rows in subdata[ ,factor_columns[j]] first,
        # or else it will cause error
        subdata <- subdata %>%
          tidyr::drop_na(factor_columns[j])

        # If the factor has only two kinds of values, run wilcoxon test
        if (length(unique(subdata[ , factor_columns[j]])) == 2) {
          # sub1 selects all the data that has the first kind of value
          sub1 <- subdata[which(subdata[, factor_columns[j]] ==
                                  unique(subdata[ , factor_columns[j]])[1]),]
          # sub2 selects all the data that has the second kind of value
          sub2 <- subdata[which(subdata[, factor_columns[j]] ==
                                  unique(subdata[ , factor_columns[j]])[2]),]
          p <- wilcox.test(sub1$value, sub2$value, paired = FALSE)$p.value
          Ps[i,j] <- p
          #d <- as.numeric(dmes(x = sub1$value, y = sub2$value)$dc)
          # Orddom is outdated. Replace it with effsize
          # The order is treatment and then control in effsize::cliff.delta()
          d <- as.numeric(effsize::cliff.delta(sub2$value, sub1$value)$estimate)
          Ps_effectsize[i, j] <- d

          #If the factor has more than two kinds of values and is continuous
          # (it's class is "numeric"), run spearman's correlation test
        } else if (length(unique(subdata[ , factor_columns[j]])) > 2 &
                   is.numeric(subdata[ , factor_columns[j]])) {
          p <- stats::cor.test(subdata[ , factor_columns[j]],
                               subdata$value, method = "spearman")$p.value
          Ps[i,j] <- p
          d <- stats::cor.test(subdata[ , factor_columns[j]],
                               subdata$value, method = "spearman")$estimate
          Ps_effectsize[i, j] <- d

          # If the factor has more than two kinds of values and
          # it's class is "factor" (not numeric), run kruskal-wallis test
        } else if (length(unique(subdata[ , factor_columns[j]])) > 2 &
                   is.factor(subdata[ , factor_columns[j]])) {
          fmla <- as.formula(
            paste("value ~ ", colnames(subdata)[factor_columns[j]], sep = ""))
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

    # Exclude columns containing NA,
    # and then turn the class of matrix into numeric
    if (length(col_NA) > 0) {
      Ps <- Ps %>%
        as.data.frame() %>%
        dplyr::select(-c(unique(col_NA))) %>%
        as.matrix()
    }
    mode(Ps) <- "numeric"

    if (length(col_NA) > 0) {
      Ps_effectsize <- Ps_effectsize %>%
        as.data.frame() %>%
        dplyr::select(-c(unique(col_NA))) %>%
        as.matrix()
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

    # Here remove the factors that are fixed to the individual
    # (eg. age, sex...)
    sel_fac <- list()
    for (i in 1:N) { # loop through all variables
      aVariable <- variables[i]
      if (verbose == TRUE) {print(i)}
      subdata <- subset(melt_data, variable == aVariable)
      facs <- c()
      for (k in seq_len(length(sel_fac_ini[[i]]))) {
        # Loop through all sel_fac
        fac <- c()
        for (m in seq_len(length(unique(subdata$Individual)))) {
          # Loop through all idividuals
          subdata_individual <- subset(subdata, Individual ==
                                         (unique(subdata$Individual)[m]))
          if (!is.logical(sel_fac_ini[[i]])) {
            if (length(unique(
              subdata_individual[ , (sel_fac_ini[[i]][k])])) == 1) {
              f <- 0 }
            else {
              f <- 1 }
          } else if (is.logical(sel_fac_ini[[i]])){
            f <- 0
          }
          fac <- c(fac, f)
        }
        if (sum(fac) == 0) { # If all individuals have the same
          #value on this factor
          ff <- NA
        } else { # If there's any individual having different value
          #on this factor
          ff <- sel_fac_ini[[i]][k]}
        facs <- c(facs, ff)
      }
      sel_fac[[i]] <- facs
      if(is.null(facs)) {
        sel_fac[[i]] <- NA
      }
    }
    names(sel_fac) <- variables
    sel_fac <- lapply(sel_fac, function(x) x[!is.na(x)])
  }
  return(list(Ps, Ps_effectsize, sel_fac))
}
