#' Calculate the p values for every factor (used for selecting factors later)

factor_p_cal <- function(melt_data, variables, factor_columns, factors, data, N, verbose) {
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
      if (verbose == T) {print(i)}
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
      if (verbose == T) {print(i)}
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
  }
  return(list(Ps, Ps_effectsize, sel_fac))
}
