# !! Tidyr, dplyr is my dependency!!
# User has to define: value_column, factor_column, paired_test
# functionu 一開始也列出adjustMethod = "fdr",
# robustCutoff = 5, QCutoff = 0.1, DCutoff = 0等等預設值，也可以給人更改
# if (missing(metaMat)) {
#stop('Error - Necessary argument "metaMat" missing.')
#}
#if (missing(featureMat)) {
#  stop('Error - Necessary argument "featureMat" missing.')
#}
#if (nrow(metaMat) != nrow(featureMat)) {
#  stop("featureMat and metaMat don't have same number of rows.")
#}

rm(list=ls())

setwd("/Users/Jessica/Documents/Lab/CORONA")

data <- read.table (file = "corona.motu2matrix.rarefied.r.Phylum.r.master.r.dose_replaced.r",
                    header = T, sep = "\t", check.names = F, row.names = NULL)

# Ask users to put all microbes at the end of the dataset,
# and define the first microbe column
value_column <- 23
predictor_names <- (colnames(data))[1: (value_column - 1)]
melt_data <- reshape::melt (data, id = predictor_names)
# Omit the rows whose value column equals to NA
melt_data %>% tidyr::drop_na(value)

# Ask users to place ordinal data in the front of data, and define the factor_column
factor_column <- 5
# Make sure that all the columns are in the right class
for (i in c((factor_column +1):((ncol(melt_data))-2))) {
  melt_data[ ,i] <- as.numeric(melt_data[ ,i])
}
for (i in c(1:factor_column)) {
  melt_data[ ,i] <- as.factor(melt_data[ ,i])
}

# Variables are all the bacteria taxanomies
variables <- unique (melt_data$variable)

# Extract all the column names from melt_data except for the last two columns
# which are variables and values
factors <- colnames(melt_data)[-c(length(colnames(melt_data)),
                                  length(colnames(melt_data))-1)]
N <- length (variables)
Ps <- matrix(NA, N, ncol(melt_data)-2)
rownames(Ps) <- variables
colnames(Ps) <- factors

col_NA <- c()

for (i in 1:N) {
  # loop through all variables
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
    } else if (length(unique(melt_data[ , j])) >= 2) {
      fmla <- as.formula(paste("value ~ ", colnames(subdata)[j], sep = ""))
      p <- as.list(kruskal.test(fmla , data = subdata))$p.value
      Ps[i,j] <- p

      # If the factor has only one value
    } else if (length(unique(melt_data[ , j])) < 2) {
      Ps[i,j] <- "NA"
      col_NA <- c(col_NA, j)
    }
  }
}

# Save the original Ps matrix to Ps_ori
Ps_ori <- Ps

# Output the original P value table
#write.table(x = Ps_ori, file = "P_value_all.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)
# Exclude columns containing NA, and then turn the class of matrix into numeric
Ps <- Ps[ , -unique(col_NA)]
mode(Ps) <- "numeric"

# Do FDR correction (correct for the number of features)
adjust_fun <- function(x) p.adjust(p = x, method = "fdr", n = N)
Ps_fdr <- apply(X = Ps, MARGIN = 2, FUN = adjust_fun)


# Make an empty list
sel_fac <- c()
for (i in 1:N) {
  # loop through all variables
  facs <- c()
  for (j in 1:(ncol(Ps_fdr))) {
    if (Ps_fdr[i, j] < 0.1) {
      fac <- colnames(Ps_fdr)[j]
      facs <- c(facs, fac)
    }
  }
  sel_fac[[i]] <- facs
}
names(sel_fac) <- variables

