# !! Tidyr, dplyr, reshape are my dependencies!!
rm(list=ls())

setwd("/Users/Jessica/Documents/Lab/CORONA")

data <- read.table (file = "corona.motu2matrix.rarefied.r.Species.r.master.r.dose_replaced.r",
                    header = T, sep = "\t", check.names = F, row.names = NULL)

# Ask users to put all microbes at the end of the dataset,
# and define the first microbe column
value_column <- 24
predictor_names <- (colnames(data))[1: (value_column - 1)]
melt_data <- reshape::melt (data, id = predictor_names)
# Omit the rows whose value column equals to NA
library(tidyr)
melt_data <- melt_data %>% tidyr::drop_na(value)
#melt_data$variable <- gsub(".", " ", melt_data$variable, fixed = TRUE)

# Ask users to put "Patient" at the first column!!!
# Ask users to place the columns that aren't numeric
# in the front of data, and define the factor_column
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

# Calculate the p values for every factor (used for selecting factors later)
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
Ps <- Ps[ , -unique(col_NA)]
mode(Ps) <- "numeric"

# Before extracting sel_fac, discard "patient" in sel_fac because it's a random factor
Ps <- Ps[ , -1]

# Extract the selected factors whose p value < 0.05
sel_fac <- c() # Make an empty list first
for (i in 1:N) {  # loop through all variables
  facs <- c()
  for (j in 1:(ncol(Ps))) {
    if (Ps[i, j] < 0.05) {
      fac <- colnames(Ps)[j]
      facs <- c(facs, fac)
    }
  }
  sel_fac[[i]] <- facs
}
names(sel_fac) <- variables

################## Model test ###############
# Only focus on the comparison between "Case" and other factors!!
Ps_model <- c() # Make a empty list first
for (i in 1:N) {
  # loop through all variables
  aVariable = variables[i]
  subdata <- subset(melt_data, variable == aVariable)
  ps_lm <- c()
  if (length(sel_fac[[i]]) >= 1) {
    for (k in 1:length(sel_fac[[i]])) {
          # setdiff(pairs[ , m], sel_fac[[i]][k]) finds the element of pairs[ , m] in sel_fac[[1]][1]
          fmla1 <- as.formula(paste("rank(value) ~ ", paste((setdiff(pairs[ , m], sel_fac[[i]][k])), collapse= "+")))
          fmla2 <- as.formula(paste("rank(value) ~ ", paste(pairs[ , m], collapse= "+")))
          m1 <- lm(data = subdata, fmla1)
          m2 <- lm(data = subdata, fmla2)
          p_lm <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)"[2]
          names(p_lm) <- paste(sel_fac[[i]][k])
          ps_lm <- c(ps_lm, p_lm)
    }
  } else if (length(sel_fac[[i]]) == 0) {
    ps_lm <- "No significant factor for model test"
  }
  Ps_model[[i]] <- ps_lm
}



