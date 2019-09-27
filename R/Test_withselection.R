# !! Tidyr, dplyr is my dependency!!
# User has to define: value_column, factor_column, paired_test
# function 一開始也列出adjustMethod = "fdr",
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

data <- read.table (file = "corona.motu2matrix.rarefied.r.Genus.r.master.r.dose_replaced.r",
                    header = T, sep = "\t", check.names = F, row.names = NULL)

# Removing "." in the column names, or will cause error later!
#colnames(data) <- gsub(".", " ", colnames(data), fixed = TRUE)

# Ask users to put all microbes at the end of the dataset,
# and define the first microbe column
value_column <- 23
predictor_names <- (colnames(data))[1: (value_column - 1)]
melt_data <- reshape::melt (data, id = predictor_names)
# Omit the rows whose value column equals to NA
melt_data %>% tidyr::drop_na(value)

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
    } else if (length(unique(melt_data[ , j])) > 2) {
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

Ps_ori <- Ps # Save the original Ps matrix to Ps_ori

# Output the original P value table
#write.table(x = Ps_ori, file = "P_value_all.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
# Exclude columns containing NA, and then turn the class of matrix into numeric
Ps <- Ps[ , -unique(col_NA)]
mode(Ps) <- "numeric"

# Do FDR correction (correct for the number of features)
#adjust_fun <- function(x) p.adjust(p = x, method = "fdr", n = N)
#Ps_fdr <- apply(X = Ps, MARGIN = 2, FUN = adjust_fun)

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

# Then do the model test (lrtest)
Ps_model <- c() # Make a empty list first

for (i in 1:N) {
  # loop through all variables
  aVariable = variables[i]
  subdata <- subset(melt_data, variable == aVariable)
  ps_lm <- c()
  if (length(sel_fac[[i]]) >= 2) {
    pairs <- combn(x = sel_fac[[i]], m = 2)
    for (k in 1:length(sel_fac[[i]])) {
      for (m in 1:ncol(pairs)) {
        if (sel_fac[[i]][k] %in% pairs[ , m] ) {
          # setdiff(pairs[ , m], sel_fac[[i]][k]) finds the element of pairs[ , m] in sel_fac[[1]][1]
          fmla1 <- as.formula(paste("rank(value) ~ ", paste((setdiff(pairs[ , m], sel_fac[[i]][k])), collapse= "+")))
          fmla2 <- as.formula(paste("rank(value) ~ ", paste(pairs[ , m], collapse= "+")))
          m1 <- lm(data = subdata, fmla1)
          m2 <- lm(data = subdata, fmla2)
          p_lm <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)"[2]
          names(p_lm) <- paste(sel_fac[[i]][k])
          ps_lm <- c(ps_lm, p_lm)
        }
      }
    }
  } else if (length(sel_fac[[i]]) == 1){
    fmla3 <- as.formula(paste("value ~ ", colnames(subdata), sep = ""))
    ps_lm <- as.list(kruskal.test(fmla3 , data = subdata))$p.value
    names(ps_lm) <- paste(sel_fac[[i]])
  }
  else if (length(sel_fac[[i]]) == 0) {
    ps_lm <- "No significant factor for model test"
  }
  Ps_model[[i]] <- ps_lm
}
# Removing "." in the feature names, or will cause error later!
names(Ps_model) <- gsub(".", " ", variables, fixed = TRUE)


# Spread the Ps_model list
unlist = unlist(Ps_model, use.names=FALSE) # output only the numbers
mode(unlist) <- "numeric" # it's numeric now!
unlist_noNA <-unlist[!is.na(unlist)] # Exclude NA in unlist
unlist_withname = unlist(Ps_model) # numbers with names
unlist_name <- names(unlist_withname) # outputs the name of every element
unlist_name_noNA <- unlist_name[-which(is.na(unlist))] # remove names of NA value


max_ps_lm <- c() # Make an empty list first
for (z in 1:length(unique(unlist_name_noNA))) { # Loop through all unique unlist_name_noNA,
  #to find the largest p value of each factor in each feature(bacteria)
 m <- max(unlist_noNA[which(unlist_name_noNA == unique(unlist_name_noNA)[z])])
 max_ps_lm[z] <- m
}
names(max_ps_lm) <- unique(unlist_name_noNA) # name the numbers

# Create a empty matrix first
max_ps_lm_m <- data.frame(matrix(nrow = length(max_ps_lm), ncol = 2))
# Specify column names and fill in the values
colnames(max_ps_lm_m) <- c("name", "max_p_lm")
# Fill in the values
max_ps_lm_m[ , 1] <- unique(unlist_name_noNA)
max_ps_lm_m[ , 2] <- max_ps_lm

# Split the names of each p_lm, ex: Actinobacteria.PPI becomes "Actinobacteria" "PPI"
out <- stringr::str_split_fixed(max_ps_lm_m$name, pattern = fixed("."), 2)
# Add the 2 new columns with split names to the matrix
max_ps_lm_m2 <- cbind(out, max_ps_lm_m)
# Remove the original unsplit name column, and then name the columns
max_ps_lm_m2 $name = NULL
colnames(max_ps_lm_m2) <- c("feature", "factor", "max_p_lm")
# Reshaping the matrix so that can do fdr correction in the next step
cast.max_ps_lm_m2 <- reshape2::dcast(data = max_ps_lm_m2, feature ~ factor)

rownames(cast.max_ps_lm_m2) <- cast.max_ps_lm_m2$feature
cast.max_ps_lm_m2 <- cast.max_ps_lm_m2[ , -1]

# Do FDR correction (correct for the number of features)
adjust_fun <- function(x) p.adjust(p = x, method = "fdr", n = N)
p_lm_fdr <- apply(X = cast.max_ps_lm_m2, MARGIN = 2, FUN = adjust_fun)

# Select the significant factors of model test
sel_fac2 <- list()
for (i in 1:nrow(p_lm_fdr)) { # loop through all variables
  facs2 <- c(names(which(p_lm_fdr[i, ]< 0.1)))
  if (length(facs2) > 0) {
    sel_fac2[[i]] <- facs2
  } else {
    sel_fac2[[i]] <- NA
  }
}

sel_fac2 <- as.list(sel_fac2)
names(sel_fac2) <- row.names(p_lm_fdr)

