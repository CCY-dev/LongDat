rm(list=ls())

setwd("/Users/Jessica/Documents/Lab/CORONA")

data <- read.table (file = "corona.motu2matrix.rarefied.r.Phylum.r.master.r.dose_replaced.r",
                    header = T, sep = "\t", check.names = F, row.names = NULL)

# Ask users to put all microbes at the end of the dataset,
# and define the first microbe column
value_column <- 23
predictor_names <- (colnames(data))[1: (value_column - 1)]
melt_data <- reshape::melt (data, id = predictor_names)

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

# Exclude columns containing NA, and then turn the class of matrix into numeric
Ps <- Ps[ , -unique(col_NA)]
mode(Ps) <- "numeric"

# Select only the factors with at least one p value smaller than 0.05
sel_fac <- colnames(Ps)[apply(Ps, 2, function(x) sum(x < 0.05) > 0)]
nonsel_fac <- colnames(Ps)[apply(Ps, 2, function(x) sum(x < 0.05) == 0)]

N <- length (variables)
# Do lmtest with the selected factors
Ps_lm <- matrix(NA, N, length(sel_fac))
for (i in 1:N) {
  # loop through all variables
  aVariable = variables [i]
  print(i)
  print(aVariable)
  subdata <- subset(melt_data, variable == aVariable)

  rownames(Ps_lm) <- variables
  colnames(Ps_lm) <- sel_fac

  for (j in 1:length(sel_fac)) {
    sel_fac_minus1 <- sel_fac[-j]
    fmla1 <- as.formula(paste("value ~ ", paste(sel_fac_minus1, collapse= "+")))
    fmla2 <- as.formula(paste("value ~ ", paste(sel_fac, collapse= "+")))
    m1 <- lm(data = subdata, fmla1)
    m2 <- lm(data = subdata, fmla2)
    p_lm <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]
    Ps_lm[i,j] <- p_lm
  }
}

# 也是和上面一樣，分成2 groups (cliff's delta) and > 2 groups (?)
for (i in 1:N) {
  # loop through all variables
  aVariable = variables [i]
  print(i)
  print(aVariable)
N <- length (variables)
deltas <- matrix(NA, N, length(sel_fac))
rownames(deltas) <- variables
colnames(deltas) <- sel_fac



############################Testing############################

# Make an empty list first
Delta <- list()
for (m in 1:(ncol(subdata)-2)) {
  print(m)
  # Only calculate effect size for the ones with more than two values, and exclude
  # any factor whose levels has fewer than 3 values
  if (length(unique(subdata[ , m])) >= 2 & !any(as.numeric(table(subdata[ , m]))<3)){
  # Split the dataset into groups by the column value
  subm <- split(x = subdata, f = subdata[ , m])
  # Generate pairs to calculate pairwise effect size
  pairs <- combn(x = seq(1:length(subm)), m = 2)
  # D stores delta value
  D <- c()
  # N stores the name of the delta value
  N <- c()
  for (i in 1:ncol(pairs)) {
    print(i)
    # Pairs is a matrix, and we want to select subm by the numbers in pairs[1, i]
    # and pairs[2, i], and since subm is a list, so use [[]] to do selection
    d <- as.numeric (orddom::orddom(x = (subm[[pairs[1,i]]])$value,
                                      y = (subm[[pairs[2,i]]])$value, paired = F)[11,1])
    # This is a really clever way utilizing loop! All the D values will be filled after
    # looping
    D <- c(D, d)
    # n is the name of each d value
    n <- paste("d",names(subm)[pairs[1,i]],names(subm)[pairs[2,i]],sep = "_")
    # This is a really clever way utilizing loop! All the N values will be filled after
    # looping
    N <- c(N, n)
  }
  # Name the d values in D with the names in N
  names(D) <- N
  # Collect all the D into Delta list
  Delta[[m]] <- D
  } else {
  # Name the lists
  Delta[[m]] <- NA
  }
}
names(Delta) <- colnames(subdata)[1:(ncol(subdata)-2)]
Delta


!any(as.numeric(table(subdata[,2])) < 3)

length(which(subdata[,m] == unique(subdata[,9])[1])) > 3
length(which(subdata[,m] == unique(subdata[,9])[2])) > 3
unique(subdata[ , 9])
length(unique(subdata[ , 9])) >= 2
if (length((subm[[pairs[1,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value) <= 3 | length((subm[[pairs[2,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value) <= 3){
  print("Yes")
}
d <- as.numeric (orddom::orddom(x = (subm[[pairs[1,2]]] [! is.na (subm[[pairs[1,2]]]$value) & ! is.na (subm[[pairs[2,2]]]$value), ])$value,
                                y = (subm[[pairs[2,2]]] [! is.na (subm[[pairs[1,2]]]$value) & ! is.na (subm[[pairs[2,2]]]$value), ])$value, paired = F)[11,1])

(subm[[pairs[1,2]]] [! is.na (subm[[pairs[1,2]]]$value) & ! is.na (subm[[pairs[2,2]]]$value), ])$value

! is.na (subm[[pairs[1,2]]]$value) & ! is.na (subm[[pairs[2,2]]]$value)

subm[[pairs[1,1]]]
length(unique(subdata[ , 6])) >= 2
subm = split(x = subdata, f = subdata[ , 6])
pairs <- combn(x = seq(1:length(subm)), m = 2)
class(pairs)
subm[[pairs[1,1]]][! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ]
length((subm[[pairs[1,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value)
(subm[[pairs[2,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value

d <- as.numeric (orddom::orddom
                   (x = (subm[[pairs[1,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value,
                    y = (subm[[pairs[2,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value, paired = F)[11,1])

orddom::orddom(x = (subm[[pairs[1,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value, y = (subm[[pairs[2,1]]] [! is.na (subm[[pairs[1,1]]]$value) & ! is.na (subm[[pairs[2,1]]]$value), ])$value, paired = T)

orddom::orddom(x=c(2,0,0,1), y=c(3,1,5,0,9), paired = F)[11,1]

d12 <- as.numeric (orddom::orddom
                   (x = (sub1c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])
                     $value, y = (sub2c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])
                     $value, paired = T)[11,1])

# Then generate
for (n in length(subm)) {
  subm_separate <- paste0("subm_", n)
  subm_separate <- subm[[n]]
}
