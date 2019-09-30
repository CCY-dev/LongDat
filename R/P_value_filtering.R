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
library(reshape)
melt_data <- reshape::melt (data, id = predictor_names)
# Omit the rows whose value column equals to NA
library(tidyr)
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
    } else if (length(unique(sub[ , j])) >= 2) {
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

# Save the original Ps matrix to Ps_ori
Ps_ori <- Ps

# Output the original P value table
write.table(x = Ps_ori, file = "P_value_all.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
# Exclude columns containing NA, and then turn the class of matrix into numeric
Ps <- Ps[ , -unique(col_NA)]
mode(Ps) <- "numeric"


# Do FDR correction (correct for the number of features)
adjust_fun <- function(x) p.adjust(p = x, method = "fdr", n = N)
Ps_fdr <- apply(X = Ps, MARGIN = 2, FUN = adjust_fun)


# Select only the factors with at least one p value smaller than 0.05
sel_fac <- colnames(Ps)[apply(Ps, 2, function(x) sum(x < 0.05) > 0)]
nonsel_fac <- colnames(Ps)[apply(Ps, 2, function(x) sum(x < 0.05) == 0)]


# Output the factors sel_fac and nonsel_fac
ff = (list("selected factors" = sel_fac, "non-selected factors" = nonsel_fac))
FF <- sapply(ff, "length<-", max(lengths(ff)))
write.file(FF, file = "Factor_filtering.txt", sep = "\t", quote = F)

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
    fmla1 <- as.formula(paste("rank(value) ~ ", paste(sel_fac_minus1, collapse= "+")))
    fmla2 <- as.formula(paste("rank(value )~ ", paste(sel_fac, collapse= "+")))
    m1 <- lm(data = subdata, fmla1)
    m2 <- lm(data = subdata, fmla2)
    # lrtest: Likelihood Ratio Test Of Nested Models
    p_lm <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]
    Ps_lm[i,j] <- p_lm
  }
}

# Output the Ps_lm table
write.table(x = Ps_lm, file = "P_value_lm.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

# Calculate cliff's delta
# First, ask user to define the column names that are paired (paired = T)
# in charactor vector
paired_test <- c("Case")

# Make empty lists first, because there's paired and unpaired cliff's delta, so
# here calculate them separately
D_unpaired <- c()
D_paired <- c()
N <- length (variables)

# Calculate cliff's delta only for the sel_fac
library(dplyr)
subdata_sel <- dplyr::select(subdata, -c(nonsel_fac))

# Generate the dataset for unpaired test
subdata_unpaired <- as.data.frame(subdata_sel[ , -match(x = paired_test,
                                                    table = colnames(subdata_sel))])

# And also generate the subset for paired test (only if the paired_test
#factors are present in sel_fac and the length of paired_test isn't zero)
if (paired_test %in% sel_fac & !length(paired_test) == 0) {
subdata_paired <- as.data.frame(subdata_sel[ , match(x = paired_test,
                                                     table = colnames(subdata_sel))])
colnames(subdata_paired) <- paired_test
subdata_paired[ , "variable"] <- subdata$variable
subdata_paired[ , "value"] <- subdata$value
}

for (i in 1:N) {
  # loop through all variables
  aVariable = variables [i]
  print(i)
  print(aVariable)
  # Make empty lists first, because there's paired and unpaired cliff's delta, so
  # here calculate them separately
  Delta_unpaired <- list()
  Delta_paired <- list()

  # Do the unpaired test first
  for (m in 1:(ncol(subdata_unpaired)-2)) {
    print(m)
    # Only calculate effect size for the ones with more than two values, and exclude
    # any factor whose levels has fewer than 3 values
    # 這一行是指factor中只要有一個少於3個值的就不會算，
    # 但是如果有一個factor有五個大於3個，只有一個小於3個，
    # 這個是否要計算？ （為什麼orddom會設最小值是3？）

    if (length(unique(subdata_unpaired[ , m])) >= 2
      & !any(as.numeric(table(subdata_unpaired[ , m]))<3)){
      # Split the dataset into groups by the column value
      subm_unpaired <- split(x = subdata_unpaired, f = subdata_unpaired[ , m])
      # Generate pairs to calculate pairwise effect size
      pairs_unpaired <- combn(x = seq(1:length(subm_unpaired)), m = 2)
      # D stores delta value
      D <- c()
      # N stores the name of the delta value
      Name <- c()
      for (i in 1:ncol(pairs_unpaired)) {
        print(i)
        # Pairs is a matrix, and we want to select subm by the numbers in pairs[1, i]
        # and pairs[2, i], and since subm is a list, so use [[]] to do selection
        library(orddom)
        d <- as.numeric (orddom::orddom(x = (subm_unpaired[[pairs_unpaired[1,i]]])$value,
                                      y = (subm_unpaired[[pairs_unpaired[2,i]]])$value, paired = F)[13,1])

        # This is a really clever way utilizing loop! All the D values will be filled after
        # looping
        D <- c(D, d)
        # n is the name of each d value
        n <- paste("d", names(subm_unpaired)[pairs_unpaired[1,i]], names(subm_unpaired)[pairs_unpaired[2,i]], sep = "_")
        # This is a really clever way utilizing loop! All the N values will be filled after
        # looping
        Name <- c(Name, n)
        }
      # Name the d values in D with the names in N
      names(D) <- Name
      # Collect all the D into Delta list
      Delta_unpaired[[m]] <- D
      } else {
        Delta_unpaired[[m]] <- "NA"
      }
    }
  # Name the lists
  names(Delta_unpaired) <- paste(aVariable, colnames(subdata_unpaired)[1:(ncol(subdata_unpaired)-2)]
                                 , sep = "_")
  #names(Delta_unpaired) <- colnames(subdata_unpaired)[1:(ncol(subdata_unpaired)-2)]
  D_unpaired <- c(D_unpaired,Delta_unpaired)

  # Then do the paired test (only if the paired_test factors are present in sel_fac
  # and the length of paired_test isn't zero)
  if (paired_test %in% sel_fac & !length(paired_test) == 0) {

    for (m in 1:(ncol(subdata_paired)-2)) {
      print(m)
      # Only calculate effect size for the ones with more than two values, and exclude
      # any factor whose levels has fewer than 3 values
      if (length(unique(subdata_paired[ , m])) >= 2
        & !any(as.numeric(table(subdata_paired[ , m]))<3)){

        # Split the dataset into groups by the column value
        subm_paired <- split(x = subdata_paired, f = subdata_paired[ , m])
        # Generate pairs to calculate pairwise effect size
        pairs_paired <- combn(x = seq(1:length(subm_paired)), m = 2)
        # D stores delta value
        D <- c()
        # N stores the name of the delta value
       Name <- c()
       for (i in 1:ncol(pairs_paired)) {
        print(i)
        # Pairs is a matrix, and we want to select subm by the numbers in pairs[1, i]
        # and pairs[2, i], and since subm is a list, so use [[]] to do selection
        d <- as.numeric (orddom::orddom(x = (subm_paired[[pairs_paired[1,i]]])$value,
                                      y = (subm_paired[[pairs_paired[2,i]]])$value, paired = T)[11,1])

        # This is a really clever way utilizing loop! All the D values will be filled after
        # looping
        D <- c(D, d)
        # n is the name of each d value
        n <- paste("d",names(subm_paired)[pairs_paired[1,i]], names(subm_paired)[pairs_paired[2,i]]
                   , sep = "_")
        # This is a really clever way utilizing loop! All the N values will be filled after
        # looping
        Name <- c(Name, n)
        }
      # Name the d values in D with the names in N
      names(D) <- Name
      # Collect all the D into Delta list
      Delta_paired[[m]] <- D
      } else {
        Delta_paired[[m]] <- "NA"
      }
    }
  }
  # Name the lists
  names(Delta_paired) <- paste(aVariable, colnames(subdata_paired)[1:(ncol(subdata_paired)-2)]
                               , sep = "_")
  D_paired <- c(D_paired,Delta_paired)
}

D_paired_df <- sapply(D_paired, "length<-", max(lengths(D_paired)))
D_unpaired_df <- sapply(D_unpaired, "length<-", max(lengths(D_unpaired)))

write.table(x = D_paired_df, file = "D_paired.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
write.table(x = D_unpaired_df, file = "D_unpaired.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

############################Testing############################

###### This is my script (without rank)
  aVariable = variables [16]
  print(aVariable)
  subdata <- subset(melt_data, variable == aVariable)
    sel_fac_minus1 <- sel_fac[-2]
    fmla1 <- as.formula(paste("value ~ ", paste(sel_fac_minus1, collapse= "+")))
    fmla2 <- as.formula(paste("value ~ ", paste(sel_fac, collapse= "+")))
    m1 <- lm(data = subdata, fmla1)
    m2 <- lm(data = subdata, fmla2)
    p_lm <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]
    p_lm

###### This is my script using rank
    fmla1 <- as.formula(paste("rank(value) ~ ", paste(sel_fac_minus1, collapse= "+")))
    fmla2 <- as.formula(paste("rank(value) ~ ", paste(sel_fac, collapse= "+")))
    m1 <- lm(data = subdata, fmla1)
    m2 <- lm(data = subdata, fmla2)
    library(lmtest)
    p_lm <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]
    p_lm
# !!!!!!!!!!! There is huge difference between with and without rank!!!!!!!!!!

###### From original CORONA script
    m1 <- lm (data = subdata, rank (value) ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    m2 <- lm (data = subdata, rank (value) ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    pModel <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]
    pModel

###### From original CORONA script excluding nonsel_fac
    m1 <- lm (data = subdata, rank (value) ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D )
    m2 <- lm (data = subdata, rank (value) ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D )
    pModel <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]
    pModel

