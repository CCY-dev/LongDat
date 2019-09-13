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
for (i in c((factor_column +1):((ncol(melt_data))-2))) {
  melt_data[,i] <- as.numeric(melt_data[,i])
}
for (i in c(1:factor_column)) {
  melt_data[,i] <- as.factor(melt_data[,i])
}

variables <- unique (melt_data$variable)
factors <- colnames(melt_data)[-c(length(colnames(melt_data)),
                                  length(colnames(melt_data))-1)]
N <- length (variables)
Ps <- matrix(NA,N,ncol(melt_data)-2)
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
    if (length(unique(subdata[ , j])) == 2) {
    sub1c <- subdata[which(subdata[,j]==unique(subdata[ , j])[1]),]
    sub2c <- subdata[which(subdata[,j]==unique(subdata[ , j])[2]),]
    p <- wilcox.test(sub1c$value, sub2c$value, paired = F)$p.value
    Ps[i,j] <- p
    } else if (length(unique(melt_data[ , j])) >= 2) {
    fmla <- as.formula(paste("value ~ ", colnames(subdata)[j], sep = ""))
    p <- as.list(kruskal.test(fmla , data = subdata))$p.value
    Ps[i,j] <- p
    } else if (length(unique(melt_data[ , j])) < 2) {
    Ps[i,j] <- "NA"
    col_NA <- c(col_NA,j)
    }
    }
}

Ps_ori <- Ps

Ps <- Ps[,-unique(col_NA)]
mode(Ps) <- "numeric"

sel_fac <- colnames(Ps)[apply(Ps, 2, function(x) sum(x < 0.05) > 0)]
nosel_fac <- colnames(Ps)[apply(Ps, 2, function(x) sum(x < 0.05) == 0)]

#接下來寫for迴圈去做lmtest fmla2 <. as.formula(paste(sel_fac))

############################Testing############################

length(unique(melt_data[ , 1]))==2
unique(melt_data[ , 1])[1]
unique(melt_data[ , 1])[2]
sub1c <- melt_data[which(melt_data[,5]==unique(melt_data[ , 5])[1]),]
sub2c <- melt_data[which(melt_data[,5]==unique(melt_data[ , 5])[2]),]

which(melt_data[,5]==unique(melt_data[ , 5])[1])
sub1c$value
sub2c$value
pair <-
p12 <- wilcox.test(sub1c$value, sub2c$value, paired = F)$p.value
p12


aVariable = variables [1]
print(1)
print(aVariable)
variables <- unique (melt_data$variable)
N <- length (variables)

variables

subdata <- subset(melt_data, variable == aVariable)

fmla <- as.formula(paste("value ~ ", colnames(subdata)[5], sep = ""))
p12 <- as.list(kruskal.test(fmla, data = subdata))$p.value
value ~ colnames(subdata)[5]
colnames(subdata)[5]
fmla
p12
Ps <- data.frame (variable = variables)
Ps

factors <- colnames(melt_data)[-c(length(colnames(melt_data)),length(colnames(melt_data))-1)]
variables
class(aVariable)

subdata <- subset(melt_data, variable == aVariable)
Ps <- matrix(NA,N,ncol(subdata)-2)
rownames(Ps) <- variables
colnames(Ps) <- factors
Ps

matrix(c("a","b","c","d","e",1),2,3, byrow = T)

