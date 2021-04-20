#' Data preprocessing
#' @param input Internal function argument.
#' @param test_var Internal function argument.
#' @param variable_col Internal function argument.
#' @param fac_var Internal function argument.
#' @param not_used Internal function argument.

data_preprocess <- function(input, test_var, variable_col, fac_var, not_used) {
data <- read.table (file = input, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
# Remove the features (bacteria) whose column sum is 0
values <- data %>% dplyr::select(all_of(variable_col:ncol(data)))
values <- as.data.frame(apply(values, 2, as.numeric))
values <- as.data.frame(values[ , colSums(values,  na.rm = T) > 0])
colnames(values) <- colnames(data %>% dplyr::select(all_of(variable_col:ncol(data))))
data <- as.data.frame(cbind(data[ , 1:(variable_col-1)], values))
mean_abundance <- round(colSums(values, na.rm = T)/nrow(data), 3)
prevalence <- c()
for (i in 1:ncol(values)) {
  prevalence[i] <- round(sum(values[ , i] > 0, na.rm = T)/nrow(data) * 100, 3)
}

variables_original <- colnames(data)[variable_col:ncol(data)]

# Remove all the symbols from variables
colnames(data)[variable_col:ncol(data)] <- fix_name_fun(colnames(data)[variable_col:ncol(data)])

# Here value column is defined by the user
predictor_names <- (colnames(data))[1: (variable_col - 1)]
melt_data <- reshape2::melt (data, id = predictor_names)
# Omit the rows whose value column equals to NA
melt_data <- melt_data %>% tidyr::drop_na(.data$value)

# Remove the not-used columns in melt_data
if (!is.null(not_used)) {
  melt_data <- melt_data %>% dplyr::select(-c(not_used))
}

# Make sure that all the columns are in the right class
# Columns mentioned in fac_var, and the second last column in melt data are factors
# Columns not in fac_var, and the last column in melt data are numerical numbers
"%notin%" <- Negate("%in%")
num_var <- c(which(1:(ncol(melt_data)-2) %notin% fac_var), ncol(melt_data))
fac_var <- c(fac_var, ncol(melt_data)-1)
for (i in fac_var) {
  melt_data[ ,i] <- as.factor(melt_data[ ,i])
}
for (i in num_var) {
  melt_data[ ,i] <- as.numeric(as.character(melt_data[ ,i]))
}


# Change the first column name of melt_data to "Individual"
colnames(melt_data)[1] <- "Individual"

# Variables are all the bacteria taxanomies
variables <- unique(melt_data$variable)

# Extract all the column names from melt_data except for the last two columns
# which are variables and values, and also exclude "Individual" and test_var
factors <- colnames(melt_data)[-c(ncol(melt_data), ncol(melt_data)-1)]
factors <- factors[-which(factors %in% c("Individual", test_var))]
factor_columns <- match(factors, colnames(melt_data))
return(list(mean_abundance, prevalence, variables_original,
            melt_data, variables, factor_columns, factors, data, values))
}
