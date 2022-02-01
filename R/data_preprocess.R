#' Data preprocessing
#' @param input Internal function argument.
#' @param test_var Internal function argument.
#' @param variable_col Internal function argument.
#' @param fac_var Internal function argument.
#' @param not_used Internal function argument.
#' @importFrom stats as.formula confint cor.test
#'             kruskal.test na.omit p.adjust wilcox.test
#' @importFrom magrittr '%>%'
#' @importFrom rlang .data
#' @import dplyr
#' @import reshape2
#' @import tidyr
#' @name data_preprocess

data_preprocess <- function(input, test_var, variable_col, fac_var, not_used) {
data <- input
values <- data %>% dplyr::select(all_of(variable_col:ncol(data)))
values <- as.data.frame(apply(values, 2, as.numeric))

# Count the number of unique values for each column
unique_values <- apply(values, 2, unique)
if (class(unique_values)[1] == "matrix") {
  unique_length <- apply(unique_values, 2, length)
} else {
  unique_length <- lapply(unique_values, length)
}

# Remove the features whose values are all 0
# Zero_ones are those who only have 0 as values (sum=0 & unique length=1)
zero_ones <- which(colSums(values,  na.rm = TRUE) == 0 &
                     unique_length == 1)
values <- values %>%
  dplyr::select(-all_of(zero_ones))
colnames(values) <- colnames(data %>%
                               dplyr::select(variable_col:ncol(data)) %>%
                               dplyr::select(-all_of(zero_ones)))
data <- as.data.frame(cbind(data[ , 1:(variable_col-1)], values))
mean_abundance <- round(colSums(values, na.rm = TRUE)/nrow(data), 3)
prevalence <- c()
for (i in seq_len(ncol(values))) {
  prevalence[i] <- round(sum(values[ , i] > 0,
                             na.rm = TRUE)/nrow(data) * 100, 3)
}

variables_original <- colnames(data)[variable_col:ncol(data)]

# Remove all the symbols from variables
colnames(data)[variable_col:ncol(data)] <-
  fix_name_fun(colnames(data)[variable_col:ncol(data)])

# Here value column is defined by the user
predictor_names <- (colnames(data))[1: (variable_col - 1)]
melt_data <- reshape2::melt (data, id = predictor_names)
# Omit the rows whose value column equals to NA
melt_data <- melt_data %>%
  tidyr::drop_na(.data$value)

# Make sure that all the columns are in the right class
# Columns mentioned in fac_var, and the second last column
# in melt data are factors
# Columns not in fac_var, and the last column in melt data are
# numerical numbers
"%notin%" <- Negate("%in%")
num_var <- c(which(1:(ncol(melt_data)-2) %notin% fac_var), ncol(melt_data))
fac_var <- c(fac_var, ncol(melt_data)-1)
for (i in fac_var) {
  melt_data[ ,i] <- as.factor(melt_data[ ,i])
}
for (i in num_var) {
  melt_data[ ,i] <- as.numeric(as.character(melt_data[ ,i]))
}

# Remove the not-used columns in melt_data
if (!is.null(not_used)) {
  melt_data <- melt_data %>%
    dplyr::select(-c(not_used))
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
