#' Plot theta values of negative binomial models versus non-zero count for count data
#' @param input A character vector. This is the path to a txt file with the first column as "Individual", and all the dependent variables (ex: bacteria)
#'         should be at the end of the table.
#' @param test_var The name of the independent variable you are testing for, should be a character vector (ex: c("Time"))
#'        identical to its column name and make sure there is no space in it.
#' @param variable_col The column number of the position where the dependent variable columns (ex: bacteria) start in the table
#' @param fac_var The column numbers of the position where the columns that aren't
#'         numerical  (e.g. characters, categorical numbers, ordinal numbers), should be a numerical vector (ex: c(1, 2, 5:7))
#' @param not_used The column position of the columns not are irrevelant and can be ignored when in the analysis.
#'        This should be a number vector, and the default is NULL.
#' @param point_size The point size for plotting in ggplot2. The default is 1.
#' @param x_interval_value The interval value for tick marks on x-axis. The default is 5.
#' @param y_interval_value The interval value for tick marks on y-axis. The default is 5.
#' @param verbose A boolean vector indicating whether to print detailed message. The default is T.
#' @export
#' @import tidyr
#' @import reshape2
#' @import glmmTMB
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @importFrom stats as.formula confint cor.test kruskal.test na.omit p.adjust wilcox.test
#' @importFrom rlang .data
#' @importFrom magrittr '%>%'
#' @name theta_plot
#' @return a ggplot object
#' @details
#' This function outputs a plot that facilitates the setting of theta_cutoff in longdat_disc()
#' and longdat_cont(). This only applies when the dependent variables are count data. Longdat_disc()
#' and longdat_cont() implements negative binomial (NB) model for count data, and if the theta (dispersion parameter) of
#' NB model gets too high, then the p value of it will be extremely low regardless of whether there is real significance
#' or not. Therefore, the highest threshold of theta value is set and any variable beyond the threshold will be excluded
#' from the test. The default value of theta_cutoff is set to 2^20 from the observation that 2^20 is a clear cutoff line
#' for several datasets. Users can change theta_cutoff value to fit their own data.
#' The "nonzero_count_vs_theta.pdf" will be in the output directory.
#' @examples
#'\dontrun{
#' # Get the path of example dataset
#' system.file("Fasting_disc.txt", package = "longdat")
#' # Paste the directory to the input below
#' test_theta <- theta_plot(input = "your_path_to/Fasting_disc.txt", test_var = "Time_point",
#'            variable_col = 7, fac_var = c(1:3))
#'}
utils::globalVariables(c("values", "NB_theta", "Nonzero_count"))

theta_plot <- function(input, test_var, variable_col, fac_var, not_used = NULL,
                      point_size = 1, x_interval_value = 5,
                       y_interval_value = 5, verbose = T) {
  if (missing(input)) {
    stop('Error! Necessary argument "input" is missing.')
  }
  if (missing(test_var)) {
    stop('Error! Necessary argument "test_var" is missing.')
  }
  if (missing(variable_col)) {
    stop('Error! Necessary argument "variable_col" is missing.')
  }
  if (missing(fac_var)) {
    stop('Error! Necessary argument "fac_var" is missing.')
  }

  if (verbose == T) {print("Start data preprocessing.")}
  data <- read.table (file = input, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  # Remove the features (bacteria) whose column sum is 0
  values <- as.data.frame(data[ , variable_col:ncol(data)])
  values <- as.data.frame(apply(values, 2, as.numeric))
  values <- values[ , colSums(values) > 0]
  data <- as.data.frame(cbind(data[ , 1:(variable_col-1)], values))
  non_zero_count <- c()
  for (i in 1:ncol(values)) {
    non_zero_count[i] <- sum(values[ , i] != 0)
  }
  # Here value column is defined by the user
  predictor_names <- (colnames(data))[1: (variable_col - 1)]
  melt_data <- reshape2::melt (data, id = predictor_names)
  # Omit the rows whose value column equals to NA
  melt_data <- melt_data %>% tidyr::drop_na(.data$value)
  # Remove all dots in the bacteria name or it will cause problem
  melt_data$variable <- gsub(".", "_", melt_data$variable, fixed = TRUE)

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
    melt_data[ ,i] <- as.numeric(melt_data[ ,i])
  }
  # Remove the not-used columns in melt_data
  if (!is.null(not_used)) {
    melt_data <- melt_data %>% dplyr::select(-c(not_used))
  }
  # Change the first column name of melt_data to "Individual"
  colnames(melt_data)[1] <- "Individual"
  # Variables are all the bacteria taxanomies
  variables <- unique (melt_data$variable)
  # Extract all the column names from melt_data except for the last two columns
  # which are variables and values, and also exclude "Individual" and test_var
  factors <- colnames(melt_data)[-c(ncol(melt_data), ncol(melt_data)-1)]
  factors <- factors[-which(factors %in% c("Individual", test_var))]
  factor_columns <- match(factors, colnames(melt_data))
  N <- length (variables)
  if (verbose == T) {print("Finish data preprocessing.")}

  ################## Negative binomial model #################
  if (verbose == T) {print("Start running negative binomial models.")}
  Theta <- as.data.frame(matrix(data = NA, nrow = N, ncol = 1))
  for (i in 1:N) { # loop through all variables
    aVariable = variables[i]
    if (verbose == T) {print(i)}
    subdata <- subset(melt_data, variable == aVariable)
    tryCatch({
      fmla2 <- as.formula(paste("value ~ (1| Individual) +", test_var))
      m2 <- glmmTMB::glmmTMB(formula = fmla2, data = subdata, family = nbinom2, REML = F)
      # Extract overdispersion theta out of model
      Theta[i, 1] <- glmmTMB::sigma(m2)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  rownames(Theta) <- gsub(".", "_", variables, fixed = TRUE)
  colnames(Theta) <- c("NB_theta")
  all_info <- cbind(Theta, non_zero_count)
  colnames(all_info)[2] <- "Nonzero_count"
  if (verbose == T) {print("Finish running negative binomial models.")}

  ################## Plot nonzero count v.s. theta ##################
  if (verbose == T) {print("Start plotting.")}
  suppressWarnings(
    plot <- ggplot2::ggplot(all_info, aes(x=Nonzero_count, y = log(NB_theta, base = 2))) +
      geom_point(size = point_size, alpha = 0.9, color = "dodgerblue2") + theme_light() +
      scale_y_discrete(limits = seq(-10, max(log(all_info$NB_theta, base = 2)) + 10, by = y_interval_value)) +
      ggtitle("Non-zero count vs negative binomial theta") +
      theme(plot.title = element_text(size = 18, face = "plain")) +
      xlab("Non-zero count") +
      ylab("log(Theta, base = 2)") +
      scale_x_continuous(breaks = seq(0, nrow(data), by = x_interval_value)) +
      expand_limits(x = c(0, nrow(data)), y = c(-10, max(log(all_info$NB_theta, base = 2)) + 10)))
  return(plot)
  print("Finished successfully!")
}
