#' Unlist confound and inverse confound tables, turn them into tables
#' @param x The list to be unlisted and turned into table

unlist_table <- function(x, N, variables) {
  unlist = unlist(x, use.names=FALSE) # output only the numbers
  mode(unlist) <- "numeric" # it's numeric now!
  unlist_noNA <-unlist[!is.na(unlist)] # Exclude NA in unlist
  unlist_withname = unlist(x) # numbers with names
  unlist_name <- names(unlist_withname) # outputs the name of every element
  # if there's NA, then remove names of NA value
  if (TRUE %in% is.na(unlist)) {
    unlist_name_noNA <- unlist_name[-which(is.na(unlist))]
  } else {
    unlist_name_noNA <- unlist_name
  }

  # Create a empty matrix first
  ps_lm_m <- data.frame(matrix(nrow = length(unlist), ncol = 2))
  # Specify column names and fill in the values
  colnames(ps_lm_m) <- c("name", "p_lm")
  # Fill in the values
  ps_lm_m[ , 1] <- unique(unlist_name)
  ps_lm_m[ , 2] <- unlist

  # Split the names of each p_lm, ex: Actinobacteria.PPI becomes "Actinobacteria" "PPI"
  library(stringr)
  out <- stringr::str_split_fixed(ps_lm_m$name, pattern = stringr::fixed("."), n = 2)
  # Add the 2 new columns with split names to the matrix
  ps_lm_m2 <- cbind(out, ps_lm_m)
  # Remove the original unsplit name column, and then name the columns
  ps_lm_m2 $name = NULL
  colnames(ps_lm_m2) <- c("feature", "factor", "p_lm")
  # Reshaping the matrix so that can do fdr correction in the next step
  cast.ps_lm_m2 <- reshape2::dcast(data = ps_lm_m2, feature ~ factor)
  if (sum(is.na(cast.ps_lm_m2$Var.2)) == N) {
    cast.ps_lm_m2 <- as.data.frame(cast.ps_lm_m2[ , -2])
  }
  rownames(cast.ps_lm_m2) <- cast.ps_lm_m2$feature
  cast.ps_lm_m2$feature <- NULL
  Ps_conf_model_unlist <- cast.ps_lm_m2
  #Reorder the rows to make it follow the order of variables
  Ps_conf_model_unlist <- as.data.frame(Ps_conf_model_unlist[match(variables,rownames(Ps_conf_model_unlist)), ])
  return(Ps_conf_model_unlist)
  }
