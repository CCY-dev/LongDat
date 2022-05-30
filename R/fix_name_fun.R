#' Replace the symbols in variable and covariate names in raw input
#' @param z A character vector. This is the character vector
#'          that needs to be changed.
#' @importFrom rlang .data
#' @name fix_name_fun

fix_name_fun <- function(z) {
  n1 <- gsub(x = z, pattern = "-", replacement = "_", fixed = TRUE)
  n2 <- gsub(x = n1, pattern = "#", replacement = "_", fixed = TRUE)
  n3 <- gsub(x = n2, pattern = "%", replacement = "percent", fixed = TRUE)
  n4 <- gsub(x = n3, pattern = "*", replacement = "_", fixed = TRUE)
  n5 <- gsub(x = n4, pattern = "(", replacement = "_", fixed = TRUE)
  n6 <- gsub(x = n5, pattern = ")", replacement = "", fixed = TRUE)
  n7 <- gsub(x = n6, pattern = "~", replacement = "to", fixed = TRUE)
  n8 <- gsub(x = n7, pattern = ".", replacement = "_", fixed = TRUE)
  n9 <- gsub(x = n8, pattern = "&", replacement = "and", fixed = TRUE)
  n10 <- gsub(x = n9, pattern = "$", replacement = "_", fixed = TRUE)
  n11 <- gsub(x = n10, pattern = "@", replacement = "at", fixed = TRUE)
  n12 <- gsub(x = n11, pattern = "!", replacement = "_", fixed = TRUE)
  n13 <- gsub(x = n12, pattern = "+", replacement = "plus", fixed = TRUE)
  n14 <- gsub(x = n13, pattern = "/", replacement = "_", fixed = TRUE)
  n15 <- gsub(x = n14, pattern = "?", replacement = "QuestionMark",
              fixed = TRUE)
  n16 <- gsub(x = n15, pattern = ",", replacement = "_", fixed = TRUE)
  n17 <- gsub(x = n16, pattern = ">", replacement = "GreaterThan", fixed = TRUE)
  n18 <- gsub(x = n17, pattern = "<", replacement = "LessThan", fixed = TRUE)
  n19 <- gsub(x = n18, pattern = "]", replacement = "_", fixed = TRUE)
  n20 <- gsub(x = n19, pattern = "[", replacement = "_", fixed = TRUE)
  n21 <- gsub(x = n20, pattern = " ", replacement = "_",fixed = TRUE)
  n22 <- gsub(x = n21, pattern = ":", replacement = "_",fixed = TRUE)
  n23 <- gsub(x = n22, pattern = ";", replacement = "_",fixed = TRUE)
  n24 <- gsub(x = n23, pattern = "=", replacement = "_",fixed = TRUE)
  n25 <- gsub(x = n24, pattern = "^", replacement = "_",fixed = TRUE)
}


