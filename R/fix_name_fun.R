#' Replace the symbols in variable and confounder names in raw input
#' @param z A character vector. This is the character vector that needs to be changed.

fix_name_fun <- function(z) {
  n1 <- gsub(x = z, pattern = "-", replacement = "_", fixed = T)
  n2 <- gsub(x = n1, pattern = "#", replacement = "_", fixed = T)
  n3 <- gsub(x = n2, pattern = "%", replacement = "_", fixed = T)
  n4 <- gsub(x = n3, pattern = " ", replacement = "_",fixed = T)
  n5 <- gsub(x = n4, pattern = "(", replacement = "_", fixed = T)
  n6 <- gsub(x = n5, pattern = ")", replacement = "", fixed = T)
  n7 <- gsub(x = n6, pattern = "%", replacement = "percent", fixed = T)
  n8 <- gsub(x = n7, pattern = ".", replacement = "_", fixed = T)
  n9 <- gsub(x = n8, pattern = "&", replacement = "and", fixed = T)
  n10 <- gsub(x = n9, pattern = "$", replacement = "_", fixed = T)
  n11 <- gsub(x = n10, pattern = "@", replacement = "at", fixed = T)
  n12 <- gsub(x = n11, pattern = "!", replacement = "_", fixed = T)
  n13 <- gsub(x = n12, pattern = "+", replacement = "add", fixed = T)
  n14 <- gsub(x = n13, pattern = "/", replacement = "_", fixed = T)
  n15 <- gsub(x = n14, pattern = "?", replacement = "_", fixed = T)
  n16 <- gsub(x = n15, pattern = ",", replacement = "_", fixed = T)
  n17 <- gsub(x = n16, pattern = ">", replacement = "_", fixed = T)
  n18 <- gsub(x = n17, pattern = "<", replacement = "_", fixed = T)
  n19 <- gsub(x = n18, pattern = "*", replacement = "_", fixed = T)
  n20 <- gsub(x = n19, pattern = ".", replacement = "_", fixed = T)
}
