#' Determine if the tests are significant based on the previous p value table
#' @param file The p value table file you get from pvalue_calc function
#' @param pair Pair is the two time points to compare at. For example, 13 means
#'             timepoint 1 and 3. Here list all the pairs to compare with c(...)
#' @return None
#' @examples sig_calc(file, pair)
#' @export

sig_calc <- function(file, pair){
  N <- length (unique(pdata$variable))
  # Write a empty table for next steps to fill in
  sig_table <- data.frame ("Feature" = pdata$variable,
                           "Space" = tags[as.numeric(grep(as.character(file), tags$V1, ignore.case = F)), 2],
                           "Effect_12" = rep (NA, N), "Effect_23" = rep (NA, N), "Effect_13" = rep (NA, N),
                           "FDR (intervention model comparison)" = pdata$qModel, "Fasting effect size (Cliff's delta)" = pdata$cd12,
                           "Fasting FDR (MWU post-hoc, N=3)" = pdata$pPost12, "Recovery effect size (Cliff's delta)" = pdata$cd23,
                           "Recovery FDR (MWU post-hoc, N=3)" = pdata$pPost23, "Study effect size (Cliff's delta)" = pdata$cd13,
                           "Study FDR (MWU post-hoc, N=3)" = pdata$pPost13, "Signal" = rep (NA, N))

  # The criteria is that the qmodel of the feature < 0.1, and that the post-hoc p value < 0.05
  for (i in pair) {
    pPost_pair <- as.character(paste0("pPost", pair))
    cd_pair <- as.character(paste0("cd", pair))
    effect_pair <- as.character(paste("Effect", pair, sep = "_"))

    sig_table[ , effect_pair] <- (
      ifelse (test = (!(is.na(pdata$qModel)) & pdata$qModel < 0.1 &
                        !(is.na(pdata[ , pPost_pair])) & pdata[ , pPost_pair] < 0.05 ),
              yes = ifelse(test = pdata[ , cd_pair] < 0,
                           yes = "Yes-Depleted",
                           no = "Yes-Enriched"),
              no = "NS"))
  }
  sig_table$Signal <-
    ifelse(test = (sig_table$Effect_12 != "NS" | sig_table$Effect_23 != "NS" | sig_table$Effect_13 != "NS"),
           yes = "OK", no = "Insignificant")
  write.table (file = "sig_table", sig_table, sep = "\t", quote = F, row.names = F)
  print("Great! Your significance result is now in your current directory.")
}
