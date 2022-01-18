#' Create input master table from metadata and feature tables for
#' longdat_disc() and longdat_cont()
#' @param metadata_table A data frame whose columns consist of
#' sample identifiers (sample_ID), individual, time point and other meta data.
#' Each row corresponds to one sample_ID. Metadata table should have the same
#' number of rows as feature table does. Please avoid using characters that
#' don't belong to ASCII printable characters for the column names.
#' @param feature_table A data frame whose columns only consist of
#'  sample identifiers (sample_ID) and features
#'  (dependent variables, e.g. microbiome). Each row corresponds to
#'  one sample_ID. Please do not include any columns other than
#'  sample_ID and features. Please avoid using characters that
#' don't belong to ASCII printable characters for the column names.
#' Also, feature table should have the same number
#'  of rows as metadata table does.
#' @param sample_ID The name of the column which stores sample identifiers.
#' Please make sure that sample_IDs are unique for each sample, and that
#' metadata and feature tables have the same sample_IDs. If sample_IDs don't
#' match between the two tables, it will fail to join them together. This
#' should be a string, e.g. "Sample_ID"
#' @param individual The name of the column which stores individual information
#' in the metadata table. This should be a string, e.g. "Individual"
#' @param keep_id A boolean vector indicating whether keep sample_ID
#' column in the output master table. The default is FALSE.
#' @export
#' @import tidyr
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom magrittr '%>%'
#' @name make_master_table
#' @return a data frame which complies with the required format of an input
#' data frame for longdat_disc() and longdat_cont().
#' @details
#'  This function joins metadata and feature tables by the sample_ID
#'  column. Users can create master tables compatible with the format of
#'  longdat_disc() and longdat_cont() input easily.This function outputs a
#'  master table with individual as the first column, followed by time point
#'  and other metadata, and then by feature columns.
#' @examples
#' test_master <- make_master_table(
#' metadata_table = LongDat_disc_metadata_table,
#' feature_table = LongDat_disc_feature_table,
#' sample_ID = "Sample_ID",
#' individual = "Individual")

#utils::globalVariables(c(""))

make_master_table <- function(metadata_table, feature_table,
                              sample_ID, individual,
                              keep_id = FALSE) {

  if (missing(metadata_table)) {
    stop('Error! Necessary argument "metadata_table" is missing.')}
  if (missing(feature_table)) {
    stop('Error! Necessary argument "feature_table" is missing.')}
  if (missing(sample_ID)) {
    stop('Error! Necessary argument "sample_ID" is missing.')}
  if (missing(individual)) {
    stop('Error! Necessary argument "individual" is missing.')}

  # Check if nrow of the two tables are the same
  if (nrow(metadata_table) != nrow(feature_table)) {
    stop("Error! Row numbers aren't consistent between the
         metadata_table and feature_table.")}

  # Check if the sample_ID is unique for every sample
  if (length(unique(metadata_table[ , sample_ID])) != nrow(metadata_table)) {
    stop("Error! There are repeated sample_IDs in the metadata_table.
    Please make sure that sample_IDs are unique for each sample.")}
  if (length(unique(feature_table[ , sample_ID])) != nrow(feature_table)) {
    stop("Error! There are repeated sample_IDs in the feature_table
    Please make sure that sample_IDs are unique for each sample.")}

  # Check if the sample_IDs match between the two tables
  diff <- setdiff(metadata_table[ , sample_ID], feature_table[ , sample_ID])
  if (length(diff) > 0) {
    stop("Error! The sample_IDs don't match between the
         metadata_table and feature_table. ")}

  # Arrange the rows by sample_ID and make sample_ID the first column
  metadata_table <- metadata_table %>%
    dplyr::arrange(sample_ID) %>%
    dplyr::select(all_of(sample_ID), dplyr::everything())
  feature_table <- feature_table %>%
    dplyr::arrange(sample_ID) %>%
    dplyr::select(all_of(sample_ID), dplyr::everything())

  # Join the two tables
  master <- inner_join(metadata_table, feature_table, by = sample_ID)

  # Discard the sample_ID column (or keep it by user's choice)
  if (keep_id == TRUE) {
    master <- master
  } else {
    master <- master %>%
      dplyr::select(-dplyr::all_of(sample_ID))
  }

  # Make Individual as the first column
  master <- master %>%
    dplyr::select(dplyr::all_of(individual), dplyr::everything())

  print("Finished creating master table successfully!")
  return(master)
}
