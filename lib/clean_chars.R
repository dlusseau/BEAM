
#' Cleans up excessive whitespaces and converts to lowercase
#' 
#' @description
#' `clean_chars` does three things to `x`:
#' 1. converts values in all character columns to lower case
#' 2. removes all leading and trailing whitespaces from values in char columns
#' 3. replaces multiple sequential whitespaces in values in char columns with a single space
#'
#' @details
#' The object `x` is modified in place (by reference), but is also returned.
#' 
#' @param x A data.table
#' @returns The modified data.table
#' @export
# 
clean_chars <- function(x) {
    char_cols <- which(sapply(colnames(x), function(col) is.character(x[[col]])))
    
    trim <- function(x) {
        x <- gsub("^ +| +$", "", x) # leading or trailing whitespaces
        x <- gsub(" +", " ", x) # multiple sequential whitespaces
        x
    }
    
    # convert all character columns to lower case
    x[, (char_cols) := lapply(.SD, tolower), .SDcols = char_cols] 
    # remove all leading and trailing whitespaces from values in char columns
    x[, (char_cols) := lapply(.SD, trim), .SDcols = char_cols]
    x
}
