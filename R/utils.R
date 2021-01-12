#' converts interval values from the \code{cut} function to scalars. Works
#' inside \code{tidyverse}'s \code{mutate} function.
#'
#' @param id value or set of values
#' @param .lower boolean to isolate the lower bound or upper bound. Default is
#' \code{TRUE} - aka .lower bound.
#'
#' @return vector of scalar values for the bound
#' @export
cut_to_numeric <- function(id, .lower = TRUE){
  if (.lower){
    initial_extract <- id %>%
      as.character() %>%
      stringr::str_extract(paste0("(\\(-*[0-9]+\\.*[0-9]*(e-*[0-9]{2}){0,1}\\,)|",
                                  "(\\(-Inf\\,)"))
  } else {
    initial_extract <- id %>%
      as.character() %>%
      stringr::str_extract(paste0("(\\,-*[0-9]+\\.*[0-9]*(e-*[0-9]{2}){0,1}\\])|",
                                  "(\\,Inf])"))
  }
  out_numeric <- initial_extract %>%
    stringr::str_extract("(-*[0-9]+\\.*[0-9]*(e-*[0-9]{2}){0,1})|(-Inf)") %>%
    as.numeric()

  return(out_numeric)
}
