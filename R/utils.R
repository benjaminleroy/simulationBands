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


#' check if a file exists with a certain name in a certain directory
#'
#' Uses \code{stringr::str_detect}, if you want exact, use \code{^} and \code{$}
#' before and after your string
#'
#' @param file_name file name or part of string name (if not unique will just
#' check that at least a single file matches the provided name)
#' @param dir file directory (default is current directory)
#' @param .logic boolean - if should return a boolean (\code{TRUE}) or should
#' return a vector a files that match \code{file_name} string (\code{FALSE}).
#'
#' @return see \code{.logic} parameter.
#' @export
check_file_exists <- function(file_name,dir = ".", .logic = T){
  files <- list.files(path = dir)
  detection <- stringr::str_detect(files, pattern = file_name)
  if (!.logic){
    return(files[detection])
  } else {
    return(any(detection))
  }

}
