# WeStCOMS
# Miscellaneous helper functions
# Tim Szewczyk


#' Lag multiple variables at once
#'
#' From https://stackoverflow.com/questions/55814028/multiple-lags-with-dplyr
#'
#' @param data Dataframe
#' @param ... Unquoted variable names to lag
#' @param n Number of lags
#'
#' @return dataframe with lags added
#' @export
#'
multijetlag <- function(data, ..., n=10){
  library(rlang)
  variable <- enquos(...)

  indices <- seq_len(n)
  combos <- crossing(indices, var =as.list(variable))

  quosures <- map2(combos$indices, combos$var,
                   ~quo(lag(!!.y, !!.x)) ) %>%
    set_names(paste(map_chr(combos$var, quo_text), combos$indices, sep = "_"))
  mutate( data, !!!quosures )

}
