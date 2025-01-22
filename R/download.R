

#' Download NetCDF Files from THREDDS Server
#'
#' This function downloads NetCDF files from a specified THREDDS server catalog
#' for a given domain and year.
#'
#' @param domain Character. The domain to download data from. Options are
#'   "westcoms2" or "etive28". Default is "westcoms2".
#' @param year Numeric. The year of the data to download. Default is 2024.
#' @param out_dir Character. The output directory where the downloaded files
#'   will be saved. Default is "~/".
#' @param cores Numeric. The number of cores to use for parallel processing.
#'   Default is 4.
#' @param overwrite Logical. Whether to overwrite existing files. Default is
#'   FALSE.
#' @param start_ymd Date. The start date for the data to download. Default is
#'   NULL, which sets it to January 1st of the specified year.
#' @param end_ymd Date. The end date for the data to download. Default is NULL,
#'   which sets it to December 31st of the specified year.
#'
#' @return None. The function downloads files to the specified output directory.
#' @export
#'
#' @examples
#' download_from_thredds(domain="westcoms2", year=2024, out_dir="~/data", cores=4, overwrite=FALSE, start_ymd=NULL, end_ymd=NULL)
download_from_thredds <- function(domain="westcoms2",
                                  year=2024,
                                  out_dir="~/",
                                  cores=4,
                                  overwrite=FALSE,
                                  start_ymd=NULL,
                                  end_ymd=NULL) {
  library(tidyverse); library(ncdf4); library(lubridate); library(glue);
  library(xml2); library(rvest)
  library(doFuture); library(progressr)

  if(!is.null(start_ymd)) year <- year(start_ymd)
  domain_name <- switch(domain,
                        westcoms2="WeStCOMS2",
                        etive28="etive28")
  catalog <- glue("https://thredds.sams.ac.uk/thredds/catalog/scoats-{domain}/Archive/netcdf_{year}/catalog.html")
  file_base <- glue("https://thredds.sams.ac.uk/thredds/fileServer/scoats-{domain}/Archive/netcdf_{year}/")
  out_dir <- glue("{out_dir}/{domain_name}/Archive/netcdf_{year}/")
  dir.create(out_dir, showWarnings=F, recursive=T)

  if(is.null(start_ymd)) {
    start_ymd <- ymd(glue("{year}-01-01"))
  }
  if(is.null(end_ymd)) {
    end_ymd <- ymd(glue("{year}-12-31"))
  }

  dates <- seq(start_ymd, end_ymd, by=1) |> format("%Y%m%d")
  if(!overwrite) {
    f <- str_split_fixed(dir(out_dir, "*.nc", recursive=T), "_", 4)[,3]
    dates <- dates[! dates %in% f]
  }

  westcom_links <- read_html(catalog) |> html_nodes("a") |> html_attr("href")
  nc_links <- map(dates, ~grep(.x, westcom_links, value=T) |> basename()) |>
    do.call('c', args=_)

  plan(multisession, workers=cores)
  handlers(global=T)
  download_nc(nc_links, file_base, out_dir)
  plan(sequential)
}





#' Download NetCDF Files
#'
#' This function downloads NetCDF files from a specified list of links.
#'
#' @param nc_links Character vector. A list of NetCDF file links to download.
#' @param file_base Character. The base URL for the files to be downloaded.
#' @param out_dir Character. The output directory where the downloaded files will be saved.
#'
#' @return None. The function downloads files to the specified output directory.
#' @export
#'
#' @examples
#' download_nc(nc_links=c("file1.nc", "file2.nc"), file_base="https://example.com/files/", out_dir="~/data/")
download_nc <- function(nc_links, file_base, out_dir) {
  library(tidyverse); library(ncdf4); library(lubridate); library(glue);
  library(xml2); library(rvest)
  library(doFuture); library(progressr)
  p <- progressor(along=nc_links)
  foreach(i=nc_links,
          .options.future=list(globals=structure(TRUE, add=c("file_base", "out_dir")))
  ) %dofuture% {
    options(timeout=3600*5)
    if( !file.exists(glue("{out_dir}{i}")) ) {
      download.file(glue("{file_base}{i}"), glue("{out_dir}{i}"), method="wget")
    }
    p(sprintf("i=%s", i))
  }
}
