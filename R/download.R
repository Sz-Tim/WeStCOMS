

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
#' @param OPeNDAP Character string to append to url to only download select
#'   variables. Default is NULL, which downloads full file
#'
#' @return None. The function downloads files to the specified output directory.
#' @export
#'
#' @examples
#' download_from_thredds(domain="westcoms2", year=2024, out_dir="~/data", cores=4, overwrite=FALSE, start_ymd=NULL, end_ymd=NULL)
download_from_thredds <- function(domain="westcoms3",
                                  year=2025,
                                  out_dir="~/",
                                  cores=4,
                                  overwrite=FALSE,
                                  start_ymd=NULL,
                                  end_ymd=NULL,
                                  OPeNDAP=NULL) {
  library(tidyverse); library(ncdf4); library(lubridate); library(glue);
  library(xml2); library(rvest)
  library(doFuture); library(progressr)

  if(!is.null(start_ymd)) year <- year(start_ymd)
  domain_name <- switch(domain,
                        "norscoms"="NORSCOMS1",
                        "westcoms3"="WeStCOMS3",
                        "westcoms2-nc"="WeStCOMS2",
                        "etive28"="etive28",
                        "swan"="WeStCOMS2")
  if(is.null(OPeNDAP)) {
    if(domain=="swan") {
      catalog <- glue("https://thredds.sams.ac.uk/thredds/catalog/scoats-westcoms2/SWAN/Archive_forecast/netcdf_{year}F/catalog.html")
      file_base <- glue("https://thredds.sams.ac.uk/thredds/fileServer/scoats-westcoms2/SWAN/Archive_forecast/netcdf_{year}F/")
      out_dir <- glue("{out_dir}/WeStCOMS2/SWAN/Archive_forecast/netcdf_{year}F/")
    } else {
      catalog <- glue("https://thredds.sams.ac.uk/thredds/catalog/scoats-{domain}/Archive/netcdf_{year}/catalog.html")
      file_base <- glue("https://thredds.sams.ac.uk/thredds/fileServer/scoats-{domain}/Archive/netcdf_{year}/")
      out_dir <- glue("{out_dir}/{domain_name}/Archive/netcdf_{year}/")
    }
  } else {
    if(domain=="swan") {
      catalog <- glue("https://thredds.sams.ac.uk/thredds/catalog/scoats-westcoms2/SWAN/Archive_forecast/netcdf_{year}F/catalog.html")
      file_base <- glue("https://thredds.sams.ac.uk/thredds/fileServer/scoats-westcoms2/SWAN/Archive_forecast/netcdf_{year}F/")
      out_dir <- glue("{out_dir}/WeStCOMS2/SWAN/Archive_forecast/netcdf_{year}F/")
    } else {
      catalog <- glue("https://thredds.sams.ac.uk/thredds/catalog/scoats-{domain}/Archive/netcdf_{year}/catalog.html")
      file_base <- glue("https://thredds.sams.ac.uk/thredds/dodsC/scoats-{domain}/Archive/netcdf_{year}/")
      out_dir <- glue("{out_dir}/{domain_name}/Archive/netcdf_{year}/")
    }
  }

  dir.create(out_dir, showWarnings=F, recursive=T)

  if(is.null(start_ymd)) {
    start_ymd <- ymd(glue("{year}-01-01"))
  }
  if(is.null(end_ymd)) {
    end_ymd <- ymd(glue("{year}-12-31"))
  }

  dates <- seq(start_ymd, end_ymd, by=1) |> format("%Y%m%d")
  if(!overwrite) {
    f <- str_split_fixed(dir(out_dir, "*.nc", recursive=T), "_", 4)[,2]
    dates <- dates[! dates %in% f]
  }

  westcom_links <- read_html(catalog) |> html_nodes("a") |> html_attr("href")
  nc_links <- map(dates, ~grep(.x, westcom_links, value=T) |> basename()) |>
    do.call('c', args=_)

  plan(multisession, workers=cores)
  handlers(global=T)
  download_nc(nc_links, file_base, out_dir, OPeNDAP)
  plan(sequential)
}





#' Download NetCDF Files
#'
#' This function downloads NetCDF files from a specified list of links.
#'
#' @param nc_links Character vector. A list of NetCDF file links to download.
#' @param file_base Character. The base URL for the files to be downloaded.
#' @param out_dir Character. The output directory where the downloaded files will be saved.
#' @param OPeNDAP Character string to append to url to only download select
#'   variables. Default is NULL, which downloads full file
#'
#' @return None. The function downloads files to the specified output directory.
#' @export
#'
#' @examples
#' download_nc(nc_links=c("file1.nc", "file2.nc"), file_base="https://example.com/files/", out_dir="~/data/")
download_nc <- function(nc_links, file_base, out_dir, OPeNDAP=NULL) {
  library(tidyverse); library(ncdf4); library(lubridate); library(glue);
  library(xml2); library(rvest)
  library(doFuture); library(progressr)

  if(is.null(OPeNDAP)) {
    p <- progressor(along=nc_links)
    foreach(i=nc_links,
            .options.future=list(globals=structure(TRUE, add=c("file_base", "out_dir")))
    ) %dofuture% {
      options(timeout=3600*5)
      if( !file.exists(glue("{out_dir}{i}")) ) {
        download.file(glue("{file_base}{i}"), glue("{out_dir}{i}"), method="auto")
      }
      p(sprintf("i=%s", i))
    }
  } else {
    p <- progressor(along=nc_links)
    foreach(i=nc_links,
            .errorhandling="pass",
            .options.future=list(globals=structure(TRUE, add=c("file_base", "out_dir", "OPeNDAP")))
    ) %dofuture% {
      options(timeout=3600*5)
      if( !file.exists(glue("{out_dir}{i}")) ) {
        download_nc_opendap(glue("{file_base}{i}{OPeNDAP}"), glue("{out_dir}{i}"))
      }
      p(sprintf("i=%s", i))
    }
  }

}

  download_nc_opendap <- function(f_src, f_dest) {
    library(ncdf4)

    nc_in <- nc_open(f_src)
    ## ---- Copy dimensions ----
    dims <- lapply(nc_in$dim, function(d) {
      ncdim_def(
        name = d$name,
        units = d$units,
        vals  = d$vals,
        create_dimvar = TRUE
      )
    })
    names(dims) <- names(nc_in$dim)

    ## ---- Copy variables ----
    vars <- lapply(nc_in$var, function(v) {
      # Get dimension objects in correct order
      dim_list <- dims[v$dimids + 1]
      ncvar_def(
        name  = v$name,
        units = v$units,
        dim   = dim_list,
        missval = v$missval,
        longname = v$longname,
        prec = v$prec
      )
    })
    names(vars) <- names(nc_in$var)

    # Create new file with copied metadata
    nc_out <- nc_create(f_dest, vars)

    ## ---- Copy global attributes ----
    for (attname in names(nc_in$gatts)) {
      ncatt_put(nc_out, 0, attname, nc_in$gatts[[attname]])
    }

    ## ---- Copy variable attributes + data ----
    for (vname in names(nc_in$var)) {
      # Copy attributes
      atts <- ncatt_get(nc_in, vname)
      for (att in names(atts)) {
        ncatt_put(nc_out, vname, att, atts[[att]])
      }

      # Copy actual data
      data <- ncvar_get(nc_in, vname)
      ncvar_put(nc_out, vname, data)
    }

    nc_close(nc_out)
    nc_close(nc_in)
  }


