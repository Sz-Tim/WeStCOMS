# HABReports Bayesian modelling
# WeStCOMS processing functions
# Tim Szewczyk



#' Load WeStCOMS mesh
#'
#' @param mesh.f Mesh file name
#' @param type Output type: `sf` (default) or `nc`
#'
#' @return Mesh object (ncdf or sf)
#' @export
#'
#' @examples
loadMesh <- function(mesh.f, type="sf") {
  switch(type,
         sf=sf::st_read(mesh.f),
         nc=ncdf4::nc_open(stringr::str_replace(mesh.f, "gpkg$", "nc"))
         )
}









#' Extract hydrodynamic variables from WeStCOMS meshes
#'
#' Given an input dataframe with site locations, dates, hours, and depths, this
#' function extracts the corresponding WeStCOMS data. Output can be specific
#' hours and depths, or summarised for the day or water column (surface to depth
#' *d*) as given in `daySummaryFn` and `depthSummaryFn`.
#'
#' @param sampling.df Dataframe with a row for each sample. Columns must include
#'   `site.id`, `obs.id`, `date` (YYYYMMDD), `hour`, `depth`, `grid`,
#'   `site.elem`, `trinode_1`, `trinode_2`, and `trinode_3`, which can be added
#'   through a spatial join with the mesh files (see vignette)
#' @param westcoms.dir Character vector with directories for WeStCOMS v1 and
#'   WeStCOMS v2 hydrodynamic files
#' @param hydroVars Character vector of hydrodynamic variables to extract
#' @param daySummaryFn Vector of functions to summarise by day, one per var; if
#'   `NULL` (default), single hour is extracted
#' @param depthSummaryFn Vector of functions to summarise across depth, one per
#'   var; if `NULL` (default), single hour is extracted
#' @param regional Logical: does sampling.df include rows to calculate regional
#'   averages? Requires a row for each element to average for each site. For
#'   example, a row for each element within a 1km radius, with the results
#'   averaged
#' @param sep Directory separation character (Windows: '\\', Unix: '/')
#' @param cores Number of cores for extracting in parallel; default is 1
#' @param progress Logical: Show progress bar?
#' @param errorhandling Passed to `foreach`: one of "stop", "remove", or "pass"
#'
#' @return dataframe with site.id, date, and hydrodynamic variables
#' @export
#'
#' @examples
extractHydroVars <- function(sampling.df, westcoms.dir, hydroVars,
                          daySummaryFn=NULL, depthSummaryFn=NULL, regional=F,
                          sep="/", cores=1, progress=T, errorhandling="remove") {
  library(doSNOW); library(foreach)

  dates <- unique(sampling.df$date)

  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max=length(dates), style=3)
  updateProgress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=updateProgress)
  out.df <- foreach(i=seq_along(dates),
                    .packages=c("ncdf4", "tidyverse", "glue", "lubridate"),
                    .export=c("sampling.df", "westcoms.dir", "hydroVars", "daySummaryFn", "depthSummaryFn", "sep"),
                    .options.snow=opts,
                    .errorhandling=errorhandling,
                    .combine=rbind) %dopar% {
    # Load appropriate file
    rows_i <- which(sampling.df$date==dates[i])
    df_i <- sampling.df[rows_i,]
    trinodes_i <- as.matrix(select(df_i, starts_with("trinode")))
    dir_i <- glue("{westcoms.dir[df_i$grid[1]]}{sep}netcdf_{str_sub(dates[i],1,4)}")
    file_i <- dir(dir_i, glue("{dates[i]}.*nc$"))
    nc_i <- nc_open(glue("{dir_i}{sep}{file_i}"))
    n_node <- nc_i$dim$node$len

    # Extract variables, take mean of nodes if necessary
    hydro_extracted <- hydro_out <- vector("list", length(hydroVars)) %>% setNames(hydroVars)
    for(v in seq_along(hydroVars)) {
      hydro_var <- ncvar_get(nc_i, hydroVars[v])
      if(dim(hydro_var)[1] == n_node) {
        hydro_extracted[[hydroVars[v]]] <- meanOfNodes(hydro_var, trinodes_i)
      } else {
        hydro_extracted[[hydroVars[v]]] <- valueOfElement(hydro_var, df_i$site.elem)
      }
    }
    waterDepth <- c(meanOfNodes(ncvar_get(nc_i, "h"), trinodes_i))
    siglay <- abs(ncvar_get(nc_i, "siglay")[1,])
    if("zeta" %in% hydroVars) {
      waterDepth <- hydro_temp[["zeta"]] + waterDepth
    } else {
      waterDepth <- meanOfNodes(ncvar_get(nc_i, "zeta"), trinodes_i) + waterDepth
    }
    nc_close(nc_i)

    # summarise as appropriate
    for(v in seq_along(hydroVars)) {
      var_dims <- dim(hydro_extracted[[hydroVars[v]]])
      if(length(var_dims)==2) {  # [elem,hour]

        if(is.null(daySummaryFn)) {
          # Extract single hour
          if(hydroVars[v]=="short_wave") {
            hydro_out[[hydroVars[v]]] <- map_dbl(1:nrow(df_i),
                                            ~integrateShortWave(hydro_extracted[[hydroVars[v]]], .x, df_i$hour[.x]))
          } else {
            hydro_out[[hydroVars[v]]] <- hydro_extracted[[hydroVars[v]]][,df_i$hour+1]
          }

        } else {
          # Apply daySummaryFn[v]
          if(hydroVars[v]=="short_wave") {
            # Ignore fn, integrate over whole day
            hydro_out[[hydroVars[v]]] <- map_dbl(1:nrow(df_i),
                                            ~integrateShortWave(hydro_extracted[[hydroVars[v]]], .x, 23))
          } else {
            hydro_out[[hydroVars[v]]] <- apply(hydro_extracted[[hydroVars[v]]], 1, daySummaryFn[[v]])
          }
        }
      }
      if(length(var_dims)==3) {  # [elem,siglay,hour]

        if(is.null(daySummaryFn)) {
          # Indexes for extraction
          siglays <- apply(cbind(waterDepth[,c(df_i$hour+1)]) %*% rbind(siglay) - df_i$depth,
                           1, function(x) which.min(abs(x)))
          ii <- cbind(1:nrow(df_i), siglays, df_i$hour+1)

          if(is.null(depthSummaryFn)) {
            # Extract single hour at a single depth
            hydro_out[[hydroVars[v]]] <- hydro_extracted[[hydroVars[v]]][ii]
          } else {
            # Extract single hour from 0-siglay[findDepth], apply depthSummaryFn[v]
            hydro_out[[hydroVars[v]]] <-
              map_dbl(1:nrow(df_i),
                      ~depthSummaryFn[[v]](hydro_extracted[[hydroVars[v]]][.x, 1:ii[.x,2], ii[.x,3]]))
          }

        } else {
          # Indexes for extraction
          siglays <- apply(df_i$depth/waterDepth, 1:2, function(x) which.min(abs(siglay - x)))
          ii <- cbind(rep(1:nrow(df_i), 24), c(siglays), rep(1:24, each=nrow(df_i)))

          if(is.null(depthSummaryFn)) {
            # Apply daySummaryFn[v] to single depth
            hydro_out[[hydroVars[v]]] <- hydro_extracted[[hydroVars[v]]][ii] %>%
              matrix(nrow=nrow(df_i)) %>%
              apply(., 1, daySummaryFn[[v]])
          } else {
            # Apply daySummaryFn[v] to depthSummaryFn[v](values)
            hydro_out[[hydroVars[v]]] <-
              map_dbl(1:nrow(ii),
                    ~depthSummaryFn[[v]](hydro_extracted[[hydroVars[v]]][ii[.x,1], 1:ii[.x,2], ii[.x,3]])) %>%
              matrix(., nrow(df_i)) %>%
              apply(., 1, daySummaryFn[[v]])

          }
        }
      }
    }
    date_out <- df_i %>% select(obs.id) %>%
      bind_cols(hydro_out)
    if(regional) {
      date_out %>% group_by(obs.id) %>%
        summarise(across(where(is.numeric), mean)) %>%
        ungroup
    } else {
      date_out
    }
  }
  close(pb)
  stopCluster(cl)
  return(out.df)
}







#' Find mean of trinodes
#'
#' Find the element mean, not weighting by distance (yet)
#'
#' @param nc.ar Hydrodynamic variable array
#' @param node.mx Trinode index matrix
#'
#' @return Mean value with same dimensions as nc.ar
#' @export
#'
#' @examples
meanOfNodes <- function(nc.ar, node.mx) {
  if (length(dim(nc.ar))==1) {
    node1 <- nc.ar[node.mx[,1]]
    node2 <- nc.ar[node.mx[,2]]
    node3 <- nc.ar[node.mx[,3]]
  } else if(length(dim(nc.ar))==2) {
    node1 <- nc.ar[node.mx[,1],,drop=F]
    node2 <- nc.ar[node.mx[,2],,drop=F]
    node3 <- nc.ar[node.mx[,3],,drop=F]
  } else if(length(dim(nc.ar))==3) {
    node1 <- nc.ar[node.mx[,1],,,drop=F]
    node2 <- nc.ar[node.mx[,2],,,drop=F]
    node3 <- nc.ar[node.mx[,3],,,drop=F]
  }
  return((node1 + node2 + node3)/3)
}





#' Find value of element
#'
#' Find the element value, not weighting by distance (yet)
#'
#' @param nc.ar Hydrodynamic variable array
#' @param elem.id Element id
#'
#' @return Value with same dimensions as nc.ar
#' @export
#'
#' @examples
valueOfElement <- function(nc.ar, elem.id) {
  switch(as.character(length(dim(nc.ar))),
         "1"=nc.ar[elem.id],
         "2"=nc.ar[elem.id,,drop=F],
         "3"=nc.ar[elem.id,,,drop=F])
}





#' Integrate daily accumulation of shortwave radiation
#'
#' @param shortwave.mx Matrix of short_wave variable [elem, hour]
#' @param rowNum Row number of shortwave.mx
#' @param endTime Upper integration time limit (max: 23)
#' @param startTime Lower integration time limit (min: 0S)
#'
#' @return Vector of values, length nrow(shortwave.mx)
#' @export
#'
#' @examples
integrateShortWave <- function(shortwave.mx, rowNum, endTime, startTime=0) {
  # linear interpolation between hours
  # short_wave units = joules/s
  hourlyJoules <- rep(0, ceiling(endTime)-floor(startTime))
  for(hour in seq_along(hourlyJoules)) {
    varStart <- shortwave.mx[rowNum,hour] # indexes = 1:24, hours = 0:23
    varEnd <- shortwave.mx[rowNum,hour+1]
    dt <- 3600
    if(hour == ceiling(endTime)) {
      propHour <- endTime %% 1
      varEnd <- dVar*propHour + varStart
      dt <- 3600*propHour
    }
    dVar <- varEnd - varStart
    hourlyJoules[hour] <- dt*min(varStart, varEnd) + (0.5*dt*abs(dVar))
  }
  return(sum(hourlyJoules))
}






#' Find 90th quantile
#'
#' Same as quantile(x, probs=0.9), but works with just a single argument
#'
#' @param x Vector
#'
#' @return
#' @export
#'
#' @examples
q90 <- function(x) {
  quantile(x, probs=0.9)
}
