
#' Create Mesh Files from NetCDF File
#'
#' This function creates mesh files from a NetCDF file. Note that there is a two
#' stage process where only a *_TEMP.gpkg file is output first, allowing the
#' manual identification of open boundary elements which are then saved to a
#' .csv and provided to the `open_elems_csv` argument
#'
#' @param nc_file Character. The path to the NetCDF file.
#' @param open_elems_csv Character. The path to the CSV file containing open boundary elements. Default is NULL.
#' @param out_dir Character. The output directory where the mesh files will be saved.
#' @param domain Character. The domain name for the mesh.
#'
#' @return Character. A message indicating the status of the mesh creation.
#' @export
#'
#' @examples
#' create_mesh_from_nc(nc_file="path/to/nc_file.nc", open_elems_csv="path/to/open_elems.csv", out_dir="~/output", domain="WeStCOMS2")
create_mesh_files_from_nc <- function(nc_file, open_elems_csv=NULL, out_dir, domain) {

  library(tidyverse); library(ncdf4); library(sf); library(glue)

  day2 <- nc_open(nc_file)
  if(!is.null(open_elems_csv)) {
    open_elems <- read_csv(open_elems_csv)
  }


  # extract dims ------------------------------------------------------------

  # Variables to extract, following Tom's naming convention for clarity (see
  # WeStCOMS2_mesh.nc) and to avoid rewriting a bunch of Java code
  nc <- list(
    uvnode=cbind(ncvar_get(day2, "lonc"), ncvar_get(day2, "latc")),
    nodexy=cbind(ncvar_get(day2, "lon"), ncvar_get(day2, "lat")),
    trinodes=ncvar_get(day2, "nv"),
    nbe=ncvar_get(day2, "nbe"),
    depthNodexy=ncvar_get(day2, "h"),
    depthUvnode=ncvar_get(day2, "h_center"),
    siglay=ncvar_get(day2, "siglay")[1,],
    siglev=ncvar_get(day2, "siglev")[1,]
  )
  nc_close(day2)
  nc$boundaryNodesAll <- which(rowSums(nc$nbe==0) > 0)
  if(!is.null(open_elems_csv)) {
    nc$boundaryNodesOpen <- sort(unique(unlist(open_elems %>% select(starts_with("trinode")))))
  }
  nc$uvnode_os <- as_tibble(nc$uvnode) %>%
    st_as_sf(coords=c("V1", "V2"), crs=4326) %>%
    st_transform(27700) %>% st_coordinates
  nc$nodexy_os <- as_tibble(nc$nodexy) %>%
    st_as_sf(coords=c("V1", "V2"), crs=4326) %>%
    st_transform(27700) %>% st_coordinates


  # assign dims -------------------------------------------------------------

  elems <- ncdim_def("elems", "element", 1:nrow(nc$uvnode))
  two <- ncdim_def("two", "two", 1:2)
  node <- ncdim_def("node", "nodes", 1:nrow(nc$nodexy))
  nele <- elems
  three <- ncdim_def("three", "three", 1:3)
  siglayers <- ncdim_def("siglayers", "sigma layers", 1:length(nc$siglay))
  siglevels <- ncdim_def("siglevels", "sigma levels", 1:length(nc$siglev))
  bnode <- ncdim_def("bnode", "boundary", 1:length(nc$boundaryNodesAll))
  if(!is.null(open_elems_csv)) {
    obcnode <- ncdim_def("obcnode", "boundary", 1:length(nc$boundaryNodesOpen))
  }



  # create vars -------------------------------------------------------------

  mesh.vars <- list(
    uvnode=ncvar_def("uvnode", "degrees", list(elems, two)),
    nodexy=ncvar_def("nodexy", "degrees", list(node, two)),
    uvnode_os=ncvar_def("uvnode_os", "m", list(elems, two)),
    nodexy_os=ncvar_def("nodexy_os", "m", list(node, two)),
    trinodes=ncvar_def("trinodes", "", list(nele, three), prec="integer"),
    nbe=ncvar_def("nbe", "", list(nele, three), prec="integer"),
    depthUvnode=ncvar_def("depthUvnode", "m", elems),
    depthNodexy=ncvar_def("depthNodexy", "m", node),
    boundaryNodesAll=ncvar_def("boundaryNodesAll", "", bnode, prec="integer"),
    siglay=ncvar_def("siglay", "proportion", siglayers),
    siglev=ncvar_def("siglev", "proportion", siglevels)
  )
  if(!is.null(open_elems_csv)) {
    mesh.vars$boundaryNodesOpen=ncvar_def("boundaryNodesOpen", "", obcnode, prec="integer")
  }



  # create mesh.nc ----------------------------------------------------------

  if(!is.null(open_elems_csv)) {
    mesh.nc <- nc_create(glue("{out_dir}/{domain}_mesh.nc"), mesh.vars)
    iwalk(names(mesh.vars), ~ncvar_put(mesh.nc, .x, nc[[.x]]))
    nc_close(mesh.nc)
  }



  # create mesh.gpkg --------------------------------------------------------

  suffix <- ifelse(is.null(open_elems_csv), "_TEMP", "")
  mesh.sf <- tibble(i=1:nrow(nc$uvnode_os),
                    depth=nc$depthUvnode,
                    lonc=nc$uvnode_os[,1],
                    latc=nc$uvnode_os[,2],
                    trinode_1=nc$trinodes[,1],
                    trinode_2=nc$trinodes[,2],
                    trinode_3=nc$trinodes[,3]) %>%
    rowwise() %>%
    mutate(coords=list(nc$nodexy_os[c(trinode_1, trinode_2, trinode_3, trinode_1),]),
           geom=list(st_polygon(list(coords)))) %>%
    ungroup %>%
    st_as_sf(crs=27700) %>%
    mutate(area=st_area(.)) %>%
    select(i, area, depth, trinode_1, trinode_2, trinode_3, geom)
  write_sf(mesh.sf, glue("{out_dir}/{domain}_mesh{suffix}.gpkg"))

  if(!is.null(open_elems_csv)) {
    mesh.footprint <- mesh.sf %>% st_union
    write_sf(mesh.footprint, glue("{out_dir}/{domain}_mesh_footprint.gpkg"))
  }

  if(is.null(open_elems_csv)) {
    return(glue("Please open {out_dir}/{domain}_mesh{suffix}.gpkg and identify open boundary elements."))
  } else {
    return(glue("Created {domain}_mesh.gpkg, {domain}_mesh_footprint.gpkg, and {domain}_mesh.nc in {out_dir}"))
  }

}
