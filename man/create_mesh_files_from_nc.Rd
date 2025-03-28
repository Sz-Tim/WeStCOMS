% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mesh.R
\name{create_mesh_files_from_nc}
\alias{create_mesh_files_from_nc}
\title{Create Mesh Files from NetCDF File}
\usage{
create_mesh_files_from_nc(nc_file, open_elems_csv, out_dir, domain)
}
\arguments{
\item{nc_file}{Character. The path to the NetCDF file.}

\item{open_elems_csv}{Character. The path to the CSV file containing open boundary elements. Default is NULL.}

\item{out_dir}{Character. The output directory where the mesh files will be saved.}

\item{domain}{Character. The domain name for the mesh.}
}
\value{
Character. A message indicating the status of the mesh creation.
}
\description{
This function creates mesh files from a NetCDF file. Note that there is a two
stage process where only a *_TEMP.gpkg file is output first, allowing the
manual identification of open boundary elements which are then saved to a
.csv and provided to the `open_elems_csv` argument
}
\examples{
create_mesh_from_nc(nc_file="path/to/nc_file.nc", open_elems_csv="path/to/open_elems.csv", out_dir="~/output", domain="WeStCOMS2")
}
