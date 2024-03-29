% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_raw_data.R
\name{read_abundance}
\alias{read_abundance}
\title{Reads in the abundance raw data file}
\usage{
read_abundance(path_abundance, contact_matrix)
}
\arguments{
\item{path_abundance}{path to abundance data .csv file}

\item{contact_matrix}{species contact matrix data frame}
}
\value{
a dataframe containing species names and abundances (in number of colonies)
}
\description{
Needs the contact matrix as an input, to crosscheck the species names and
their order
}
