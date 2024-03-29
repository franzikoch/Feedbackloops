% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_Jacobians_Functions.R
\name{assemble_jacobian}
\alias{assemble_jacobian}
\title{Construct a community matrix}
\usage{
assemble_jacobian(interaction_table, species_list, ij_col, ji_col)
}
\arguments{
\item{interaction_table}{dataframe created by the interaction strengths function}

\item{species_list}{list of species names in the community, can be taken from the abundance table}

\item{ij_col}{name of the column that contains effects of species i on species j}

\item{ji_col}{name of the column that contains effects of species j on species j}
}
\value{
A community matrix
}
\description{
Turns interaction strengths from the tabular format (output of \emph{interaction_strengths()})
into a community matrix.
}
\details{
To do this, species names from the abundance table is used to name the
rows and columns of the matrix. The corresponding interaction strengths
are then picked from the dataframe, based on the species name.
}
