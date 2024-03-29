% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Randomization_functions.R
\name{randomize_minimal}
\alias{randomize_minimal}
\title{Minimal randomisation procedure}
\usage{
randomize_minimal(df, ij_col, ji_col)
}
\arguments{
\item{df}{interaction table (created by interaction_strenghts())}

\item{ij_col}{column of \eqn{a_{ij}} values to randomise (choose scaled or unscaled)}

\item{ji_col}{column of \eqn{a_{ji}} values to randomise (scaled or unscaled)}
}
\value{
Interaction table with two additional columns \emph{$a_ij_min} and \emph{$a_ji_min},
containing randomised interaction strengths
}
\description{
In this version, pairs are kept together and their above-below diagonal orientation
is also preserved.
}
\details{
In the interaction table this means the following: \eqn{a_{ij}} and \eqn{a_{ji}}
values that appear in the same row in the original table, will also be in the same
row in the pairwise randomised table. In contrast to the pairwise (weak) randomisation,
they cannot switch from  \eqn{a_{ij}} to \eqn{a_{ji}} or vice versa. Only the order of rows
in the table is randomised.

The function returns a new interaction table that contains two new columns $a_ij_min and
$a_ji_min that contains the same values as the original columns but in a different order.

To get a minimally randomised Jacobian matrix, use assemble_jacobian() and specify
the new columns.
}
