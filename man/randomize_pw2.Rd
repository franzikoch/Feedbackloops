% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Randomization_functions.R
\name{randomize_pw2}
\alias{randomize_pw2}
\title{Second version of pairwise randomisation procedure}
\usage{
randomize_pw2(df)
}
\arguments{
\item{df}{interaction table (created by interaction_strengths)}
}
\value{
the same interaction table but with two additional columns $F_ij_B_pw2
and $F_ji_B_pw2, that contain randomised interaction strengths
}
\description{
In this version, pairs are randomly reshuffled but F_ji and F_ij values cannot
be switched!! Thus values always stay on one side of the matrix diagonal, which
should avoid the creation of strong intransitive loops ? Not really well tested though
}
