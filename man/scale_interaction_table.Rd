% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling.R
\name{scale_interaction_table}
\alias{scale_interaction_table}
\title{Scale interaction strengths within an interaction table}
\usage{
scale_interaction_table(interaction_table, r_factor)
}
\arguments{
\item{interaction_table}{table of pairwise interactions}

\item{r_factor}{factor that is mulitplied with the mean interaction strength
to replace missing self-regulation values}
}
\value{
interaction table data.frame with two new columns a_ij_scaled and a_ji_scaled
}
\description{
Off-diagonal values are scaled by dividing them by the absolute value of
their corresponding diagonal element. The resulting scaled interaction strengths
are thus defined as: \eqn{\bar{a_{ij}} = \frac{a_{ij}}{|a_{ii}|}} (Neutel and Thorne
2016).
}
\details{
For each pairwise interaction, the function identifies the corresponding
self-regulation term needed to scale an interspecific interaction strength.
If the self-regulation value is missing, it is replaced by the mean
interaction strength multiplied  with a specified factor.
}
