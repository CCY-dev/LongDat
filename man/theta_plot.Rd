% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theta_plot.R
\name{theta_plot}
\alias{theta_plot}
\title{Plot theta values of negative binomial models versus non-zero count
for count data}
\arguments{
\item{input}{A data frame with the first column as "Individual"
and all the columns of dependent variables (features, e.g. bacteria)
        at the end of the table.}

\item{test_var}{The name of the independent variable you are testing for,
should be a character vector (e.g. c("Time"))
       identical to its column name and make sure there is no space in it.}

\item{variable_col}{The column number of the position where the dependent
variable columns (e.g. bacteria) start in the table}

\item{fac_var}{The column numbers of the position where the columns
that aren't numerical  (e.g. characters, categorical numbers,
ordinal numbers), should be a numerical vector (e.g. c(1, 2, 5:7))}

\item{not_used}{The column position of the columns not are irrelevant and
can be ignored when in the analysis.
       This should be a number vector, and the default is NULL.}

\item{point_size}{The point size for plotting in 'ggplot2'. The default is 1.}

\item{x_interval_value}{The interval value for tick marks on x-axis.
The default is 5.}

\item{y_interval_value}{The interval value for tick marks on y-axis.
The default is 5.}

\item{verbose}{A boolean vector indicating whether to print detailed
message.
The default is TRUE.}
}
\value{
a 'ggplot' object
}
\description{
Plot theta values of negative binomial models versus non-zero count
for count data
}
\details{
This function outputs a plot that facilitates the setting of theta_cutoff
in longdat_disc() and longdat_cont(). This only applies when the dependent
variables are count data. longdat_disc() and longdat_cont() implements
negative binomial (NB) model for count data,
and if the theta (dispersion parameter) of
NB model gets too high, then the p value of it will be extremely low
regardless of whether there is real significance
or not. Therefore, the highest threshold of theta value is set and any
 variable beyond the threshold will be excluded
from the test. The default value of theta_cutoff is set to 2^20 from the
observation that 2^20 is a clear cutoff line
for several datasets. Users can change theta_cutoff value to fit
their own data.
}
\examples{
test_theta_plot <- theta_plot(input = LongDat_disc_master_table,
 test_var = "Time_point", variable_col = 7, fac_var = c(1:3))
}
