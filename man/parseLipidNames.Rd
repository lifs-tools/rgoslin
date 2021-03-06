% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgoslin.R
\name{parseLipidNames}
\alias{parseLipidNames}
\title{Parse multiple lipid names and return a data frame with the results.}
\usage{
parseLipidNames(lipidNames)
}
\arguments{
\item{lipidNames}{The vector of lipid names to parse.}
}
\value{
Data frame where each row reports the parsing result of each element in lipidNames.
}
\description{
\code{parseLipidNames} reads the provided lipid names vector and returns structural information as a data frame.
Will return a cell with the "Grammar" column set to "NOT_PARSEABLE" if none of the parsers was able to parse the provided name successfully.
If any error was raised, returns an empty data frame.
}
\examples{
parseLipidNames(c("PC 32:1","LPC 34:1","TG(18:1_18:0_16:1)"))
}
