% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svar_functions.R
\name{findQs}
\alias{findQs}
\title{findQs}
\usage{
findQs(k, f)
}
\arguments{
\item{k}{number of predictors}

\item{f}{zero restriction matrix}
}
\value{
the function returns a list which contains three elements.
        Q is a list containting the linear restrictions for each column
        index is ordered index
        flag is the status of the system, where "1" means Overidentified,
        "0" means exactly indentified and "-1" means underidentified
}
\description{
find the Q matrices
}

