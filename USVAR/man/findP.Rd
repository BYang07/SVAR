% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svar_functions.R
\name{findP}
\alias{findP}
\title{findP}
\usage{
findP(C, B, Q, p, k, ordering)
}
\arguments{
\item{C}{initial short run impact matrix, usually from a cholesky}

\item{B}{matrix of coefficients (including intercept estimates)}

\item{Q}{a list containting the linear restrictions for each column}

\item{p}{number of lags}

\item{k}{number of dependent variabels}

\item{ordering}{original column ordering in the matrix of restrictions}
}
\value{
the function returns an orthogonal rotation matrix
}
\description{
find the orthogonal rotation matrix P
}

