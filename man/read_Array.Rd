\name{read_Array}
\alias{read_Array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read Array}
\description{
This function return a matrice from joining median signal of
expression files
 
}
\usage{
read_Array(namesArray, files, type, namesGenes, dye = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{namesArray}{Column names : samples names }
  \item{files}{ path and names of miccro array files }
  \item{type}{type of miccro array 4 }
  \item{namesGenes}{Gene names for miccroarray : row names : gene names }
  \item{dye}{ Number of dye (mono or bicolor array}
}
%\details{
 % ~~ If necessary, more details than the description above ~~
%}
\value{
  matrix
}

\author{Edouard Hirchaud }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{
 % }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
