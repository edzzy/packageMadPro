\name{normQuantile}
\alias{normQuantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normQuantile(mat, pngDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{ ~~Describe \code{mat} here~~ }
  \item{pngDir}{ ~~Describe \code{pngDir} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(mat,pngDir){

	nomEchan = colnames(mat)
	matN<-normalizeBetweenArrays(as.matrix(mat),method="quantile")
	profMed<-apply(mat,1,median)
	profMedN<-apply(matN,1,median)	

	increm=c(1:dim(mat)[2])
########## Plots ###########################
  	diagonal=c(min(profMed,na.rm=T), max(profMed,na.rm=T))
  	
  	# Graphs avant normalisation pour tous les echantillons
  	dev.set(dev.next())
  	graph=sapply(increm, traceGraph, mat, profMed, diagonal, NULL ,nomEchan,"1-Avant", pngDir)
		
  	# Graphs apres normalisation pour tous les echantillons
  	dev.set(dev.next())
  	graph=sapply(increm, traceGraph, matN, profMedN, diagonal, NULL,nomEchan,"2-Apres", pngDir)
	
  	#Graph valeurs brutes/valeurs normalisees
  	dev.set(dev.next())
  	graph=sapply(increm, traceAvantApres, mat, matN, nomEchan, pngDir)
return(matN)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
