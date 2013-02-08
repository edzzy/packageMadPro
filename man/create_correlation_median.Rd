\name{create_correlation_median}
\alias{create_correlation_median}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create_correlation_median(data, output_name, seuil = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{output_name}{ ~~Describe \code{output_name} here~~ }
  \item{seuil}{ ~~Describe \code{seuil} here~~ }
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
function (data,output_name,seuil=0.8)
##Correlation par rapport au profil median
{
  bottomarg = nchar(max(colnames(data))) #nombre de ligne pour la marge du bas
  prof_med = apply(data,1,median)
  correl=apply(data,2,cor,prof_med)
  pcol<-rep(1,ncol(data))
  badCor<-which(correl<seuil)
  pcol[badCor]=2
  jpeg(filename = output_name, width = 1300, height = 900, quality = 100, bg = "white", res = NA)
  par(mar=c(bottomarg +5,5,3,3))
  plot(correl, type="l",xlab="", ylab="Correlation par rapport au profil median", ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,xaxt="n")
  points(correl, col=pcol,xlab="",pch=19,cex=1.5, ylab="Correlation par rapport au profil median", ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,xaxt="n")
  axis(1,1:dim(data)[2],labels=colnames(data),las="2",cex.axis=1.5)
  abline(h=seuil,col="red")
  abline(v=badCor,col="green")	
  dev.off()
	if(! is.null(badCor)){
		return(names(badCor))
	}
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
