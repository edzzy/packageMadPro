\name{graphMmobile}
\alias{graphMmobile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graphMmobile(filename, value, seuil = NULL, pas = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ ~~Describe \code{filename} here~~ }
  \item{value}{ ~~Describe \code{value} here~~ }
  \item{seuil}{ ~~Describe \code{seuil} here~~ }
  \item{pas}{ ~~Describe \code{pas} here~~ }
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
function(filename,value,seuil=NULL,pas=50){

	curveMobile<-rep(1/pas,pas)
	fil<-filter(value,curveMobile)
	png(filename=filename,width = 1300, height = 900, bg = "white", res = NA )
	par(mar=c(0,0,0,0),oma=c(0,0,0,0))
	plot(value,pch=19,cex=0.8,cex.lab=1.5,cex.axis=3,col="slateblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
	lines(fil,lwd=5,col="orange2")
	axis(side=2,tcl=0.5,label=TRUE,pos=c(0,0),cex.axis=3)
	if(is.null(seuil)){
	abline(h=seuil,lwd=2,pch=19,cex=0.8,cex.lab=1.5,cex.axis=1.5,col="red")
	}
	dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
