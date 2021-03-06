\name{create_stats_des}
\alias{create_stats_des}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create_stats_des(data, output_name, minMax = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{output_name}{ ~~Describe \code{output_name} here~~ }
  \item{minMax}{ ~~Describe \code{minMax} here~~ }
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
function (data,output_name,minMax=TRUE){
##Statistiques descriptives
##DESCRIPTION
#fonction qui tracent les valeurs particulieres pour chaque echantillons d'une matrice d'expression
#dans un fichier image (png)
##ARG
#data est une matrice de donnees d'expression de type data.frame ou matrix mais dont le type dois etre numerique ou double.
#output est le nom du fichier de sortie de l'image : chaine de caractere
#minMax est un boolen qui permet ou non de tracer les valeurs max et min ainsi que les 10eme plus grande et plus petite valeur
#par defaut minMax est TRUE
##RETURN
#cette fonction ne retourne aucune valeur.

    #Je tri les valeurs par colonnes en ordre croissant
    dataSortedBySample=apply(data,2,sort)
    nbgenes = dim(data)[1]
    nbech = dim(data)[2]
    bottomarg = nchar(max(colnames(data))) #nombre de ligne pour la marge du bas
    #1er, 10e ... plus grande valeur
    v1=dataSortedBySample[nbgenes,]
    v10=dataSortedBySample[nbgenes-9,]
    v100=dataSortedBySample[nbgenes-99,]
    v1000=dataSortedBySample[nbgenes-999,]
    if(nbgenes>9999){v10000=dataSortedBySample[nbgenes-9999,]}
    if(nbgenes>19999){v20000=dataSortedBySample[nbgenes-19999,]}
    v1p=dataSortedBySample[1,]
    v10p=dataSortedBySample[1+9,]
    v100p=dataSortedBySample[1+99,]
    v1000p=dataSortedBySample[+999,]
    if(nbgenes>9999){v10000p=dataSortedBySample[1+9999,]}
    if(nbgenes>19999){v20000p=dataSortedBySample[1+19999,]}

    x = c(1,nbech+1) #x pour le plot invisible
    rangev = range(v1,v1p) #pour y : 'range' returns a vector containing the minimum and maximum of all the given arguments.
    max = max(v1)
    milieu = max/2 #milieu de l'axe des y 	

    
    jpeg(filename = output_name, width = 1300, height = 900, quality = 100, bg = "white", res = NA)
    par(mar=c(bottomarg +5,5,3,25))
    plot(x,rangev, log="y", type="n",xaxt="n", xlab="", ylab="Intensity",cex.lab=1.5,cex.axis=1.5)
    # type="n" '"n"' for no plotting
	if(minMax){
		lines(v1, col="black")
		lines(v10, col="red")
		legend=c("valeur max","10eme plus grande valeur","100eme plus grande valeur","1000eme plus grande valeur")
		col=c("black","red","blue","green")
	}else{
		legend=c("100eme plus grande valeur","1000eme plus grande valeur")
		col=c("blue","green")
	}


	lines(v100, col="blue")
    lines(v1000, col="green")
    lty<- rep(1,4)
    if(exists("v10000")){
    	lines(v10000, col="orange")
    	legendtmp<-c("10000 plus grande valeur")
    	coltmp<-c("orange")
    	ltytmp<-1
    
    	legend<-c(legend,legendtmp)
    	col<-c(col,coltmp)
    	lty<-c(lty,ltytmp)
    }
    if(exists("v20000")){
    	lines(v20000, col="red")
    	legendtmp<-c("20000 plus grande valeur")
    	coltmp<-c("red")
    	ltytmp<-1
    
    	legend<-c(legend,legendtmp)
    	col<-c(col,coltmp)
    	lty<-c(lty,ltytmp)
    }
    
	if(minMax){
		lines(v1p, col="blue",lty=5)
		lines(v10p, col="orange",lty=5)
	}
	
    lines(v100p, col="green",lty=5)
    lines(v1000p, col="red",lty=5)
    legendp<-NULL	
    if(exists("v20000p")){
    	lines(v20000p, col="grey",lty=5)
    	legendp<-c("20000eme petite valeur")
        coltmp<-c("grey")
        ltytmp<-5                                            
       col<-c(col,coltmp)
       lty<-c(lty,ltytmp)
    }

    if(exists("v10000p")){
      	lines(v10000p, col="purple",lty=5)
      	legendp<-c(legendp,"10000eme petite valeur")
        coltmp<-c("purple")
        ltytmp<-5
                                                 
       col<-c(col,coltmp)
       lty<-c(lty,ltytmp)
	 }
	if(minMax){
		colp<-c("red","green","orange","blue")
		legendp<-c(legendp,"1000eme petite valeur", "100eme petite valeur", "10eme petite valeur", "valeur min")	
	}else{
		colp<-c("red","green")
		legendp<-c(legendp,"1000eme petite valeur", "100eme petite valeur")	
	}
    legend<-c(legend,legendp)
    col<-c(col,colp)
    ltyp<-rep(5,4)
    lty<-c(lty,ltyp)
	 
    axis(1,1:nbech,labels=colnames(data),las="2",cex.axis=1.5)
    par(xpd=TRUE)	
    
    legend(nbech+1.5,milieu, yjust=1, legend=legend,lwd=1,cex=1.5,col=col,lty=lty)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
