\name{LOWESS}
\alias{LOWESS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LOWESS(nom_fichier, data, pngDir, profil.median = "NA", graph = 1, projet = stop("nom de projet manquant"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nom_fichier}{ ~~Describe \code{nom_fichier} here~~ }
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{pngDir}{ ~~Describe \code{pngDir} here~~ }
  \item{profil.median}{ ~~Describe \code{profil.median} here~~ }
  \item{graph}{ ~~Describe \code{graph} here~~ }
  \item{projet}{ ~~Describe \code{projet} here~~ }
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
function(nom_fichier,data,pngDir,profil.median="NA",graph=1,projet=stop("nom de projet manquant"))
{
	nomGenes=rownames(data)
	nomEchan = colnames(data)
	
	# Calcul du profil median
	if(profil.median=="NA"){
	   profil.median = apply(data,1,median,na.rm=TRUE)
	}
	
	######### Tri par ordre croissant des profils medians  #######
	
	#On ordonne la matrice par ordre croissant des profils median
	ordre.profil = order(profil.median)
	
	#On range la matrice de depart
	matOrdonne = data[ordre.profil,]
	
	#On range le vecteur de profil median
	profilOrdonne = profil.median[ordre.profil]
	
	#On range le vecteur des noms de genes
	nomGenes = as.character(nomGenes)
	nomGenes = nomGenes[ordre.profil]	
	########################################################
	
	############# Lowess #######################	
	#Matrice lowess
	# On passe les donnees en log
	#pngDir<-"../,"02-lowess/"
	lowessCurve = NULL
	lowessCurve = cbind(lowessCurve,apply(log(matOrdonne),2,my.lowess,log(profilOrdonne),f=0.01))
	
	# Calcul de la matrice normalisee
	increm=c(1:dim(data)[2])
	matNormlog = sapply(increm, calculeNorm, log(matOrdonne), lowessCurve, log(profilOrdonne))
	matNorm=exp(matNormlog) 	
	
  profilMedNorm = apply(matNorm,1,median,na.rm=TRUE)
	
	
  if(graph==1){	
  	########## Plots ###########################
  	diagonal=c(min(profilOrdonne,na.rm=T), max(profilOrdonne,na.rm=T))
  	
  	# Graphs avant normalisation pour tous les echantillons
  	dev.set(dev.next())
  	graph=sapply(increm, traceGraph, matOrdonne, profilOrdonne, diagonal, lowessCurve,nomEchan,"1-Avant", pngDir)
		
  	# Graphs apres normalisation pour tous les echantillons
  	dev.set(dev.next())
  	graph=sapply(increm, traceGraph, matNorm, profilMedNorm, diagonal, NULL,nomEchan,"2-Apres", pngDir)
	
  	#Graph valeurs brutes/valeurs normalisees
  	dev.set(dev.next())
  	graph=sapply(increm, traceAvantApres, matOrdonne, matNorm, nomEchan, pngDir)
	} ## If graph
	##############################
	
	### Exportation des resultats	
	rownames(matNorm)=nomGenes
	colnames(matNorm) = nomEchan
	matNorm<-round(matNorm,2)
	
	nomFile = paste(pngDir,projet,"-",nom_fichier, ".normalisation.txt", sep="")
	write.table(matNorm,nomFile,sep="\t",row.names=TRUE, col.names=NA)
	print("Fonction LOEWESS!!")
  return (matNorm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
