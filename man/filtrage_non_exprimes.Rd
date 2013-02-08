\name{filtrage_non_exprimes}
\alias{filtrage_non_exprimes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filtrage_non_exprimes(mat.cdt, nbclasses, annotEch, seuil, nbval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.cdt}{ ~~Describe \code{mat.cdt} here~~ }
  \item{nbclasses}{ ~~Describe \code{nbclasses} here~~ }
  \item{annotEch}{ ~~Describe \code{annotEch} here~~ }
  \item{seuil}{ ~~Describe \code{seuil} here~~ }
  \item{nbval}{ ~~Describe \code{nbval} here~~ }
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
function (mat.cdt,nbclasses,annotEch,seuil,nbval)
{
## Fonction de filtrage #########################
### Parametres de la fonction : 
##### 1 : matrice ordonnée dans l'ordre du clustering
##### 2 : nombre de classes d'échantillon. =1 si on faire un filtrage global sur l'ensemble des valeurs
##### 3 : fichier d'annotation des échantillon. Peut-être NA si filtrage global
##### 4 : seuil choisi
##### 5 : pourcentage de valeur qui doivent passer le seuil
 if(nbclasses==1){
     classes.matrix=list(mat.cdt)
     summary(classes.matrix)
  }
  else{
    #annot_ech=read.delim(annotFile)  
    classech=annotEch
    
    classes.matrix=create_level_matrix(mat.cdt,classech)
    summary(classes.matrix)
  }
  # On regarde le nombre de valeurs dépassant le seuil dans chaque classe d'échantillon, et on garde la sonde si on a la moitié des valeurs
  current_filter= rep(FALSE, dim(classes.matrix[[1]])[1])
  
  n=nbval/100
  for(i in 1:nbclasses){
     filtre = kOverA(dim(classes.matrix[[i]])[2]*n,seuil)
     result_filter = genefilter(classes.matrix[[i]],filtre)
     current_filter= result_filter | current_filter
  }
  return(current_filter)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
