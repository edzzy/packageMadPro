\name{madPro}
\alias{madPro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
madPro(savingData = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{savingData}{ ~~Describe \code{savingData} here~~ }
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
function(savingData=TRUE){
require(limma)
require("genefilter")
print("Debut Setup de l'analyse")
test<- dir(pattern="pSetup.txt$")
if(length(test) == 0)
  stop("Le fichier pSetup.txt n'est pas present")
  
pSetup<-read.delim("pSetup.txt",row.names=1,header=FALSE,stringsAsFactors=FALSE)
projet=pSetup[1,]
nom_fichier <- pSetup[2,] 
#agilent=AG ou gpr=GPR nimbelgen=NG
typeArray<- pSetup[3,] 
#nombre de sondes pour la test matrice aleatoire.
alea<-pSetup[4,]
alea<-as.numeric(alea)
fileEch<-pSetup[5,]
fileGene<-pSetup[6,]
dirFile<-pSetup[7,]
nbval<-as.numeric(pSetup[8,])
nbclasses<-as.numeric(pSetup[9,])
filterParam<-as.numeric(pSetup[10,])
print ("Fin Setup de l'analyse")

test<- dir(pattern=fileEch)
answer<-paste("le fichier :",fileEch,"n'est pas present",sep=" ")
if(length(test) == 0)
  stop(answer)
###Annotation des echantillons########  
print ("début annotation des echantillons")
frameSample<-read.delim(fileEch)
namesArray<-createNamesArray(frameSample)
frameFac<-createFactor4matrix2png(frameSample,namesArray)
if(is.null(dim(frameFac))){
	frameFac<-as.data.frame(frameFac)
  }
print ("fin annotation des echantillons")
 	
namesFiles<-frameSample$nomFichiers
if(fileGene != ""){
  test<- dir(pattern=fileGene)
  answer<-paste("le fichier :",fileGene,"n'est pas present",sep=" ")
  if(length(test) == 0)
    stop(answer)
	namesGenes<-read.delim(fileGene,header=FALSE)
	namesGenes=unlist(namesGenes[,1])
  }else{
	namesGenes<-c()
  }

#######################################################
create_pathway(projet)

# on lit la matrice
if(typeArray == "GPR"){
  
  files<- dir(path=dirFile,pattern=".*\\.gpr$")
  #on change de repertoire
 
  if (all(files != namesFiles))
    stop("Non correspondance entre les noms du fichiers d'annotation et les noms reels")
    files<-paste(dirFile,files,sep="/")
	dataMA<-read_GPR(files,namesArray)
	
  }else if(typeArray=="AG"){
  
  files<- dir(path=dirFile,pattern="U.*\\.txt$")
  #on change de repertoire
  
  if (all(files == namesFiles) == FALSE)
    stop("Non correspondance entre les noms du fichiers d'annotation et les noms reels")
    files<-paste(dirFile,files,sep="/")
	 dataMA<-read_Agilent(namesArray=namesArray,files=files,namesGenes=namesGenes,dye=2,type="AG")
	 dataMA<-delete.control.in.matrix(dataMA)
	 fileName<- paste(projet,"-",nom_fichier,".rawdata.txt",sep="")
	 fileName<-paste(dirFile,fileName,sep="/")
	 	write.table(dataMA,fileName,col.names=NA,sep="\t");
  }else if(typeArray=="NG"){
       
  files<- dir(path=dirFile,pattern=".*\\.pair$")
  #on change de repertoire
  
  if (all(files != namesFiles))
    stop("Non correspondance entre les noms du fichiers d'annotation et les noms reels")
    files<-paste(dirFile,files,sep="/") 
	read_NG(namesEch,files)
  }

fd<-paste(projet,"-01-stat-descriptive/",sep="")
##connaître la taille de la matrice
print("Dimension de la matrice : ")
#print(dim(dataMA))

#######GRAPHIQUES AVANT NORMALISATION####
##creation du boxplot
print("Raw BoxPlot")

output_name<-paste(fd,projet,"-",nom_fichier,".raw.boxplot.jpeg",sep="")
create_boxplot(dataMA,output_name)

##Statistiques descriptives
print("Statistiques descriptives")
output_name = paste(fd,projet,"-",nom_fichier,".raw-stat-des.jpeg",sep="")
create_stats_des(dataMA,output_name)

##Correlation
print("Correlation")
output_name = paste(fd,projet,"-",nom_fichier,".raw-correlation.jpeg",sep="")
create_correlation_median(dataMA,output_name)

####NORMALISATION  (LOWESS)##########
print("Normalisation")
fd<-paste(projet,"-02-normalisation/",sep="")
LOWESS(nom_fichier=nom_fichier,data=dataMA,pngDir=fd,profil.median="NA",graph=1,projet=projet)
fd<-paste(fd,"images/",sep="")
imgNames<-paste(fd,projet,"-norImg.png",sep="")
#images comprenant tous les plot de normalisation 4 échantillons par image (avant apres avant/apres)
montage_cmd<-paste("montage ",fd,"*.png -tile 3x4 -geometry 100\% -background none ", imgNames,sep="")
system(montage_cmd)
#On efface la matrice de donnees brutes de la memoire
#write.table(dataMA,paste(projet,"-",nom_fichier,".txt",sep=""),sep="\t",col.names=TRUE)
#rm(dataMA)
#on recupere la matrice normalisee

#dataN = read.table(paste(nom_fichier,".lowess.txt",sep=""), sep="\t", header=T, row.names=1, comment.char="")
### verifier si le script fonctionne toujours avec des sortis gpr
dataN = read.delim(paste(projet,"-02-normalisation/",projet,"-",nom_fichier,".normalisation.txt",sep=""), sep="\t", header=TRUE, comment.char="",row.names=1)
print("Fin Normalisation")

###GRAPHIQUE APRES NORMALISATION#####
##creation du boxplot
print("Boxplot Norm")
output_name = paste(fd,projet,"-",nom_fichier,".lowess.boxplot1.jpeg",sep="")
create_boxplot(dataN,output_name)

##Statistiques descriptives
print("Statistiques descriptives Norm")
output_name = paste(fd,projet,"-",nom_fichier,".lowess-stat-des.jpeg",sep="")
create_stats_des(dataN,output_name)

##Statistiques descriptives
print("Correlation Norm")
output_name = paste(fd,projet,"-",nom_fichier,".lowess-correlation.jpeg",sep="")
create_correlation_median(dataN,output_name)


## Tout effacer
#rm(list=ls(all=TRUE))

#creation d'une matrice reduite
print("matrice aleatoire")
sampMatrix<-sampleMatrix(dataN,alea)
sampleMName=paste(projet,"-03-clusterAleatoire/",projet,"-",nom_fichier,alea,"sample.txt",sep="")
writeMatrice(sampMatrix,sampleMName)
print("debut clustering matrice aleatoire")
###Cluster######
commandCluster<-paste(" cluster -f ",sampleMName," -l  -cg m -g 1 -e 1  -m c",sep="")
system(commandCluster)
print("fin clustering matrice aleatoire")
###Filtrage######

nameCDT<-paste(projet,"-03-clusterAleatoire/",projet,"-",nom_fichier,alea,"sample.cdt",sep="")
print("filtrage matrice aleatoire")
matSamplecdt<-read.delim(nameCDT,sep="\t")
#ordonne la matrice normalisé
nG<-as.vector(matSamplecdt[,2])
nG<-nG[-1:-2]
nG<-as.character(nG)
nA<-colnames(matSamplecdt)
nA<-nA[-1:-4]
rm(matSamplecdt)
###############
sampMatrix<-sampMatrix[,nA]
sampMatrix<-sampMatrix[nG,]
#####################Utilisation de matrix2png pour le clustering de la matrice aleatoire
facteurName<-paste(projet,"-facteur.txt",sep="")
write.table(frameFac,facteurName,col.names=NA,sep="\t",row.names=TRUE);
mapName<-paste(projet,"-mapFac.txt",sep="")
mapFac<-createMap4matrix2png(frameFac)
write.table(mapFac,mapName,sep="\t",row.names=FALSE,quote=FALSE);
frameFac<-as.data.frame(frameFac)[,nA]
ffm2p<-createMatrix2png(frameFac)
frameNames<-paste(projet,"-03-clusterAleatoire/",projet,"-frameFac.txt",sep="")
write.table(ffm2p,frameNames,sep="\t",row.names=FALSE,quote=FALSE);
outfile<-paste(projet,"-03-clusterAleatoire/",projet,"-colorSample.png",sep="")
colorName<-paste("matrix2png -size 24:24  -r -s -dmap ", mapName, " -data " , frameNames ,"  > ", outfile,sep="")
system(colorName)
outfile<-paste(projet,"-03-clusterAleatoire/",projet,"-colorSample2.png",sep="")
colorName<-paste("matrix2png -size 24:24  -dmap ", mapName, " -data " , frameNames ,"  > ", outfile,sep="")
system(colorName)
##############################################
#choix du seuil de filtrage
med.sample<-apply(sampMatrix,2,median)
seuil=mean(med.sample)
#annotation du parametre choisi pour le filtrage dans le sens du clustering
#Si il n'y a qu'un parametre il faut tester (conflit data.frame et vecteur)
if(is.null(dim(frameFac))){
annotFilter<-as.factor(as.character(unlist(frameFac)))
  }else{
annotFilter<-as.factor(as.character(unlist(frameFac[filterParam,])))
  }
#filtre<-result_filter(filtrage_non_exprimes(sampMatrix,nbclasses,annotFilter,seuil,nbval)
#filtre<-as.numeric(result_filter)
result_filter<-filtrage_non_exprimes(sampMatrix,nbclasses,annotFilter,seuil,nbval)
filtre<-as.numeric(result_filter)
########Graph visualisation du filtre choix des sondes
filename=paste(projet,"-03-clusterAleatoire/",projet,"-filtre.png",sep="")
graphMmobile(filename,filtre)
######Visualisation du signal median et du seuil
profmed<-apply(sampMatrix,1,median)
profmed<-log10(profmed)
lseuil<-log10(seuil)
filename=paste(projet,"-03-clusterAleatoire/",projet,"-signalMedian.png",sep="")
graphMmobile(filename,profmed,lseuil)
#########################
####Filtrage matrice totale
print("filtrage matrice totale")
if(is.null(colnames(frameFac))){
	dataN<-dataN[,match(colnames(dataN),names(frameFac))]
  }else{
	dataN<-dataN[,match(colnames(dataN),colnames(frameFac))]
  }
seuil<-mean(apply(dataN,2,median))
result_filter=filtrage_non_exprimes(dataN,nbclasses,annotFilter,seuil,nbval)

m.filtered=dataN[result_filter == TRUE,]
filterName<-paste(projet,"-04-filtre/",projet,"-matrix_filtree.txt",sep="")
write.table(m.filtered,filterName,col.names=NA,row.names=TRUE,sep="\t")
################Clustering de la matrice totale
print("debut clustering matrice totale")
commandCluster<-paste(" cluster -f ",filterName," -l  -cg m -g 1 -e 1  -m c",sep="")
system(commandCluster)

print("fin clustering matrice totale")

filterNamecdt<-paste(projet,"-04-filtre/",projet,"-matrix_filtree.cdt",sep="")
matcdt<-read.delim(filterNamecdt,sep="\t")
#ordonne la matrice normalisé filtre selon le cdt
nG<-as.vector(matcdt[,2])
nG<-nG[-1:-2]
nG<-as.character(nG)
nA<-colnames(matcdt)
nA<-nA[-1:-4]
rm(matcdt)
m.filtered<-m.filtered[,nA]
m.filtered<-m.filtered[nG,]
frameFac<-as.data.frame(frameFac)[,nA]
#if(is.null(colnames(frameFac))){
#frameFac<-frameFac[match(colnames(m.filtered),names(frameFac))]
#}else{
#frameFac<-frameFac[,match(colnames(m.filtered),colnames(frameFac))]
#}
########Utilisation de matrix2png pour visualation des parametres
ffm2p<-createMatrix2png(frameFac)
frameNames<-paste(projet,"-04-filtre/",projet,"-frameFacF.txt",sep="")
write.table(ffm2p,frameNames,sep="\t",row.names=FALSE,quote=FALSE);
outfile<-paste(projet,"-04-filtre/",projet,"-colorSample.png",sep="")
colorName<-paste("matrix2png -size 24:24  -r -s -dmap ", mapName, " -data " , frameNames ,"  > ", outfile,sep="")
system(colorName)

outfile<-paste(projet,"-04-filtre/",projet,"-colorSample2.png",sep="")
colorName<-paste("matrix2png -size 24:24   -dmap ", mapName, " -data " , frameNames ,"  > ", outfile,sep="")

system(colorName)

##########test stat
print("test stat")
if(is.null(colnames(frameFac))){
	result_ttest<-pairRows.t.test(m.filtered,as.factor(frameFac),padj.method="none",rank=FALSE)
  }else{
	for(i in 1:1){
		result_ttest<-pairRows.t.test(m.filtered,frameFac[i,],padj.method="none",rank=FALSE)
	}
  }	
if(savingData == TRUE){
	#sauvegarde des donnees.
	#donnees brutes
	save(frameFacMA,dataMA,file=paste(projet,"-dataRaw.Rdata",sep=""))

	#donnees normalise
	save(frameFacN,dataN,sampMatrix,mapFac,file=paste(projet,"-dataNorm.Rdata",sep=""))

	#donnees filtree
	save(frameFacF,mapFac,ffm2p,m.filtered,result_ttest,file=paste(projet,"-dataFiter.Rdata",sep=""))
  }


  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
