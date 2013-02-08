\name{pairRows.t.test}
\alias{pairRows.t.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pairRows.t.test(mat, f, padj.method = "none", rank = FALSE, pas = 50, sSctrict = 0.001, sLache = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{ ~~Describe \code{mat} here~~ }
  \item{f}{ ~~Describe \code{f} here~~ }
  \item{padj.method}{ ~~Describe \code{padj.method} here~~ }
  \item{rank}{ ~~Describe \code{rank} here~~ }
  \item{pas}{ ~~Describe \code{pas} here~~ }
  \item{sSctrict}{ ~~Describe \code{sSctrict} here~~ }
  \item{sLache}{ ~~Describe \code{sLache} here~~ }
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
function(mat,f,padj.method="none",rank=FALSE,pas=50,sSctrict=0.001,sLache=0.0001){
	f<-as.matrix(f)
	require(genefilter)
	l=levels(as.factor(f))
	print(paste("levels :",l,sep=" "))
	nbClasse= length(l)
	print(paste("nbClasse :",nbClasse,sep=" "))
	n<-2
	mat<-as.matrix(mat)
	if(rank==TRUE){
			mat<-apply(mat,1,rank)
			mat<-t(mat)
	}
	finale<-data.frame(row.names(mat))
	for(i in 1:(nbClasse-1)){

		for(j in (i+1):nbClasse){
	
			CC<-cbind(mat[,which(f==l[i])],mat[,which(f==l[j])])
			f1<-apply(mat[,which(f==l[i])],1,mean)
			f2<-apply(mat[,which(f==l[j])],1,mean)
			MC<-cbind(f1,f2)
			UDMC<-apply(as.matrix(MC),1,UD)
			ff<-c(f[which(f==l[i])],f[which(f==l[j])])
			ff<-as.factor(ff)
			tt<-rowFtests(CC,ff,var.equal=FALSE)
			pval<-tt$p.value
			pval<-p.adjust(pval,method=padj.method)
			p<-cbind(pval,UDMC)			
			comparaison<-paste("p.val",l[i],l[j],sep="_")
			filename=paste(projet,"-05-student/",projet,"-",comparaison,".png",sep="")
			namesFiles=paste(projet,"-05-student/",projet,"-p3",comparaison,".txt",sep="")
			p3<-p[which(pval<0.001)]
			ud3<-UDMC[which(pval<0.001)]
			nameUD=paste(projet,"-06-annotation/",projet,"-p3",comparaison,sep="")
			getUpDown(ud3,nameUD) #genere 2 fichiers un avec les up et un avec les downs
			names(p3)<-finale[,1][which(pval<0.001)]
			write.table(p3,namesFiles,col.names=NA,sep="\t",quote=FALSE)
			namesFiles=paste(projet,"-05-student/",projet,"-p4",comparaison,".txt",sep="")
			p4<-p[which(pval<0.0001)]
			ud4<-UDMC[which(pval<0.0001)]
			nameUD=paste(projet,"-06-annotation/",projet,"-p4",comparaison,sep="")
			getUpDown(ud4,nameUD) #genere 2 fichiers un avec les up et un avec les downs
			names(p4)<-finale[,1][which(pval<0.0001)]
			write.table(p4,namesFiles,col.names=NA,sep="\t",quote=FALSE)
			nMat<-nrow(mat)
			np3<-length(p3)
			np4<-length(p4)
			#nombre de gene attendu pour une pval <  0.001
			att3<-round(nMat*sSctrict,0)	
			#nombre de gene attendu pour une pval <  0.0001
			att4<-round(nMat*sLache,0)
			FDR3<-round((att3/np3)*100,1)
			FDR4<-round((att4/np4)*100,1)
	        fd3<-c(np3,att3,FDR3)
	        fd4<-c(np4,att4,FDR4)
			FDRn<-c("Observe","Attendu","FDR")
			matFDR<-cbind(fd3,fd4)
			row.names(matFDR)=FDRn
			filenameFDR=paste(projet,"-05-student/",projet,"-FDR",comparaison,".txt",sep="")
			write.table(matFDR,filenameFDR,col.names=NA,sep="\t",quote=FALSE)
			pv<-as.numeric(p[,1])
			graphMmobile(filename,-log10(pv),pas=pas)
			finale<-data.frame(finale,p)
			colnames(finale)[n]=comparaison
			n<-n+1
     }
  }
  comp<-colnames(finale)
  comp<-comp[-1]
  nameGene<-finale[,1]
  finale<-finale[,-1]
  finale<-as.matrix(finale)
  rownames(finale)<-nameGene
  colnames(finale)<-comp
  #rownames(finale)<-finale[,1]
  #finale<-final[,-1]
  #finale<-as.matrix(finale)
  namesFiles<-"";
	for(i in 1:nbClasse){
		namesFiles<-paste(namesFiles,l[i],sep="-")
	}
  namesFiles<-paste(projet,"-05-student/",projet,namesFiles,"-student.txt",sep="")
  write.table(finale,namesFiles,col.names=NA,sep="\t",quote=FALSE)
  return(finale)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
