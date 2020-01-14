#Load data
HIST<-read.table("mer_histo.txt")

#determine error cutoff as coverage with first frequency increase (MIN)
LEN<-length(HIST[,1])
MIN<-min(which(HIST[2:LEN,2]-HIST[1:(LEN-1),2]>0))

#calculate per-base error rate based on MIN (ERR)
ERR<-sum(as.numeric(HIST[1:MIN,1]*HIST[1:MIN,2]))/sum(as.numeric(HIST[,1]*HIST[,2]))/21

#determine coverage with maximum frequency
COVP<-which.max(HIST[MIN:10000,2])+MIN-1

# estimate coverage
COVRANGE<-c(-1,0,1)+COVP
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
SCALE<-scale01(HIST[COVRANGE,2])
COV<-round((COVP+(SCALE[3]-SCALE[1])/2), digits=1)

GSIZE<-round(sum(as.numeric(HIST[MIN:LEN,1]*HIST[MIN:LEN,2]))/COV/1000000)
GSIZEMIN<-round(sum(as.numeric(HIST[MIN:LEN,1]*HIST[MIN:LEN,2]))/(COV+1)/1000000)
GSIZEMAX<-round(sum(as.numeric(HIST[MIN:LEN,1]*HIST[MIN:LEN,2]))/(COV-1)/1000000)

# Plot histogram data and info
pdf("mer_histo.pdf",width=7,height=5)
MAX<-5*COV
plot(HIST[MIN:MAX,], type="l", xlab="multiplicity", ylab="number of kmer with given multiplicity", cex=0.8)
abline(v=COV, col="red")
text(x=COV,y=0, labels=paste("estimated coverage =",COV), pos=3, offset=0.2, cex=0.8, col="red")
text(x=3*COV,y=0, labels=paste("estimated genome size = ", GSIZE, " (", GSIZEMIN,"-",GSIZEMAX, ") Mb", sep=""), 
     pos=3, offset=5, cex=0.8)
dev.off()

#write text output
writeLines(c(
paste("total kmers = ", sum(as.numeric(HIST[,1]*HIST[,2])),sep=""),
paste("error-free kmers = ", sum(as.numeric(HIST[MIN:LEN,1]*HIST[MIN:LEN,2])),sep=""), 
paste("error-threshold multiplicity = ", MIN ,sep=""), 
paste("peak multiplicity = ", COVP, sep=""), 
paste("estimated coverage = ", COV, sep=""),
paste("estimated genome size = ", GSIZE, " (", GSIZEMIN,"-",GSIZEMAX, ") Mb", sep="")  ), 
con="mer_histo.out")
