### Tree analysis
library(phangorn)
setwd("~/Documents/Parasponia/Genomic_sequencing/OrthoFinder/third.run/Eurosids/Results_Jun21")

files<-list.files(path="./NJ_Trees", pattern="OG.+_tree.txt")
dir.create("NJ_Orthologies")

blacklist <-as.character(read.table("/home/robin/Documents/Parasponia/Genomic_sequencing/OrthoFinder/TorRG33x02_asm01_ann01.blacklist", sep="\t")[,1])

#Function to test if all descendants match to a vector of tips
match.descendants<- function(desc,tips){
  mean(desc %in%tips)==1    
}

#Function to get nodes comprising paralogs
get.para<- function(tree,node,tips){
  tip<-node
  while (Ancestors(node,x=tree,type="parent")!=0  && mean(unlist(Descendants(Ancestors(node,x=tree,type="parent"),x=tree,type="tips")) %in% tips)==1) {
    node<-Ancestors(node,x=tree,type="parent")
  }
  para<-unlist(Descendants(node,x=tree,type="tips"))
  para<-para[para%in%tip!=T]
  list(Node=node,Para=para)
}

#Function to translate tip numbers to names, concatenating multiple names into a single string
get.names<-function(x,names) {
  if (anyNA(x)){
    x } else {
  paste(names[x], collapse=",")
    }
}

#function to generate orthology data from a tree file
get.orthologies <- function(treefile){
t<-read.tree(file=paste("NJ_Trees",treefile, sep="/" ))
t$tip.label<-gsub("^[[:alpha:]]+_","",t$tip.label) #remove suffix
t<-drop.tip(t,blacklist)
if (length(t$tip.label)>2){
t<-multi2di(midpoint(t))}
t$tip.label<-gsub("Pan_Pan","Pan",t$tip.label)
t$tip.label<-gsub("Tor_Tor","Tor", t$tip.label)
Pan<-grep("Pan",t$tip.label)
Tor<-grep("Tor", t$tip.label)
if (length(Pan)!=0 || length(Tor)!=0) {
if (length(Pan)!=0) {
  PanPara<-lapply(Pan,get.para, tree=t, tips=Pan)
  PanOrtho<-sapply(sapply(lapply(PanPara,'[[',"Node"),FUN=Siblings,x=t),FUN=Descendants, type="tips", x=t)
  PanOrtho[sapply(PanOrtho,match.descendants,Tor)!=T]<-NA
  PanOrtho[sapply(PanOrtho,length)==0]<-NA
  PanOrtho<-sapply(PanOrtho,get.names,names=t$tip.label)
  PanPara<-sapply(lapply(PanPara,'[[',"Para"),get.names,names=t$tip.label)
} else {
PanPara=NULL
PanOrtho=NULL
}
if (length(Tor)!=0) {
  TorPara<-lapply(Tor,get.para, tree=t, tips=Tor)
  TorOrtho<-sapply(sapply(lapply(TorPara,'[[',"Node"),FUN=Siblings,x=t),FUN=Descendants, type="tips", x=t)
  TorOrtho[sapply(TorOrtho,match.descendants,Pan)!=T]<-NA
  TorOrtho[sapply(TorOrtho,length)==0]<-NA
  TorOrtho<-sapply(TorOrtho,get.names,names=t$tip.label)
  TorPara<-sapply(lapply(TorPara,'[[',"Para"),get.names,names=t$tip.label)
} else {
TorPara=NULL
TorOrtho=NULL
}

write(t(cbind(rep(gsub("_tree.+","",treefile),length(c(Pan,Tor))),t$tip.label[c(Pan,Tor)],c(PanOrtho,TorOrtho),c(PanPara,TorPara))), file=paste("NJ_Orthologies", gsub("_tree","_PanTor.orthology",treefile),sep="/" ), ncolumns=4, sep="\t")
} 
}

for (F in files)
{get.orthologies(F)}






