library(ape)
library("seqinr")
setwd("~/Documents/Parasponia/Genomic_sequencing/OrthoFinder/third.run/Eurosids/Results_Jun21")

files<-list.files(path="./Alignments", pattern="OG.[[:digit:]]+.fa")
dir.create("NJ_Trees")

generate_nj_tree<-function(file) {
align<-read.alignment(paste("Alignments/",file, sep=""), format="fasta")
align$nam<-gsub("^[[:alpha:]]+_","",align$nam)
dist<-dist.alignment(align)
if (length(dist)<2){
  tree<-read.tree(text = paste("(",align$nam[1],":",dist/2,",",align$nam[2],":",dist/2,");",sep=""))
} else { 
tree<-njs(dist)
}
write.tree(tree,paste("NJ_Trees/", gsub(".fa$","_tree.txt",file), sep=""))
}


for (F in files){
  generate_nj_tree(F)
}
