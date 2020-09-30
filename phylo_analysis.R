library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
library("ggmsa")
library(seqinr)
seqs <- readDNAStringSet("ff", format = "fasta")


aligned <- AlignSeqs(seqs)
writeXStringSet(aligned,
                  file="aligned.fasta")
prot<- read.alignment("aligned.fasta", format = "fasta")
  
D <- dist.alignment(prot, matrix = "similarity")
temp <- as.data.frame(as.matrix(D))
  #table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)+ scale_color_viridis()
  
tre <- nj(D)
tre <- ladderize(tre)

tre.new <- tre
# change tip labels to full alignment names
tre.new$tip.label <- aligned@ranges@NAMES
  
  
ggtre1 <- ggtree(tre) + geom_tiplab(size=3)
  
  
  
ggtre2 <-msaplot(p=ggtre1,"aligned.fasta", offset=0.5, width=2)
  
 
plot(ggtre1)
    
  
plot(ggtre2, cex = 0.6)
    
  
  