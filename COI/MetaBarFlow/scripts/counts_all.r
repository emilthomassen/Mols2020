
library(dplyr)
library(utils)

#Set taxa to use for filtering
#noId <- c("uncultured fungus", "uncultured eukaryote","uncultured Metazoa")
#grapId <- c("Metazoa")

#Function to filter out contaminant taxa
#taxFilter <- function(taxonomy, noId=NULL){
#  if(is.null(noId))
#    stop("no input")
#  idx <- c()
#  if(!is.null(noId))
#    idx <- c(idx, unlist(sapply(noId, function(x) which(taxonomy$species == x) )))
#  return(taxonomy[-idx,])
#}

##Function to select target group
#taxFilter2 <- function(taxonomy, grapId=NULL){
#  if(is.null(grapId))
#    stop("no input")
#  idx <- c()
#  if(!is.null(grapId))
#    idx <- c(idx, unlist(sapply(grapId, function(x) which(taxonomy$kingdom != x | is.na(taxonomy$kingdom)) )))
#  return(taxonomy[-idx,])
#}


#Import OTU table and add column w. total read counts per OTU. 
otu.table <- read.delim("YOUR_DATA_PATH/DADA2_nochim.table",check.names=FALSE)
n<-ncol(otu.table)
otu.table$count=rowSums(otu.table[2:n])

colnames(otu.table)[1] <- "qseqid"

## Tax uden filter hvis jeg ikke har kontaminanter
#Import taxonomy file. 
tax <- read.delim("YOUR_DATA_PATH/classified_corrected_02_09_22.txt", stringsAsFactors = FALSE)

totalreads<-sum(otu.table$count)

#Add final otu names to the otu table
species<-data.frame(qseqid=tax$qseqid,kingdom=tax$kingdom,phylum=tax$phylum,class=tax$class,order=tax$order,family=tax$family,genus=tax$genus,dung_associated=tax$dung_associated)


otu.table_species<-merge(otu.table,species,by="qseqid")

with_id<-sum(otu.table_species$count)

k1<-sum(otu.table_species$count[which(otu.table_species$kingdom=="Fungi")])
k2<-sum(otu.table_species$count[which(otu.table_species$kingdom=="Metazoa")])
k3<-sum(otu.table_species$count[which(otu.table_species$kingdom=="Viridiplantae")])
k4<-with_id-sum(k1,k2,k3)

p1<-sum(otu.table_species$count[which(otu.table_species$phylum=="Annelida")])
p2<-sum(otu.table_species$count[which(otu.table_species$phylum=="Arthropoda")])
p3<-sum(otu.table_species$count[which(otu.table_species$phylum=="Chordata")])
p4<-sum(otu.table_species$count[which(otu.table_species$phylum=="Mollusca")])
p5<-sum(otu.table_species$count[which(otu.table_species$phylum=="Nematoda")])
p6<-sum(otu.table_species$count[which(otu.table_species$phylum=="Rotifera")])
p7<-sum(otu.table_species$count[which(otu.table_species$phylum=="Tardigrada")])
p8<-k2-sum(p1,p2,p3,p4,p5,p6,p7)

c1<-sum(otu.table_species$count[which(otu.table_species$class=="Arachnida")])
c2<-sum(otu.table_species$count[which(otu.table_species$class=="Collembola")])
c3<-sum(otu.table_species$count[which(otu.table_species$class=="Diplopoda")])
c4<-sum(otu.table_species$count[which(otu.table_species$class=="Insecta")])
c5<-sum(otu.table_species$count[which(otu.table_species$class=="Malacostraca")])
c6<-p2-sum(c1,c2,c3,c4,c5)

o1<-sum(otu.table_species$count[which(otu.table_species$order=="Blattodea")])
o2<-sum(otu.table_species$count[which(otu.table_species$order=="Coleoptera")])
o3<-sum(otu.table_species$count[which(otu.table_species$order=="Diptera")])
o4<-sum(otu.table_species$count[which(otu.table_species$order=="Hemiptera")])
o5<-sum(otu.table_species$count[which(otu.table_species$order=="Hymenoptera")])
o6<-sum(otu.table_species$count[which(otu.table_species$order=="Lepidoptera")])
o7<-sum(otu.table_species$count[which(otu.table_species$order=="Phthiraptera")])
o8<-sum(otu.table_species$count[which(otu.table_species$order=="Psocoptera")])
o9<-sum(otu.table_species$count[which(otu.table_species$order=="Thysanoptera")])
o10<-c4-sum(o1,o2,o3,o4,o5,o6,o7,o8,o9)

dass1<-sum(otu.table_species$count[which(otu.table_species$phylum=="Arthropoda" & otu.table_species$dung_associated=="yes")])
dass2<-sum(otu.table_species$count[which(otu.table_species$order=="Coleoptera" & otu.table_species$dung_associated=="yes")])
dass3<-sum(otu.table_species$count[which(otu.table_species$order=="Diptera" & otu.table_species$dung_associated=="yes")])
dass4<-dass1-sum(dass2,dass3)

level<-c(rep("total",each=2),rep("kingdoms",each=4),rep("plyla",each=8),rep("classes",each=6),rep("orders",each=10),rep("dung",each=4))
group<-c("total","with_id","Fungi","Metazoa","Viridiplantae","NA","Annelida","Arthropoda","Chordata","Mollusca","Nematoda","Rotifera","Tardigrada","Other","Arachnida","Collembola","Diplopoda","Insecta","Malacostraca","Other","Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Phthiraptera","Psocoptera","Thysanoptera","Other","Dung-ass. Arthropods","Dung-ass. Coleoptera", "Dung-ass. Diptera", "Other dung-ass.")
numbers<-c(totalreads,with_id,k1,k2,k3,k4,p1,p2,p3,p4,p5,p6,p7,p8,c1,c2,c3,c4,c5,c6,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,dass1,dass2,dass3,dass4)

numbers_df<-data.frame("Level"=level,"Group"=group,"Count"=numbers)

saveRDS(numbers_df, file="YOUR_DATA_PATH/counts_summary")

  
write.table(numbers_df, file="YOUR_DATA_PATH/counts_summary.txt", quote=FALSE, sep='\t', col.names = NA,row.names=TRUE)


