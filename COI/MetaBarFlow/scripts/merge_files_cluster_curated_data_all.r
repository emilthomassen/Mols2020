
library(dplyr)
library(utils)

#Set taxa to use for filtering
noId <- c("uncultured fungus", "uncultured eukaryote","uncultured Metazoa", "unidentified Nematode")
grapId <- c("Metazoa")

#Function to filter out contaminant taxa
taxFilter <- function(taxonomy, noId=NULL){
  if(is.null(noId))
    stop("no input")
  idx <- c()
  if(!is.null(noId))
    idx <- c(idx, unlist(sapply(noId, function(x) which(taxonomy$species == x) )))
  return(taxonomy[-idx,])
}

#Function to select target group
taxFilter3 <- function(taxonomy, grapId=NULL){
  if(is.null(grapId))
    stop("no input")
  idx <- c()
  if(!is.null(grapId))
    idx <- c(idx, unlist(sapply(grapId, function(x) which(taxonomy$kingdom != x | is.na(taxonomy$kingdom)) )))
  return(taxonomy[-idx,])
}


#Import OTU table and add column w. total read counts per OTU. 
otu.table <- read.delim("YOUR_DATA_PATH/DADA2_nochim.table",check.names=FALSE)
n<-ncol(otu.table)
otu.table$count=rowSums(otu.table[2:n])

colnames(otu.table)[1] <- "qseqid"

## Tax uden filter hvis jeg ikke har kontaminanter
#Import taxonomy file. 
taxUnfiltered <- read.delim("YOUR_DATA_PATH/classified_corrected_02_09_22.txt", stringsAsFactors = FALSE)

#Create final.id column (Temporary!!!!!!!!!!)

#final.id <- taxUnfiltered$species
#taxUnfiltered<-cbind(taxUnfiltered, final.id, stringsAsFactors=FALSE)

#Remove contaminant taxa
tax <- taxFilter(taxUnfiltered, noId)

#Select chosen group
tax <- taxFilter3(tax,grapId)


#Make a new column showing whether the final taxonomic id is obtained
tax <- tax %>% mutate(tax.id=if_else(is.na(tax$score.id)==TRUE,"no","yes")) 

#Make a new column with a final otu name, which will be the species name if the identification was to species level, and otherwise will be the qseqid

tax <- tax %>% mutate(final.otu=if_else(tax.id=="yes",score.id,qseqid))

#Remove NAs temporary
tax[is.na(tax)==TRUE]="NOCLASS"

#Correct taxonomy to only valid levels
no=0 
s=0 
g=0 
f=0
o=0 
c=0 
p=0 
k=0

for (i in 1:nrow(tax)){
  if(tax$final.otu[i]==tax$kingdom[i]){
    k<-k+1
    tax$phylum[i]=paste(as.character(tax$kingdom[i]),"phylum",as.character(tax$qseqid[i]))
    tax$class[i]=paste(as.character(tax$kingdom[i]),"class",as.character(tax$qseqid[i]))
    tax$order[i]=paste(as.character(tax$kingdom[i]),"order",as.character(tax$qseqid[i]))
    tax$family[i]=paste(as.character(tax$kingdom[i]),"family",as.character(tax$qseqid[i]))
    tax$genus[i]=paste(as.character(tax$kingdom[i]),"genus",as.character(tax$qseqid[i]))
    tax$species[i]=paste(as.character(tax$kingdom[i]),"sp.",as.character(tax$qseqid[i]))
  }
  else if (tax$final.otu[i]==tax$phylum[i]) {
    p<-p+1
    tax$class[i]=paste(as.character(tax$phylum[i]),"class",as.character(tax$qseqid[i]))
    tax$order[i]=paste(as.character(tax$phylum[i]),"order",as.character(tax$qseqid[i]))
    tax$family[i]=paste(as.character(tax$phylum[i]),"family",as.character(tax$qseqid[i]))
    tax$genus[i]=paste(as.character(tax$phylum[i]),"genus",as.character(tax$qseqid[i]))
    tax$species[i]=paste(as.character(tax$phylum[i]),"sp.",as.character(tax$qseqid[i]))
  }
  else if (tax$final.otu[i]==tax$class[i]) {
    c<-c+1
    tax$order[i]=paste(as.character(tax$class[i]),"order",as.character(tax$qseqid[i]))
    tax$family[i]=paste(as.character(tax$class[i]),"family",as.character(tax$qseqid[i]))
    tax$genus[i]=paste(as.character(tax$class[i]),"genus",as.character(tax$qseqid[i]))
    tax$species[i]=paste(as.character(tax$class[i]),"sp.",as.character(tax$qseqid[i]))
  }
  else if (tax$final.otu[i]==tax$order[i]) {
    o<-o+1
    tax$family[i]=paste(as.character(tax$order[i]),"family",as.character(tax$qseqid[i]))
    tax$genus[i]=paste(as.character(tax$order[i]),"genus",as.character(tax$qseqid[i]))
    tax$species[i]=paste(as.character(tax$order[i]),"sp.",as.character(tax$qseqid[i]))
  }
  else if (tax$final.otu[i]==tax$family[i]) {
    f<-f+1
    tax$genus[i]=paste(as.character(tax$family[i]),"genus",as.character(tax$qseqid[i]))
    tax$species[i]=paste(as.character(tax$family[i]),"sp.",as.character(tax$qseqid[i]))
  }
  else if (tax$final.otu[i]==tax$genus[i]) {
    g<-g+1
    tax$species[i]=paste(tax$genus[i], "sp.", as.character(tax$qseqid[i]))
  }
  else if (tax$final.otu[i]==tax$species[i]) {
    s<-s+1
  }
  else {
    no<-no+1
    tax$kingdom[i]=tax$qseqid[i]
    tax$phylum[i]=tax$qseqid[i]
    tax$class[i]=tax$qseqid[i]
    tax$order[i]=tax$qseqid[i]
    tax$family[i]=tax$qseqid[i]
    tax$genus[i]=tax$qseqid[i]
    tax$species[i]=tax$qseqid[i]
  }
  
}

df_summary<-data.frame("species_level"=s, "genus_level"=g, "family_level"=f, "order_level"=o, "class_level"=c, "phylum_level"=p, "kingdom_level"=k, "no_classification"=no)
write.table(df_summary, file="YOUR_DATA_PATH/summary_classifications_all.txt", sep = "\t")

#Assign final.id as final.otu

tax$final.otu<-tax$final.id




#Add final otu names to the otu table
species<-data.frame(qseqid=tax$qseqid,final.otu=tax$final.id)

otu.table_species<-merge(otu.table,species,by="qseqid")

otu.table_species$qseqid<-as.factor(otu.table_species$qseqid)

#Aggregate based upon final.otu

#Species
 
  reads <- aggregate( . ~ final.otu, data = otu.table_species, sum)
  #Now reads$qseqid is useless, as the ids have been "summed" 
  reads$qseqid<-NULL

    
  #Determine maximum sequence similarity for each final otu and add this to the taxonomy file (after removing "duplicate" rows w. #same final otu name).
  best.id <- aggregate(pident ~ final.otu, data = tax, max)
  names(best.id)[names(best.id) == 'pident'] <- 'max_id'
  tax_uniq<-distinct(tax, final.otu, .keep_all = TRUE)
  tax_uniq_maxid<-merge(best.id,tax_uniq,by="final.otu")
  
  #Make obitab file from taxonomy file. NB! id is just from the top row of the classifid file, so does not necessarily match the max_id value.
  obitab<-data.frame(id=tax_uniq_maxid$qseqid,
                     "best_identity:ncbi"=tax_uniq_maxid$max_id,
                     kingdom_name=tax_uniq_maxid$kingdom,
                     phylum_name=tax_uniq_maxid$phylum,
                     class_name=tax_uniq_maxid$class,
                     order_name=tax_uniq_maxid$order,
                     family_name=tax_uniq_maxid$family,
                     genus_name=tax_uniq_maxid$genus,
                     species_name=tax_uniq_maxid$species,
                     scientific_name=tax_uniq_maxid$final.otu,
                     final.otu=tax_uniq_maxid$final.otu,
                     dung.associated=tax_uniq_maxid$dung_associated,
                     type=tax_uniq_maxid$type,
                     check.names=FALSE)
  
  merged<-merge(reads,obitab,by="final.otu")
    
  
write.table(merged, file="YOUR_DATA_PATH/merged_table_all.txt", quote=FALSE, sep='\t', col.names = NA,row.names=TRUE)



