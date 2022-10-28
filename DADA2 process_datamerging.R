## https://github.com/benjjneb/dada2/issues/95
## with lots of help from Bia Segovia and Geneveive Lajoie

##### WORKSPACE SET UP #####
## load packages
library(tidyverse)
library(phyloseq)
library(dada2)
#library(vegan)
#library(DECIPHER)
#library(purrr)
library(data.table)
library(zoo)


## set working directory 
theme_set(theme_bw())
setwd("~/desktop/R_Project/Merging_RDS/")


#####read Fasta (sequence info) file######
s <- readDNAStringSet("Westbeach_rep_set_16s.fasta")
View(as.data.frame(s))
write.csv(s, file="st.csv", row.names=T)

#### read csv files after manually revise OTU table: replace ASV# with real sequences ####
#### save as RDS file for merging two datasets ####
csv1 <- read.csv(file="Calvert_2018_readyfor_merge.csv")
csv2 <- read.csv(file="Calvert_2019_readyfor_merge.csv")
csv3 <- read.csv(file="Quadra_2017_readyfor_merge.csv")
csv4 <- read.csv(file="Westbeach_2015_readyfor_merge.csv")
write_rds(csv1, "Calvert_2018_readyfor_merge.RDS")
write_rds(csv2, "Calvert_2019_readyfor_merge.RDS")
write_rds(csv3, "Quadra_2017_readyfor_merge.RDS")
write_rds(csv4, "Westbeach_2015_readyfor_merge.RDS")

##### READ IN DATA #####
## load processed sequence files (2018 and 2019 Calvert seaweed data)
st1 <- readRDS("Calvert_2018_readyfor_merge.RDS")
st2 <- readRDS("Calvert_2019_readyfor_merge.RDS")
st3 <- readRDS("Quadra_2017_readyfor_merge.RDS")
st4 <- readRDS("Westbeach_2015_readyfor_merge.RDS")
#for viewing
View(as.data.frame(unclass(sample_data(test))))
View(as.data.frame(unclass(otu_table(test))))
otu <- as.data.frame(unclass(otu_table(project_data)))
asv <- colSums(otu)
View(asv)
View(as.data.frame(unclass(tax_table(test))))

## renames rows to be the sample ID
samp1 <- st1[,-1]
rownames(samp1) <- st1[,1]
samp2 <- st2[,-1]
rownames(samp2) <- st2[,1]
samp3 <- st3[,-1]
rownames(samp3) <- st3[,1]
samp4 <- st4[,-1]
rownames(samp4) <- st4[,1]

## turn into matrix to merge
samp1 = as.matrix(samp1)
samp2 = as.matrix(samp2)
samp3 = as.matrix(samp3)
samp4 = as.matrix(samp4)

##### MERGE DATASETS ####
## merge both datatsets 
colnames(samp1) <- gsub("N", "", colnames(samp1))
colnames(samp2) <- gsub("N", "", colnames(samp2))
colnames(samp3) <- gsub("N", "", colnames(samp3))
colnames(samp4) <- gsub("N", "", colnames(samp4))

#create phyloseq otu_table
otus <- otu_table(t(samp4), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) < 5) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.1) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

write_rds(seqtab.nosingletons, "Calvert_2018_readyfor_merge_filtered.RDS")
write_rds(seqtab.nosingletons, "Calvert_2019_readyfor_merge_filtered.RDS")
write_rds(seqtab.nosingletons, "Quadra_2017_readyfor_merge_filtered.RDS")
write_rds(seqtab.nosingletons, "Westbeach_2015_readyfor_merge_filtered.RDS")
st1 <- readRDS("Calvert_2018_readyfor_merge_filtered.RDS")
st2 <- readRDS("Calvert_2019_readyfor_merge_filtered.RDS")
st3 <- readRDS("Quadra_2017_readyfor_merge_filtered.RDS")
st4 <- readRDS("Westbeach_2015_readyfor_merge_filtered.RDS")
st.all <- mergeSequenceTables(st1, st2, st3, st4, tryRC=TRUE)
st.all <- mergeSequenceTables(samp1, samp2, samp3, samp4)
## remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE,verbose = FALSE)
dim(seqtab.nochim)


# Collapse ASVs that only differ by end base pair
seqtab.nochim<-collapseNoMismatch(seqtab.nochim, verbose=T) # Check Issue 716 from dada2 github
# https://github.com/benjjneb/dada2/issues/626
dim(seqtab.nochim)

## save seqtab file as RDS
write_rds(seqtab.nochim, "seqtab_nochim_seaweedcore_all.rds")


##### ASSIGN TAXONOMY #####
## SILVA needs to be downloaded in the directory to run
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/R_project/silva_for_dada2/silva_nr_v132_train_set.fa.gz", 
                       multithread=TRUE)
## add species information to taxonomy data
taxa_sp <- addSpecies(taxa, "~/Desktop/lab_member_files/siobhan/laby_seagrass/silva_species_assignment_v132.fa")

# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa),taxa),"~/Desktop/calvert_taxonomy_SILVAv132_16s_addspecies.txt", 
            row.names=FALSE, quote=F, sep="\t")

# read saved taxonomy table
taxa <- read.delim("fucuscore_taxonomy_SILVAv132_16s_addspecies.txt")
rownames(taxa) <- taxa[,1]
taxa <- taxa[,-1]

## add species information to taxonomy data
#taxa_sp <- addSpecies(taxa, "~/Desktop/lab_member_files/siobhan/laby_seagrass/silva_species_assignment_v132.fa")

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa) <- c("domain", "phylum", "class", "order", "family", "genus", "species") 


# Save taxonomy table
write.table(data.frame("row_names"=rownames(taxa),taxa),"~/Desktop/fucus_core_taxonomy_noNAs_SILVAv132_16s.txt", 
            row.names=FALSE, quote=F, sep="\t")


##### READ IN FILES #####
## read in sequence table
seqtab.nochim = readRDS("seqtab_nochim_fucuscore_all.rds")
## read in taxonomy table
tax = read.table(file="fucus_core_taxonomy_noNAs_SILVAv132_16s.txt", sep='\t', header = T)
## read in metadata
metadata = read.csv("metadata_fucuscore.csv")

##### REMOVE UNWANTED TAXA FROM TAXONOMY FILE #####
## remove chloroplasts and other non-target taxa
tax_clean = subset(tax, 
                   tax$domain !="Eukaryota" & 
                     tax$order !="Chloroplast" &
                     tax$family !="Mitochondria" &
                     tax$family !="Chloroplast" &
                     tax$domain != "Unassigned")
## save file with taxonomy that was removed 
tax_nontargets = subset(tax, 
                        tax$domain =="Eukaryota" | 
                          tax$order =="Chloroplast" |
                          tax$family =="Mitochondria" |
                          tax$family =="Chloroplast" |
                          tax$domain == "Unassigned")

## save non-target file
write_rds(tax_nontargets, "tax_nontarget_data.rds")


##### PROPAGATE THE TAXONOMY INTO EMPTY LOWER RANKS #####
# propagates taxonomy from left
tax_propagated = tax_clean %>%
  t() %>% #transpose (moves taxonomy from column names to row names)
  na.locf() %>% #fill the NAs with the values from the cell to the left (higher taxonomic rank)
  t() # transpose back to have column names be taxonomy and row names be ASV

# make into dataframe
tax_propagated = as.data.frame(tax_propagated)

# extract rownames into column
setDT(tax_propagated, keep.rownames = TRUE)[] 

# change column names
colnames(tax_propagated) <- c("asv_id","asv_sequence","Domain","Phylum","Class","Order","Family","Genus","Species") 

write.csv(tax_propagated, "tax_for_FastTree.csv")

##### MATCH ASV AND TAXAONOMY #####
## rotate asv table and make into dataframe
asv.t = t(as.matrix(seqtab.nochim))
asv.td = as.data.frame(asv.t)
# extract rownames into column
setDT(asv.td, keep.rownames = TRUE)[] 
## rename sequence column
colnames(asv.td)[grepl("rn", colnames(asv.td))] <- "asv_sequence"

## merge with taxonomy by sequence to remove asvs and taxonomy that don't have a pair in the other dataset
asvtax = inner_join(tax_propagated, asv.td)

# set row names for phyloseq
row.names(asvtax) <- asvtax$asv_id
row.names(asvtax)
head(asvtax)
## remove taxonomy and asv_id column
asv.clean = asvtax[,-c(1:9)]
head(asv.clean)

# set row names for phyloseq
row.names(tax_propagated) <- tax_propagated$asv_id
## remove extra columns
tax_propagated = tax_propagated[,-c(1:2)]
row.names(tax_propagated)
head(tax_propagated)

##### MATCH METADATA AND ASV #####
## keep only metadata with matching ASVs
asv.for.meta = t(as.matrix(asv.clean))
## make into dataframe
asv.for.meta = as.data.frame(asv.for.meta)
# extract rownames into column
setDT(asv.for.meta, keep.rownames = TRUE)[] 
## rename sequence column
colnames(asv.for.meta)[grepl("rn", colnames(asv.for.meta))] <- "sample_id"
## keep sample_id
asv.for.meta = asv.for.meta[,1]

## rename sequence column
colnames(metadata)[grepl("X", colnames(metadata))] <- "sample_id"


## inner join with metadata
metadata.clean = left_join(asv.for.meta, metadata)

## rename rows
metadata.clean =  metadata.clean %>% 
  tibble::column_to_rownames("sample_id") 



##### SAVE INDIVIDUAL FILES AS RDS ####
write_rds(asv.clean, "ASV_fucus_core.rds")
write_rds(tax_propagated, "TAXONOMY_fucus_core.rds")
write_rds(metadata.clean, "METADATA_fucus_core.rds")

asv.clean = readRDS("ASV_fucus_core.rds")
tax_propagated = readRDS("TAXONOMY_fucus_core.rds")
metadata.clean = readRDS("METADATA_fucus_core.rds")

##### CREATE PHYLOSEQ OBJECT #####
## make into matrix
otu_mat = as.matrix(asv.clean)
## the row names are the ASV DNA sequences, which matches with the ASV file's names
tax_mat = as.matrix(tax_propagated)


## transform into phyloseq-ready objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)  
samples = sample_data(metadata.clean)

## make the phyloseq object
fucus = phyloseq(OTU, TAX, samples)
fucus
rank_names(fucus)
sample_variables(fucus)

phyltax = tax_table(fucus)
print(phyltax)

head(laby@otu_table)

## save phyloseq object as RDS
write_rds(fucus, "fucus_core_merged_cleaned_phyloseq.RDS")


## go for downstream analysis for next steps
