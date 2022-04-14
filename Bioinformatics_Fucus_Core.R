library(tidyverse)
# Load the reshape2 package for converting between long and wide format data
library(reshape2)
# Load the stringr package for improved filtering of text
library(stringr)
# Load the ape package for reading and modifying phylogenetic trees
library(ape)
# Load the phyloseq package for microbial community analysis
library(phyloseq)
# Load the data.table package for better metadata manipulation
library(data.table)
# Load the viridis package for colour palettes for continuous data
library(viridis)
# Load the qualpalr package for colour palettes for qualitative data
library(qualpalr)
# load the ggplot2 package for visualization of data
library(ggplot2)
#load the vegan library
library(vegan)
library(dplyr)        
# filter and reformat data frames
library(tibble)      
# Needed for converting column to row names

#### environment settings ####
#set working directory
theme_set(theme_bw())
setwd("~/desktop/R_Project/Merging_RDS/downstream/")

####read phyloseq RDS from dada2###
project_data <- readRDS("fucus_core_merged_cleaned_phyloseq_new_1500.RDS")

####subset data if needed####
#project_data <- project_data %>%
  #subset_samples(project != "saccharina_hatchery") %>%
  #subset_samples(sample_type %in% c("fucus", "rock", "water"))



#########################
#### Basic alpha div plot ####
#chao1
pdf("Chao1_fuscus_merged_site_project.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, height = 9)# define plot width and height. completely up to user.




###########NMDS - Beta Diversity##############
set.seed(24)
NMDS.bray <- ordinate(
  physeq = project_data,
  method = "NMDS",
  distance = "bray"
) # you can choose different methods and distance metrics, see the ordinate help page for details. this function works with "phyloseq" class objects.

#### making beta div plots ####
#we get more plotting control if we don't use the phyloseq plotting functions for ordination plots, and instead add the results of the ordination to our existing metadata
NMDS <- as.data.frame(unclass(sample_data(project_data)))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2


pdf("NMDS_fucus.pdf"
    , width = 10 # Default is 7
    , height = 10 # Change to 10; make it taller
)
p <- ggplot(NMDS.sort, aes(x=NMDS.bray1, y=NMDS.bray2, color = sample_type)) 
# change the first argument to NMDS.sort if the optional command was ran
p + geom_point(size=3, alpha=0.85) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  labs(title="NMDS_fucus_all_betadiversity") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  xlim(-1.5, 1.5) + ylim(-1.5, 1.5)
dev.off()

#######PERMANOVA ANALAYSIS######
View(as.data.frame(unclass(sample_data(project_data))))
project_bray <- phyloseq::distance(project_data, method = "bray")
sample_df <- data.frame(sample_data(project_data))
adonis2(project_bray ~ site_project,
        data= sample_df , permutations=9999, by = "margin")




#####2018/2019 Calvert Island Fucus core#####

theme_set(theme_bw())
setwd("~/desktop/R_Project/Merging_RDS/downstream/bubble_plot/2018_2019calvert/")
###upload data
New_plots<- read.csv("Averages_for_Dot_plot_granulo_site_removed_artificials.csv")
head(New_plots)

#### Go to wide/fat format
Long_np<-gather(New_plots,Taxa,Relative_Abundance, ASV1,	ASV4,	ASV11,	ASV16,	ASV40,	ASV59,	ASV120,	ASV130,	ASV202,	ASV232,	ASV579,	ASV997,	ASV1419, ASV1423, factor_key=TRUE)
head(Long_np)
#remove Rows having NAs
Long_np <- Long_np[rowSums(is.na(Long_np)) == 0,]


#Order your factors for plot 
Long_np$Samples <- factor(Long_np$Samples, level = c("2018_calvert",	"2019_calvert",	"2017_quadra",	"2015_calvert_WB",	"2018_calvert_NB",	"2018_calvert_PB",	"2018_calvert_WB_high",	"2018_calvert_WB_low",	"2018_calvert_WB_wall",	"2019_calvert_NB",	"2019_calvert_PB",	"2019_calvert_WB_high",	"2019_calvert_WB_low",	"2015_calvert_WB_high",	"2015_calvert_WB_low"))
Long_np$Samples <- factor(Long_np$Samples, level = c("All",	"2018_calvert",	"2019_calvert",	"2017_quadra",	"2015_calvert_WB"))

###Make dot plot 
pdf("buble_plot_fucus_core_granulo_site_removed_artificials.pdf"
    , width = 11 # Default is 7
    , height = 10 # Change to 10; make it taller
)
Long_np$Relative_Abundance <- as.numeric(as.character(Long_np$Relative_Abundance))
Dot_plot_Survey_core<- ggplot(Long_np, aes(y= Taxa, x = Samples)) +geom_point(aes(size = Relative_Abundance))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_size(range=c(0,8)) +  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))  +
  theme(axis.text.x=element_text(size=12, angle= 45, hjust= 1), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=8, face="bold")) + ylab("Core Taxa") + xlab("Location") + ggtitle("")

Dot_plot_Survey_core

dev.off()

#### CORE MEMBER DETECTION - Simple Frequency Threshold Method####
#install.packages("microbiome")
#install.packages("knitr")
library(microbiome)
library(knitr)
# Transform to compositional abundances from phyloseq dataset
pseq.rel <- microbiome::transform(project_data, "compositional")

# Pick the core (e.g. >0.1% relative abundance (optional) in >50% of the samples)
a <- core(pseq.rel, detection = 0.01/100, prevalence = 50/100)
View(as.data.frame(unclass(otu_table(a))))
write.table(otu_table(a), file = "~/desktop/50prevalance.csv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)



#########Indicator species analysis (IndVal)#########
library(dplyr)
library(vegan)
library(labdsv)
library(indicspecies)
library(phyloseq)

### read metadata and OTU table
metadata <- read.csv(file="genus_level_final.sam.csv")  ### after adding columns manually
View (metadata)
otu_table<- read.csv(file="otu_rev.csv")
View (otu_table)



### subsamples from metadata using subset function
#example#metadata <- subset(metadata, project %in% "calvert_2018")

### confirm subsamples
levels(metadata$sample_type2)
metadata$sample_type2 = factor(metadata$sample_type2)
levels(metadata$sample_type2)
table(metadata$sample_type2)

master_table <- left_join(metadata, otu_table, by = "sample")
fucus <- master_table
#head(fucus)

### Set factors you are interested in i.e. want to get core for fucus but using seawater as control, then use the column with info on sample_type

## use one of the interesting info
sample_type <-fucus$sample_type2
class(sample_type)
levels(sample_type)
### Creating an object to store abundances only so you can run the analysis (i.e. remove 7 first columns of metadata with dplyr)
fucus_abund <- fucus %>% dplyr::select(-(1:8)) 
#head(fucus_abund)
#### Multipatt analysis: indval with fucus and water ####
multipatt.fucus <- multipatt(fucus_abund, sample_type, control = how(nperm=999))
summary(multipatt.fucus)
## get output and save it
indaval_output <- capture.output(summary(multipatt.fucus, indvalcomp=TRUE))
write.table(as.data.frame(indaval_output), file = "file path", quote=F, row.names=F, col.names=T, sep="n")

