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


#set working directory
theme_set(theme_bw())
setwd("~/desktop/R_Project/Merging_RDS/downstream/")

###Read phyloseq RDS files to work
project_data <- readRDS("fucus_core_merged_cleaned_phyloseq_final.RDS") 

#### rarefy data ####
set.seed(24) 
project_data.rarefied <- rarefy_even_depth(project_data, sample.size = 1500)
### Save rarifying data###
saveRDS(project_data.rarefied, "fucus_core_merged_cleaned_phyloseq_final_1500.RDS")

project_data <- readRDS("fucus_core_merged_cleaned_phyloseq_final_1500.RDS") #saved file after rarefaction for NMDS/PERMANOVA

####subset data if needed####
# subset data - NMDS A
project_data <- project_data %>%
  subset_samples(project != "saccharina_hatchery") %>%
  subset_samples(sample_type %in% c("fucus", "rock", "water"))

# subset data - NMDS B
project_data <- project_data %>%
  subset_samples(project != "saccharina_hatchery") %>%
  subset_samples(sample_type %in% c("fucus")) %>%
  subset_samples(site != "WB")

# subset data - only FUCUS
project_data <- project_data %>%
  subset_samples(sample_type %in% c("fucus"))

# subset data - only Quadra
project_data <- project_data %>%
  subset_samples(project %in% c("quadra_2017"))

#for viewing
View(as.data.frame(unclass(sample_data(project_data))))
View(as.data.frame(unclass(otu_table(project_data))))
otu <- as.data.frame(unclass(otu_table(project_data)))
asv <- colSums(otu)
View(asv)
View(as.data.frame(unclass(tax_table(project_data))))

##count taxa function##
ntaxa(project_data)


#### Basic alpha div plot ####
#chao1
pdf("Chao1_fuscus_merged_site_project.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, height = 9)# define plot width and height. completely up to user.

#### basic alpha div plot ####
#### step1: calculate alpha diversity and add it as a column in the metadata ####
project_data.chao1 <- estimate_richness(project_data, split = TRUE, measures = c("Chao1")) #estimate richness
sample_data(project_data)$chao1 <- project_data.chao1$Chao1 #add to metadata (the rows are in the same order already)
sample_data(project_data)$chao1 <- as.numeric(sample_data(project_data)$chao1)
#### step 2: use calculated alpha diversity to make basic plot with ggplot ####
#this plot lets you customize things a bit more than the plot_richness function, if desired
#p <- ggplot(sample_data(project_data), aes(x = factor(project, level = c('aboral', 'oral', 'cecum', 'Rock_surface', 'Seawater')), y=chao1, color=site))
p <- ggplot(sample_data(project_data), aes(x=site_project, y=chao1, colour = "grey"))
p + geom_boxplot(outlier.size=-1) + 
  #facet_grid(~ FACTOR_2, drop=TRUE, scales="free", space="free") + #drop, scales, and space are used to ensure that each boxplot bar is equally sized, and that empty levels created by the facet grid are dropped from the final plot
  labs(title="Alpha Diversity (Chao1)", x="site", y="Chao1") + 
  theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.text=element_text(size=15)) +
  scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()


which(is.na(tax_table(project_data)) == TRUE)
project_data[is.na(project_data)] <- 0


###########Beta Diversity##############
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
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$project, NMDS$sample_type),]

## manually designate colours 15
cbPalette <- c("#ca2700", "#972a3a", "#ff5e37", "#d1ff9b", "#95f926", "#58e72a", "#b4943e", "#4500aa", "#7835d3", "#ffb3cf", "#d9a8d1", "#8b30bb", "#933b8c", "#131212", "#7fdcf6")

###############NMDS PLOT2: Sample Type###############
pdf("NMDS_fucus_rock_water2.pdf"
    , width = 10 # Default is 7
    , height = 10 # Change to 10; make it taller
)
p <- ggplot(NMDS.sort, aes(x=NMDS.bray1, y=NMDS.bray2, color = sample_type)) 
# change the first argument to NMDS.sort if the optional command was ran
p + geom_point(size=3, alpha=0.85) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  labs(title="NMDS_fucus_all_betadiversity") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  xlim(-1.5, 0) + ylim(-0.8, 0.8)
dev.off()


####### START PERMANOVA ANALAYSIS######
View(as.data.frame(unclass(sample_data(project_data))))

project_bray <- phyloseq::distance(project_data, method = "bray")
project_jaccard <- phyloseq::distance(project_data, method = "jaccard")
project_wunifrac <- phyloseq::distance(project_data, method = "wunifrac")

sample_df <- data.frame(sample_data(project_data))


adonis(project_bray ~ Type, data=sample_df, method="bray")
adonis(project_jaccard ~ Type, data=sample_df, method="jaccard")
adonis(project_wunifrac ~ Type, data=sample_df, method="wunifrac")

adonis2(project_bray ~ site_project,
        data= sample_df , permutations=9999, by = "margin")
adonis2(project_jaccard ~ Type,
        data= sample_df , permutations=9999, by = "margin")
adonis2(project_wunifrac ~ Type,
        data= sample_df , permutations=9999, by = "margin")

res.adonis.rarefied <- adonis(project_bray ~ Type, data=sample_df, method="bray")
#res.adonis.rarefied <- adonis(project_bray.rarefied ~ Month, data=sample_df, method="bray")
res.adonis.rarefied #run for full description of results
summary(res.adonis.rarefied) #run for results summary (this is less informative if I remember correctly)

#install adnois.pair
#devtools::install_github("hadley/devtools")
#install.packages("devtools")
library(devtools)
library(EcolUtils)


## running 10,000 permutations instead of the default 999###
adonis2(project_bray ~ Type,
        data= sample_df , permutations=999, by = "margin")
options(scipen=999)

pair_type <- adonis.pair(project_bray, sample_df$site_project,
                         nper = 999, corr.method = "BH")

pair_type2 <- adonis.pair(project_bray, sample_df$month,
                          nper = 999, corr.method = "BH")

capture.output(file="adnois_pair_type2.txt", pair_type2)

#### analysis: beta dispersion test ####
beta.Type <- betadisper(project_bray, sample_df$Type) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
permutest(beta.Type) #do permutation test of above #this is more customizable, see documentation


beta.Month <- betadisper(project_bray.rarefied, sample_df$Month)
permutest(beta.Month)

#### saving output of these tests ####
# you can use the "capture.output()" function to save the output of any of the testing and summary commands above, like so
capture.output(file="beta_disp.output.txt", permutest(beta.Phylum))
capture.output(file="beta_disp.output.txt", permutest(beta.Order))
# where the form of the command is always:
capture.output(file="name_of_output_file", command)



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



#########################
#########Indival#########

### INDVAL Indicator Species Analysis ###

library(dplyr)
library(vegan)
library(labdsv)
library(indicspecies)
library(phyloseq)


#write metadata from RDS file
#metadata <- as.data.frame(unclass(sample_data(project_data)))
#write.csv(metadata, file="metadata.csv", row.names=T)

trans<- read.csv(file="genus_level_final.otu.csv")
trans <- t(trans)
View (trans)
write.csv(trans, file="otu_rev.csv", row.names=T)

### from HERE after revising metadata####
metadata <- read.csv(file="genus_level_final.sam.csv")  ### after adding columns manually
View (metadata)
otu_table<- read.csv(file="otu_rev.csv")
View (otu_table)



### FUCUS all data
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
### FUCUS all data w/o Quadra
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, project != "quadra_2017")
### FUCUS Quadra
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, project %in% "quadra_2017")
## select April, May, June
metadata <- subset(metadata, month %in% c("April", "May", "June"))

metadata <- read.csv(file="genus_level_final.sam.csv")
######SPECIFIC##### FUCUS WestBeach-low 2015
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% "WB_low_2015")
######SPECIFIC##### FUCUS WestBeach-high 2015
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% "WB_high_2015")
######SPECIFIC##### FUCUS WestBeach-low 2018
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% c("WB", "WB_low_2018"))
######SPECIFIC##### FUCUS WestBeach-high 2018
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% c("WB", "WB_high_2018"))
######SPECIFIC##### FUCUS Pruth 2018
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% "PB_2018")
######SPECIFIC##### FUCUS NorthBeach 2018
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% "NB_2018")
######SPECIFIC##### FUCUS WestWall 2018
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% c("WB", "WB_wall_2018"))

metadata <- read.csv(file="genus_level_final.sam.csv")
######SPECIFIC##### FUCUS WestBeach-low 2019
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% c("WB_low_2019"))
######SPECIFIC##### FUCUS WestBeach-high 2019
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% c("WB_high_2019"))
######SPECIFIC##### FUCUS Pruth 2019
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% "PB_2019")
######SPECIFIC##### FUCUS NorthBeach 2019
metadata <- subset(metadata, sample_type2 %in% c("fucus", "env"))
metadata <- subset(metadata, site_project %in% "NB_2019")


metadata <- read.csv(file="genus_level_final.sam.csv")
######SPECIFIC##### FUCUS WestBeach-low 2015
metadata <- subset(metadata, site_project %in% c("WB_low_2015", "WB_high_2015"))
metadata <- subset(metadata, sample_type3 != "env")



levels(metadata$sample_type2)
metadata$sample_type2 = factor(metadata$sample_type2)
levels(metadata$sample_type2)
table(metadata$sample_type2)

master_table <- left_join(metadata, otu_table, by = "sample")
View(as.data.frame(master_table))
write.csv(master_table, file="genus_prevalent.csv", row.names=T)
## remove NA##
#na.omit(master_table)
#na.omit(master_table)
## replace NA with 0##
#master_table[is.na(master_table)] <- 0


### Loading your data (this should contain metadata + taxa abundances)
#fucus <- read.table("your_sample(rows)_taxa(columns)_table.txt", header = T)
fucus <- master_table
#head(fucus)

### Set factors you are interested in i.e. want to get core for fucus but using seawater as control, then use the column with info on sample_type

## use one of the interesting info
sample_type <-fucus$sample_type2

class(sample_type)
levels(sample_type)

### Creating an object to store abundances only so you can run the analysis (i.e. remove 7 first columns of metadata with dplyr)
fucus_abund <- fucus %>% dplyr::select(-(1:10)) 
head(fucus_abund)

#### Multipatt analysis: indval with fucus and water ####
multipatt.fucus <- multipatt(fucus_abund, sample_type, control = how(nperm=999))
summary(multipatt.fucus)

## get output and save it
indaval_output <- capture.output(summary(multipatt.fucus, indvalcomp=TRUE))
write.table(as.data.frame(indaval_output), file = "~/desktop/R_Project/Merging_RDS/downstream/indival/Indval_new/Fucus_Quadra_APR-JUN.txt", quote=F, row.names=F, col.names=T, sep="n")

#########################
#########################

#### CORE MEMBER DETECTION by a simple frequency threshold####
#install.packages("microbiome")
#install.packages("knitr")
library(microbiome)
library(knitr)
library(RColorBrewer)
library(reshape)
library(hrbrthemes)
textcol <- "black"

project_data <- readRDS("fucus_core_merged_cleaned_phyloseq_new2_1500.RDS")
# subset data by location
project_data <- project_data %>%
  subset_samples(project != "quadra_2017") %>%
  subset_samples(sample_type %in% c("fucus"))

project_data <- project_data %>%
  # subset_samples(project == "quadra_2017") %>%
  subset_samples(sample_type %in% c("fucus"))

View(as.data.frame(unclass(sample_data(project_data))))
# Transform to compositional abundances
pseq.rel <- microbiome::transform(project_data, "compositional")
pseq.rel <- aggregate_taxa(pseq.rel, "Genus")

#SAVE
#genus_aggro <- as.data.frame(unclass(otu_table(a)))
write.csv(a, file="taxa_aggrogated.csv", row.names=T)

trans<- read.csv(file="taxa_aggrogated.csv")
trans <- t(trans)
View (trans)
write.csv(trans, file="taxa_aggrogated_trans.csv", row.names=T)

# Pick the core (>0.1% relative abundance in >50% of the samples)
a <- core(pseq.rel, detection = 0.001/100, prevalence = 0.001/100)
a <- core(pseq.rel, detection = 0, prevalence = 50/100)
View(as.data.frame(unclass(otu_table(a))))
#View(as.data.frame(unclass(sample_data(a))))


pdf("core_heatmap_all_fucus_quadra.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 7, height = 9)# define plot width and height. completely up to user.
# Core with absolute counts and horizontal view:
# and minimum population prevalence (given as percentage)
detections <- round(10^seq(log10(0.0001), log10(0.1), length = 50), 10)
p <- plot_core(a, plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .5, horizontal = FALSE) +
  theme_ipsum() +
  theme_grey(base_size=10)+
  theme(legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol),
        axis.text.y=element_text(size=10, vjust=0.2,colour=textcol, face="bold"),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))

print(p)

dev.off()


##scatter plot
theme_set(theme_bw())
setwd("~/desktop/R_Project/Merging_RDS/downstream/")

library(ggplot2)
library(hrbrthemes)
cbPalette <- c("#d5532c", "#6394d3", "grey40", "grey92","#436422", "black")
cbPalette <- c("#436422", "black", "grey")
# mtcars dataset is natively available in R
# head(mtcars)
iris<- read.csv(file="fucus_prevalence_indvaldetect_10frequency.csv")
iris<- read.csv(file="fucus_pre_multi_indval_abund.csv")
iris<- read.csv(file="fucus_pre_multi_indval_abund_overall_indval.csv")
head(iris)
# A basic scatterplot with color depending on Species

#make the order reverse to show dot visible
iris <- iris[order(iris$core, decreasing=TRUE), ]

# A basic scatterplot with color depending on Species
pdf("fucus_prevalence_indvaldetect_specificity_10.pdf"
    , width = 7 # Default is 7
    , height = 8 # Change to 10; make it taller
)

p <- ggplot(iris, aes(x=specificity, y=frequency, color=core_result)) + 
  scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  geom_point(size=4, alpha=0.8) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))  +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=8, face="bold")) + ylab("Prevalence") + xlab("IndVal") + ggtitle("")
print(p)
dev.off()