# R script for 14th Bioinformatic Lecture (20th July, 2021)

##################### I. Basic study section #########################
## Please select lines and push "Run" to execute following command one by one.

# Basic calculations
1+1		       #addition
4*4		       #multiplication
23/4  	      #division

# Object
a <- 1	      #Save "1" as an object "a" 
b <- 2      	#Save "2" as an object "b" 
a+b            #Add a and b

c <- 1:6  		#Save numerics 1-6 as an object "c"
c	            #Print the object c (same as print(c))

d <- "Hello"    #Save characters "happy" as an object "d"
d 	            # Print the object d

d <- Hello      #Error will be returned if "" is not given.

# Function
sum(a,b)            #Return the sum value of the objects a and b
mean(c)	           	#Return mean value of the object c
sample(c, size=2)	  #Randomly sample numerics from the object c. You can specify the number of sampling
round(3.1415)      	#Round 3.1415

round(mean(c))	#Round the value obtained from mean(c)

saikoro <- function(i){
  die <- 1:6
  dice <- sample(die, size=i, replace=TRUE)
  print(dice)
}

saikoro(i=4)	#Execute the function "saikoro" with the arument i=4

# Package
install.packages("vegan")
library("vegan")


# Programming
x <- 7		#Define x object as a numeric 7

if (x > 5){
  print("large")
} else {
  print("small")
}

# Characters are merged by the function c() to convert to a vector
for (i in c("Loop", "is", "like", "this")){
  print(i)
}

# Clean Environment before going to the analysis
rm(list = ls())

######################################################################
#
#



#################### II. Data analysis section #######################
## Please select lines and push "Run" to execute following command one by one.

############### 1. Environment settings #############################

# Install packages
install.packages("dplyr")     # data manipulation 
install.packages("ggplot2")   # draw figure
#install.packages("vegan")     # ecological statistics
install.packages("pheatmap")  # clustering and heatmap plot
install.packages("tidyverse")
install.packages("tidygraph")
install.packages("ggraph")

# Load packages
library("dplyr")
library("ggplot2")
library("vegan")
library("pheatmap")
library("tidyverse")
library("tidygraph")
library("ggraph")

# Set directory
getwd()   #Print current directory
# If you use Window VID, check if your current directory is [1] "M:/Desktop/11metagenome"
# If "NO", execute the command below.
setwd("M:/Desktop/14metagenome") 

######################################################################
#

############### 2. Load abundance table ##############################

# Load abundance file
raw_df <- read.table("tara_phycodna.txt", header = TRUE)

env_df <- raw_df[,1:4]
phyco_df <- raw_df[,5:994]

# View data frame
# View(raw_df)

######################################################################
#
  
############### 3. Preprocessing for the analyses (remove rare phylotypes) ##################
###
# This section removes phylotypes emerged less than 10 samples from the dataframe
####

# Delete phylotypes emerged <10 samples across all samples
start <- 1
end <- ncol(raw_df)
delete <- NULL  # Declare delete vector
for (i in start:end){
  if (sum(raw_df[,i]!=0) < 10){
    delete <- c(delete, i)
  }
}
reduced_df <- raw_df[,-delete]  #Delete rows recorded as "delete" vector


sum(raw_df[,i] == 0)) ??

# Check if empty OTU exists
if (min(sum(raw_df[,i]!=0)) >= 10){
  print("There is no phylotypes with the coverage less than 10")
} else {
  print("There are OTU(s) which needs to be removed")
}

# Save 
write.table(reduced_df, file="reduced_df.tsv", row.names = TRUE, col.names = TRUE)

######################################################################
#


############### 4. Make Transpose data frame #########################
####
# Then, Transpose data frame for ecological analyses
####

reduced_df <- as.data.frame(t(reduced_df))
env_df <- reduced_df[,1:4]
phyco_df <- reduced_df[,5:297]

# Hereafter, we use "reduced_df" for OTU-based analysis and "reduced_df2" for sample-based analysis

######################################################################
#


############### 4. Rerefaction analysis ##############################

# Show help window
?rarecurve   #Show help window

# Plot rarefaction curve for raw data
rarefaction_raw <- rarecurve(t(raw_df), step = 100, xlab = "Number of reads obtained", ylab = "Number of OTUs observed")

# Plot rarefaction curve for reduced data
rarefaction_reduced <- rarecurve(reduced_df2, step = 100, xlab = "Number of reads obtained", ylab = "Number of OTUs observed")

######################################################################
#


############### 6. Calculate diversity estimates #####################

# Calculate diversity
richness <- specnumber(phyco_df, MARGIN = 1)
shannon <- diversity(phyco_df, MARGIN = 1, index="shannon", base = exp(1)) # base = natural logarithms

# Plot diversity estimates
barplot(richness, xlab = "Sample", ylab = "Number of OTUs")
barplot(shannon, xlab = "Sample", ylab = "Shannon Index")

######################################################################
#


############### 7. NMDS ordination ###################################

# Make distance matrix across samples
df_bray <- vegdist(phyco_df, method="bray") # Bray-Curtis distance
print(df_bray)  #Print distance matrix

#Calculate NMDS coordinates
nm_xy <- metaMDS(df_bray)

# Extract x and y grids
nm_x <- nm_xy$points[,1]
nm_y <- nm_xy$points[,2]

# Make df summarizing the redults
metadf <- data.frame(Name=rownames(reduced_df2),
                     NMDS1 = nm_x, 
                     NMDS2 = nm_y, 
                     Richness = richness,
                     Shannon = shannon,
                     Biome = reduced_df$Biome,
                     Latitude = reduced_df$Latitude)

View(metadf)

# Plot NMDS ordination
nmplot1 <- ggplot(metadf, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(size=3) + 
  #geom_text(aes(x=NMDS1+0.01, label=Name), size=5, hjust=0)
  theme(legend.position = "bottom")
print(nmplot1)

# Plot NMDS ordination with diversity information
nmplot2 <- ggplot(metadf, aes(x=NMDS1, y=NMDS2, color=Richness)) + 
  geom_point(size=3) + 
  #geom_text(aes(x=NMDS1+0.01, label=Name), size=5, hjust=0) + 
  scale_colour_gradient(low="orange", high="blue")
print(nmplot2)

# Plot NMDS ordination with Biome information
nmplot3 <- ggplot(metadf, aes(x=NMDS1, y=NMDS2, fill=Biome)) + 
  geom_point(shape=21, size=3) + 
  #geom_text(aes(x=NMDS1+0.01, label=Name), size=5, hjust=0) + 
  scale_colour_gradient(low="orange", high="blue")
print(nmplot3)

######################################################################
#


############### 8. OTU-cooccurrence analysis ##################################

# Make function to analysze df (it takes time!!)
make_heatmap_clustering <- function(data){
  dist <- vegdist(data, method="bray", na.rm=TRUE) #Bray distance
  dist_matrix <- as.matrix(dist)
  pheatmap(
    dist_matrix,
    main = "",           #Title of the figure
    scale = "none",
    cluster_col= TRUE,   # Draw cluster in column
    cluster_row= TRUE,   # Draw cluster in row
    clustering_method = "average",  #Clusterung method
    fontsize = 8,
    legend=TRUE,
    color = colorRampPalette(c("royalblue", "firebrick3", "white"))(50),
  )
}

# Execute function for OTU data
make_heatmap_clustering(phyco_df)

######################################################################
#


############### 9. Co-occurrence network analysis ##################################

# Convert "charactor" matrix to "numeric" matrix
phyco_df_n <- apply(phyco_df, 2, as.numeric) # 1 or 2で 行列が入れ替わる

# Calcurate correlatiuon between sample
cor.mat <- cor(phyco_df_n, use = "everything", method = "spearman")
# Also try: test <- as.matrix(vegdist(dft, "bray"))

# Matrixの上三角行列を削除
cor.mat[upper.tri(cor.mat, diag = TRUE)] <- NA

# Long-Formatに変換して相関係数0.7以上に絞り込み
d <- cor.mat %>%
  as.data.frame() %>%
  mutate(ids1 = colnames(phyco_df)) %>%
  gather(key = ids2, value = cor, -ids1) %>%
  filter(!is.na(cor) & cor >= 0.7)
head(d)

# percentage of connection against all possible connections
nrow(d)*100/(ncol(phyco_df)^2)


# df -> tbl_graphオブジェクトへの変換
g <- as_tbl_graph(d, directed = FALSE) # if directed graph, chose TRUE

# Add node degree
g <-  g %>% mutate(degree = centrality_degree())


plot1 <- g %>%
  ggraph(layout = "kk") +
  geom_edge_link(aes(width=cor), alpha = 0.6, colour = "gray50") +
  scale_edge_width(range = c(0.1, 1)) +
  geom_node_point(aes(fill=degree), size=3, shape = 21, alpha=0.7) +
  #scale_color_manual(values = colnode_circle) +
  #scale_fill_manual(values = colnodes) +
  scale_fill_gradient(low = "yellow", high = "forestgreen", guide = "colourbar") +
#geom_node_point(aes(size = degree, colour = community, fill=community), shape = 21, alpha=0.7) +
  #geom_node_text(aes(label = name, colour=community,), repel = TRUE) +
  theme_graph(background = "white")
plot1

