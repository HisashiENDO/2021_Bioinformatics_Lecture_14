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
#install.packages("vegan")     # ecological statistics
install.packages("ggplot2")    # draw figure
install.packages("pheatmap")   # clustering and heatmap plot
install.packages("tidyverse")  # manipulate tidy data
install.packages("tidygraph")  # manipulate tidy graph
install.packages("ggraph")     # manipulate ggraph plot
install.packages("Rcpp")       # Integrate R and C++


# Load packages
library("vegan")
library("ggplot2")
library("pheatmap")
library("tidyverse")
library("tidygraph")
library("ggraph")
library("Rcpp")

# Set directory
getwd()   #Print working directory
# If you use Window VDI, check if your current directory is [1] "M:/Desktop/11metagenome"
# If "NO", execute the command below.
setwd("M:/Desktop/14metagenome") 

######################################################################
#

############### 2. Load abundance table ##############################

# Load frequency file
tara_df <- read.table("tara_df.txt", header = TRUE)

# View the profile
#View(tara_df)

######################################################################
#

############### 3. Plot basic data ##############################

# Plot temperature and latitude diagram with generic function
plot(tara_df$Temperature, tara_df$Latitude)

# Plot temperature and latitude diagram with ggplot function
ggplot(tara_df, aes(x = Temperature, y = Latitude, color = Biome)) +
  geom_point()

######################################################################
#


############### 4. Preprocessing for the analyses ##################

## This section removes phylotypes emerged less than 10 samples from the original data

# Separate environmental and frequency information
meta_df <- tara_df[,1:5]
virus_df <- tara_df[,6:995]

# Delete phylotypes emerged <10 samples across all samples
start <- 1
end <- ncol(virus_df)
delete <- NULL  # Declare delete vector
for (i in start:end){
  if (sum(virus_df[,i] != 0) < 10){
    delete <- c(delete, i)
  }
}
virus_df_select <- virus_df[,-delete]  #Delete rows recorded as "delete" vector

# Save 
write.table(virus_df_select, file="virus_df_select.txt", row.names = TRUE, col.names = TRUE)

######################################################################
#



############### 5. Calculate diversity estimates #####################

# Calculate diversity
Richness <- specnumber(virus_df, MARGIN = 1)
Shannon <- diversity(virus_df, MARGIN = 1, index="shannon", base = exp(1)) # base = natural logarithms

# Add diversity estimates to the meta data
meta_df2 <- cbind(meta_df, Richness, Shannon)

# Plot diversity across biomes 
ggplot(meta_df2, aes(x=Biome, y=Richness, fill=Biome)) + geom_boxplot()
ggplot(meta_df2, aes(x=Biome, y=Shannon, fill=Biome)) + geom_boxplot()

######################################################################
#


############### 6. NMDS ordination ###################################

# Make distance matrix across samples
df_bray <- vegdist(virus_df, method="bray") # Bray-Curtis distance
#print(df_bray)  #Print distance matrix

#Calculate NMDS coordinates
nm_xy <- metaMDS(df_bray)

# Extract x and y grids
nm_x <- nm_xy$points[,1]
nm_y <- nm_xy$points[,2]

# Add NMDS coordinates to the meta data
meta_df3 <- data.frame(meta_df2, 
                       NMDS1 = nm_x, 
                       NMDS2 = nm_y
                       )

View(meta_df3)  # View df

# Plot NMDS ordination
ggplot(meta_df3, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(size=2) + 
  geom_text(aes(x=NMDS1+0.01, label=Sample), size=2, hjust=0) +
  theme(legend.position = "bottom")

# Plot NMDS ordination with temperature
ggplot(meta_df3, aes(x=NMDS1, y=NMDS2, fill=Temperature)) + 
  geom_point(shape=21, size=3) + 
  geom_text(aes(x=NMDS1+0.01, label=Sample), size=2, hjust=0) +
  scale_fill_gradient(low="blue", high="red")

# Plot NMDS ordination with diversity information
ggplot(meta_df3, aes(x=NMDS1, y=NMDS2, fill=Richness)) + 
  geom_point(shape=21, size=3) + 
  geom_text(aes(x=NMDS1+0.01, label=Sample), size=2, hjust=0) +
  scale_fill_gradient(low="blue", high="red")

######################################################################
#


############### 7. Heatmap and clustering analysis ##################################

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
    color = colorRampPalette(c("red", "orange", "white"))(20)
  )
}

# Execute function for virus frequency data
make_heatmap_clustering(virus_df)

######################################################################
#


############### 8. Co-abundance network analysis ##################################

# Calcurate pairwise correlations between viruses
cor.mat <- cor(virus_df_select, use = "everything", method = "spearman")

# Remove upper symmetric matrix
cor.mat[upper.tri(cor.mat, diag = TRUE)] <- NA

# View correlation matrix
#View(cor.mat)

# Make dataset used for the network plot
cor.df <- cor.mat %>%
  as.data.frame() %>%         # Convert from matrix to data frame
  mutate(ids1 = rownames(cor.mat)) %>%           # Add column to the data frame
  gather(key = ids2, value = cor, -ids1) %>%     # Transform to long-format 
  filter(!is.na(cor) & cor >= 0.7)               # Remove pairs with the correlation <0.7 
#View(cor.df)

# Convert data frame to tbl_graph object
g <- as_tbl_graph(cor.df, directed = FALSE)

# Add node degree
g <-  g %>% mutate(degree = centrality_degree())

# Plot co-abundance network
netplot <- g %>%
  ggraph(layout = "kk") +
  scale_edge_width(range = c(0.1, 1)) +
  geom_edge_link(aes(width=cor), colour = "gray50") +
  geom_node_point(aes(fill=degree), size=3, shape = 21, alpha=0.7) +
  scale_fill_gradient(low = "yellow", high = "forestgreen", guide = "colourbar") +
  theme_graph(background = "white")
plot(netplot)


# !!! If you got trouble on usable font, please try the following commands.
# This process will be skiped in the lecture as it tales long time..
# library(extrafont)
# font_import()

######################################################################
#


############### 9. Supplementary code  ##################################

# Save plot
fig <- ggplot(tara_df, aes(x = Temperature, y = Latitude, color = Biome)) + 
  geom_point()
ggsave(filename = "fig.pdf", plot= fig, width=120, height=100, units="mm", dpi = 300 )

# Save data object
write.table(virus_df_select, file="virus_df_select.csv", sep = ",", row.names = TRUE, col.names = TRUE)

######################################################################
#
