setwd('Desktop/SRP/')
#####################
#############functions
AddCommon <- function(data){
  library('org.Sc.sgd.db')
  xx= as.list(org.Sc.sgdCOMMON2ORF)
  yy = sapply(xx, function(x){x[1]})
  common=names(xx)
  orf = as.vector(yy)
  names(common)=orf
  
  #add systematic names to dataset
  data$gene <- as.vector(common[as.vector(data$NAME)])
  genemissing_positions <- is.na(data$gene)
  data$gene[genemissing_positions] = as.vector(data$NAME[genemissing_positions])
  rmlist <- data$gene[duplicated(data$gene)] #get dups
  data <-data[!(data$gene %in% rmlist),] #remove dups
  return(data)
}


#####################
############load data
shapes <- read.delim(file="yoshi_shape_description.txt",header=TRUE,sep="\t")
df <- read.delim(file="yoshi_sig_data.txt",header=TRUE,sep="\t")
df <- AddCommon(df)
#####################
#####################

#####################
##########load libraries
library(tidyverse)
#####################
#####################

#####################
###########process

#get SRP genes

MyGenes <- c("SRP54","SRP14","SRP72","SRP101","SEC65")
srp_positions <- df$gene %in% MyGenes
df_srp <- df[srp_positions,]


#convert to tidy
df_srp <- df_srp[,-1] #get rid of essential genes
srp_tidy <- df_srp %>% 
  gather(feature,value,-gene) %>% 
  filter(value <0.01)


#do a join with descriptions
names(shapes) = c("feature","description")

final <- left_join(srp_tidy,shapes,by="feature")



