# ROBUSTNESS ANALYSIS WITH MORE THAN 14 ROUNDS


rm(list=ls(all=TRUE))
myWD <- "/Users/mariapereda/Dropbox/UPM/investigacion/Mis_trabajos_en_curso/PGG_modelito_v2_bayesiano/paraMiPaper/Herding_in_PGG/Code/Robustness"

setwd(myWD)

library("RColorBrewer") 
library('ggplot2')
library('extraDistr')
library('MGLM')
library('reshape2')
library('XNomial')
library('dtw')
library('gridExtra')
library('matrixcalc')

set.seed(5)


########### MODEL A: Bayesian agents with non-informative prior ############################################
run_model1 <- function(myreplications=2, N=100, r=14) {
  #N Number of agents
  
  print(myreplications)
  
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- r             # Number of rounds
  vector_alfa_inicial <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5) #Dirichlet Jeffreys Prior
  
  # Structure to store the distributions of decisions per round
  distributions <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(distributions)[1]<-"round"
  colnames(distributions)[2]<-"contribution"
  colnames(distributions)[3]<-"frequency"
  colnames(distributions)[4]<-"mean"
  colnames(distributions)[5]<-"replicationNumber"
  indice<-1 #just to store the data
  medias_por_ronda <- c(0,0,0,0,0,0)
  
  while(myreplications>0){
    vector_alfa= vector_alfa_inicial
    for (t in 1:final_round){
      # Agents decide
      x<-rcat(N, vector_alfa)
      
      #Posterior distribution update
      vector_c <- c(0,0,0,0,0,0)
      for (c_indice in 1:length(vector_alfa_inicial)) {
        categorias_datos <- as.numeric(names(table(x)))
        cuenta_datos <- as.numeric(table(x))
        vector_c[c_indice] <- if(any(categorias_datos == c_indice)){cuenta_datos[which(categorias_datos == c_indice)]}else{0}
      }
      alfas_posterior <- vector_alfa + vector_c
      
      # Saving distributions of decisions per round
      for (categoria in 0:5){
        distributions[indice+categoria,]$round <- t
        distributions[indice+categoria,]$contribution <- decisiones[1+categoria]
        distributions[indice+categoria,]$frequency <- vector_c[1+categoria]/N
        distributions[indice+categoria,]$mean <- sum(vector_c*decisiones)/N
        distributions[indice+categoria,]$replicationNumber <- myreplications
      }
      indice = indice +6
      
      #Posterior predictive distribution update
      pesos_post_predictive <- c(0,0,0,0,0,0)
      for (i in 1:length(pesos_post_predictive)) {
        pesos_post_predictive[i]<- (alfas_posterior[i])/(sum(vector_alfa)+length(x))
      }
      
      
      #Alfas updtate according to the Posterior predictive distribution
      vector_alfa <- pesos_post_predictive
    }
    myreplications<-myreplications-1
    print(myreplications)
  }
  
  return(distributions)
}


############# N=100 ###########
# I generate datasets of each model with n replications to cluster behavioural types.
replications<-100
dataModel1 <- run_model1(replications,100,30)
dataModel1_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel1_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel1_ <- rbind(dataModel1_,dataModel1[dataModel1$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel1_ <- rbind(dataModel1_,dataModel1[dataModel1$replicationNumber==d,]$frequency)
}
colnames(dataModel1_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


#I use data with name ending in _ to cluster

dist_mat <- dist(dataModel1_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

# I look at the groups in the dedrogram to obtain k clusters

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 6)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 6)
table(clusters)


clusterdata <- dataModel1[dataModel1$replicationNumber %in% which(clusters==6),]

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(32,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]


lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("Examples_Model1_100_cluster6_hierarchicalClust_average.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

#### Summary of 15 replications for the Supplementary Material

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model1_100_cluster6_hierarchicalClust_average_30rounds.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## I print a replication as an example
unique(clusterdata$replicationNumber)
rep=12

setEPS()
postscript("Model1_100_cluster2_rep12.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel1[dataModel1$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel1[dataModel1$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel1[dataModel1$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel1[dataModel1$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


# Scatterplot all
setEPS()
postscript("Scatter_Bayesianos100_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel1, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()

############# N=1000 ###########
replications<-100
dataModel1_1000 <- run_model1(replications,1000,30)
dataModel1_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel1_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel1_1000_ <- rbind(dataModel1_1000_,dataModel1_1000[dataModel1_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel1_1000_ <- rbind(dataModel1_1000_,dataModel1_1000[dataModel1_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel1_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dist_mat <- dist(dataModel1_1000_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#I look at the groups in the dedrogram to obtain k clusters.

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 2)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 2)
table(clusters)


# I generate examples of each cluster
clusterdata <- dataModel1_1000[dataModel1_1000$replicationNumber %in% which(clusters==2),]

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(32,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("Examples_Model1_1000_cluster1_hierarchicalClust_average.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()


#### Summary of 15 replications for the Supplementary Material

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model1_1000_cluster2_hierarchicalClust_average.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## I print a replication as an example
rep=38

setEPS()
postscript("Model1_1000_cluster1_rep38.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


# Scatterplot all
setEPS()
postscript("Scatter_Bayesianos100_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel1, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()



################################## ################################## ################################## ################################## 

set.seed(17)

########### MODEL B: Bayesian agents with data-based informative prior ############################################

data<-read.csv("/Users/mariapereda/Dropbox/UPM/investigacion/Mis_trabajos_en_curso/PGG_modelito_v2_bayesiano/paraMiPaper/datos/2017A2ECOCOOPANSA01ONL0000ESPMAD.csv") #DATA FROM ZENODO https://zenodo.org/record/2590686
#Create participant IDs combining participant_ID and treatment
data$id <- paste(data$treatment,data$participant_id)
dataset <- data[data$player_alive==1, ]
#select active participants and first round decisions:
subset_data<-subset(data, data$player_alive==1 & data$round_number==1)
participants_number<-nrow(subset_data)
dataprior<- tapply(subset_data$round_number, subset_data$player_contribution, sum)/participants_number #average of all experiments


run_model2 <- function(myreplications=1, N=100, r=14) {
  #N Number of agents
  
  ### PGG 14 rondas
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- r             # Number of rounds
  vector_alfa_inicial <- as.numeric(dataprior) 
  
  # Structure to store the distributions of decisions per round
  distributions <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(distributions)[1]<-"round"
  colnames(distributions)[2]<-"contribution"
  colnames(distributions)[3]<-"frequency"
  colnames(distributions)[4]<-"mean"
  colnames(distributions)[5]<-"replicationNumber"
  indice<-1 #just to store the data
  medias_por_ronda <- c(0,0,0,0,0,0)
  
  while(myreplications>0){
    vector_alfa= vector_alfa_inicial
    for (t in 1:final_round){
      # Agents decide
      x<-rcat(N, vector_alfa)
      
      #Posterior distribution update
      vector_c <- c(0,0,0,0,0,0)
      for (c_indice in 1:length(vector_alfa_inicial)) {
        categorias_datos <- as.numeric(names(table(x)))
        cuenta_datos <- as.numeric(table(x))
        vector_c[c_indice] <- if(any(categorias_datos == c_indice)){cuenta_datos[which(categorias_datos == c_indice)]}else{0}
      }
      alfas_posterior <- vector_alfa + vector_c
      
      # Saving distributions of decisions per round
      for (categoria in 0:5){
        distributions[indice+categoria,]$round <- t
        distributions[indice+categoria,]$contribution <- decisiones[1+categoria]
        distributions[indice+categoria,]$frequency <- vector_c[1+categoria]/N
        distributions[indice+categoria,]$mean <- sum(vector_c*decisiones)/N
        distributions[indice+categoria,]$replicationNumber <- myreplications
      }
      indice = indice +6
      
      #Posterior predictive distribution update
      pesos_post_predictive <- c(0,0,0,0,0,0)
      for (i in 1:length(pesos_post_predictive)) {
        pesos_post_predictive[i]<- (alfas_posterior[i])/(sum(vector_alfa)+length(x))
      }
      
      
      #Alfas updtate according to the Posterior predictive distribution
      vector_alfa <- pesos_post_predictive
    }
    myreplications<-myreplications-1
    print(myreplications)
  }
  return(distributions)
}

############# N=100 ###########
replications<-100
dataModel2 <- run_model2(replications,100,30)
dataModel2_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel2_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel2_ <- rbind(dataModel2_,dataModel2[dataModel2$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel2_ <- rbind(dataModel2_,dataModel2[dataModel2$replicationNumber==d,]$frequency)
}
colnames(dataModel2_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


############# N=1000 ###########
replications<-100
dataModel2_1000 <- run_model2(replications,1000,30)
dataModel2_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel2_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel2_1000_ <- rbind(dataModel2_1000_,dataModel2_1000[dataModel2_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel2_1000_ <- rbind(dataModel2_1000_,dataModel2_1000[dataModel2_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel2_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


#Scatterplots
setEPS()
postscript("Scatter_Bayesianos_DataPrior_100_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel2, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()

setEPS()
postscript("Scatter_Bayesianos_DataPrior_1000_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel2_1000, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


## CLUSTER ANALYSIS
############# N=100 ###########

#I use data with name ending in _ to cluster
dist_mat <- dist(dataModel2_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#I look at the groups in the dedrogram to obtain k clusters.

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 6)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 6)
table(clusters)


clusterdata <- dataModel2[dataModel2$replicationNumber %in% which(clusters==6),]

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(32,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("Examples_Model2_100_cluster1_hierarchicalClust_average_30rounds.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

#### Summary of 15 replications for the Supplementary Material

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model2_100_cluster6_hierarchicalClust_average_30rounds.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## I print a replication as an example
rep=sample_indexes[10]

setEPS()
postscript("Model2_100_cluster2_rep77.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel2[dataModel2$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel2[dataModel2$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel2[dataModel2$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel2[dataModel2$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


############# N=1000 ###########

#I use data with name ending in _ to cluster
dist_mat <- dist(dataModel2_1000_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#I look at the groups in the dedrogram to obtain k clusters.

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 2)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 2)
table(clusters)


clusterdata <- dataModel2_1000[dataModel2_1000$replicationNumber %in% which(clusters==1),]

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(32,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("Examples_Model2_1000_cluster1_hierarchicalClust_average.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

#### Summary of 15 replications for the Supplementary Material

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model2_1000_cluster1_hierarchicalClust_average.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## I print a replication as an example
rep=sample_indexes[1]

setEPS()
postscript("Model2_1000_cluster2_rep36.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()



################################## ################################## ################################## ################################## 

set.seed(19)

########### MODEL C: Bayesian and non-bayesian agents with data-based informative prior ############################################

# Lets define free rider: person who donates 0 in at least a 80% times: 11 rounds
# Lets define full cooperator: person who donates 10 in at least a 80% times: 11 rounds


fulls <- dataset[dataset$player_contribution==0,] #full defectors
percentage_full_def <- sum(as.numeric(table(fulls$id))>11)/participants_number #Donate 0 more than 80% rounds

coops <- dataset[dataset$player_contribution==10,] #full cooperators
percentage_full_cop<-sum(as.numeric(table(coops$id))>11)/participants_number #Donate 10 more than 80% rounds


run_model3 <- function(myreplications=1, N=100, r=14) {
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- r             # Number of rounds
  #N <- 100                      # Number of agents
  n_prime=N-round(percentage_full_def,2)*100-round(percentage_full_cop,2)*100 # Bayesian agents
  round(percentage_full_def*N,2)*100 #Full defector
  round(percentage_full_cop*N,2)*100 #Full cooperators
  vector_alfa_inicial <- as.numeric(dataprior) 
  
  # Structure to store the distributions of decisions per round
  distributions <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(distributions)[1]<-"round"
  colnames(distributions)[2]<-"contribution"
  colnames(distributions)[3]<-"frequency"
  colnames(distributions)[4]<-"mean"
  colnames(distributions)[5]<-"replicationNumber"
  indice<-1 #just to store the data
  medias_por_ronda <- c(0,0,0,0,0,0)
  
  while(myreplications>0){
    vector_alfa= vector_alfa_inicial
    for (t in 1:final_round){
      # Agents decide
      x<-rcat(n_prime, vector_alfa)
      x<-c(x,rep(1,round(percentage_full_def,2)*100), rep(6, round(percentage_full_cop,2)*100))  # Category 1 y donate 0, category 6 is donate 10
      
      #Posterior distribution update
      vector_c <- c(0,0,0,0,0,0)
      for (c_indice in 1:length(vector_alfa_inicial)) {
        categorias_datos <- as.numeric(names(table(x)))
        cuenta_datos <- as.numeric(table(x))
        vector_c[c_indice] <- if(any(categorias_datos == c_indice)){cuenta_datos[which(categorias_datos == c_indice)]}else{0}
      }
      alfas_posterior <- vector_alfa + vector_c
      
      # Saving distributions of decisions per round
      for (categoria in 0:5){
        distributions[indice+categoria,]$round <- t
        distributions[indice+categoria,]$contribution <- decisiones[1+categoria]
        distributions[indice+categoria,]$frequency <- vector_c[1+categoria]/N
        distributions[indice+categoria,]$mean <- sum(vector_c*decisiones)/N
        distributions[indice+categoria,]$replicationNumber <- myreplications
      }
      indice = indice +6
      
      #Posterior predictive distribution update
      pesos_post_predictive <- c(0,0,0,0,0,0)
      for (i in 1:length(pesos_post_predictive)) {
        pesos_post_predictive[i]<- (alfas_posterior[i])/(sum(vector_alfa)+length(x))
      }
      
      
      #Alfas updtate according to the Posterior predictive distribution
      vector_alfa <- pesos_post_predictive
    }
    myreplications<-myreplications-1
    print(myreplications)
  }
  return(distributions)
}

############# N=100 ###########
replications<-100
dataModel3 <- run_model3(replications,100,30)
dataModel3_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel3_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel3_ <- rbind(dataModel3_,dataModel3[dataModel3$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel3_ <- rbind(dataModel3_,dataModel3[dataModel3$replicationNumber==d,]$frequency)
}
colnames(dataModel3_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

############# N=1000 ###########
replications<-100
dataModel3_1000 <- run_model3(replications,1000,30)
dataModel3_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel3_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel3_1000_ <- rbind(dataModel3_1000_,dataModel3_1000[dataModel3_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel3_1000_ <- rbind(dataModel3_1000_,dataModel3_1000[dataModel3_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel3_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


#Scatterplots
setEPS()
postscript("Scatter_Bayesianos_DataPrior_Fulls_100_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()

setEPS()
postscript("Scatter_Bayesianos_DataPrior_Fulls_1000_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3_1000, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


## CLUSTER ANALYSIS
############# N=100 ###########

#I use data with name ending in _ to cluster
dist_mat <- dist(dataModel3_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#I look at the groups in the dedrogram to obtain k clusters.

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 4)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 4)
table(clusters)


clusterdata <- dataModel3[dataModel3$replicationNumber %in% which(clusters==4),]

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(32,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("Examples_Model3_100_cluster1_hierarchicalClust_average.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

#### Summary of 15 replications for the Supplementary Material

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model3_100_cluster4_hierarchicalClust_average_30rounds.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


### Plots for 100 rounds
clusterdata <- dataModel3[dataModel3$replicationNumber %in% which(clusters==4),]
sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 100), breaks=seq(1,100,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=10),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model3_100_cluster4_hierarchicalClust_average_100rounds.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()

################################





## I print a replication as an example
rep=sample_indexes[1]

setEPS()
postscript("Model3_100_cluster2_rep90.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3[dataModel3$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel3[dataModel3$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel3[dataModel3$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel3[dataModel3$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


############# N=1000 ###########

#I use data with name ending in _ to cluster
dist_mat <- dist(dataModel3_1000_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#I look at the groups in the dedrogram to obtain k clusters.

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 4)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 4)
table(clusters)


clusterdata <- dataModel3_1000[dataModel3_1000$replicationNumber %in% which(clusters==4),]

#### Summary of 15 replications for the Supplementary Material

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model3_1000_cluster3_hierarchicalClust_average.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## I print a replication as an example
rep=sample_indexes[1]

setEPS()
postscript("Model3_1000_cluster4_rep43.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()





################################## ################################## ################################## ################################## 

set.seed(21)

########### MODEL D: Bayesian and non-bayesian agents with data-based noninformative prior ############################################

run_model4 <- function(myreplications=1, N=100,r=14) {
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- r             # Number of rounds
  #N <- 100                     # Number of agents
  n_prime=N-round(percentage_full_def,2)*100-round(percentage_full_cop,2)*100 # Bayesian agents
  round(percentage_full_def*N,2)*100 #Full defector
  round(percentage_full_cop*N,2)*100 #Full cooperators
  vector_alfa_inicial <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5) #Dirichlet Jeffreys Prior
  
  # Structure to store the distributions of decisions per round
  distributions <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(distributions)[1]<-"round"
  colnames(distributions)[2]<-"contribution"
  colnames(distributions)[3]<-"frequency"
  colnames(distributions)[4]<-"mean"
  colnames(distributions)[5]<-"replicationNumber"
  indice<-1 #just to store the data
  medias_por_ronda <- c(0,0,0,0,0,0)
  
  while(myreplications>0){
    vector_alfa= vector_alfa_inicial
    for (t in 1:final_round){
      # Agents decide
      x<-rcat(n_prime, vector_alfa)
      x<-c(x,rep(1,round(percentage_full_def,2)*100), rep(6, round(percentage_full_cop,2)*100))  # Category 1 y donate 0, category 6 is donate 10
      
      #Posterior distribution update
      vector_c <- c(0,0,0,0,0,0)
      for (c_indice in 1:length(vector_alfa_inicial)) {
        categorias_datos <- as.numeric(names(table(x)))
        cuenta_datos <- as.numeric(table(x))
        vector_c[c_indice] <- if(any(categorias_datos == c_indice)){cuenta_datos[which(categorias_datos == c_indice)]}else{0}
      }
      alfas_posterior <- vector_alfa + vector_c
      
      # Saving distributions of decisions per round
      for (categoria in 0:5){
        distributions[indice+categoria,]$round <- t
        distributions[indice+categoria,]$contribution <- decisiones[1+categoria]
        distributions[indice+categoria,]$frequency <- vector_c[1+categoria]/N
        distributions[indice+categoria,]$mean <- sum(vector_c*decisiones)/N
        distributions[indice+categoria,]$replicationNumber <- myreplications
      }
      indice = indice +6
      
      #Posterior predictive distribution update
      pesos_post_predictive <- c(0,0,0,0,0,0)
      for (i in 1:length(pesos_post_predictive)) {
        pesos_post_predictive[i]<- (alfas_posterior[i])/(sum(vector_alfa)+length(x))
      }
      
      
      #Alfas updtate according to the Posterior predictive distribution
      vector_alfa <- pesos_post_predictive
    }
    myreplications<-myreplications-1
    print(myreplications)
  }
  return(distributions)
}

############# N=100 ###########
replications<-100
dataModel4 <- run_model4(replications,100,30)
dataModel4_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel4_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel4_ <- rbind(dataModel4_,dataModel4[dataModel4$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel4_ <- rbind(dataModel4_,dataModel4[dataModel4$replicationNumber==d,]$frequency)
}
colnames(dataModel4_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

############# N=1000 ###########
replications<-100
dataModel4_1000 <- run_model4(replications,1000,30)
dataModel4_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel4_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel4_1000_ <- rbind(dataModel4_1000_,dataModel4_1000[dataModel4_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel4_1000_ <- rbind(dataModel4_1000_,dataModel4_1000[dataModel4_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel4_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


#Scatterplots
setEPS()
postscript("Scatter_Bayesianos_JeffreysPrior_Fulls_100_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()

setEPS()
postscript("Scatter_Bayesianos_JeffreysPrior_Fulls_1000_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4_1000, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


## CLUSTER ANALYSIS
############# N=100 ###########

#I use data with name ending in _ to cluster
dist_mat <- dist(dataModel4_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#I look at the groups in the dedrogram to obtain k clusters.

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 6)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 6)
table(clusters)


clusterdata <- dataModel4[dataModel4$replicationNumber %in% which(clusters==1),]

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(32,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("Examples_Model4_100_cluster2_hierarchicalClust_average.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

#### Summary of 15 replications for the Supplementary Material

clusterdata <- dataModel4[dataModel4$replicationNumber %in% which(clusters==4),]
sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model4_100_cluster4_hierarchicalClust_average_30rounds.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## I print a replication as an example
rep=sample_indexes[2]

setEPS()
postscript("Model4_100_cluster1_rep61.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4[dataModel4$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel4[dataModel4$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel4[dataModel4$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel4[dataModel4$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


############# N=1000 ###########

#I use data with name ending in _ to cluster
dist_mat <- dist(dataModel4_1000_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#I look at the groups in the dedrogram to obtain k clusters.

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 2)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 2)
table(clusters)


clusterdata <- dataModel4_1000[dataModel4_1000$replicationNumber %in% which(clusters==2),]

#### Summary of 15 replications for the Supplementary Material

sample_indexes <- sample(unique(clusterdata$replicationNumber), min(15,length(unique(clusterdata$replicationNumber))))
sample <- clusterdata[clusterdata$replicationNumber %in% sample_indexes, ]

lista_plots <- list()
replication_number<-unique(sample$replicationNumber)
indice<-1
for (i in replication_number){
  lista_plots[[indice]]<-ggplot(sample[sample$replicationNumber==i,], aes(x=round, y=contribution)) + 
    geom_point(size=20, shape=15, aes(colour = frequency)) + 
    #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
    scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
    scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
    scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model4_1000_cluster2_hierarchicalClust_average.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## I print a replication as an example
rep=sample_indexes[1]

setEPS()
postscript("Model4_1000_cluster2.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 30), breaks=seq(1,30,1)) +
  geom_line(data=dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()



#### Export data for the repositoy
write.csv(dataModel1, "dataModel1.csv", row.names=TRUE)
write.csv(dataModel1_1000, "dataModel1_1000.csv", row.names=TRUE)
write.csv(dataModel2, "dataModel2.csv", row.names=TRUE)
write.csv(dataModel2_1000, "dataModel2_1000.csv", row.names=TRUE)
write.csv(dataModel3, "dataModel3.csv", row.names=TRUE)
write.csv(dataModel3_1000, "dataModel3_1000.csv", row.names=TRUE)
write.csv(dataModel4, "dataModel4.csv", row.names=TRUE)
write.csv(dataModel4_1000, "dataModel4_1000.csv", row.names=TRUE)

