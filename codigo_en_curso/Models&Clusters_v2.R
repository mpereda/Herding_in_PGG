# Models and clusters





# ver sección "1 de noviembre me quedo aquí"



################################## Noviembre 2023 #############

## Actualizo la técnica de clustering a clustering jerárquico aglomerativo con average linkage (hay que hacer sufifientes grupos)

rm(list=ls(all=TRUE))
myWD <- "/Users/mariapereda/Dropbox/UPM/investigacion/Mis_trabajos_en_curso/PGG_modelito_v2_bayesiano/paraMiPaper/histogramsOfHistograms"

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
run_model1 <- function(myreplications=2, N=100) {
  #N Number of agents
  
  print(myreplications)
  
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- 14             # Number of rounds
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
# Genero datasets de cada modelo con n replicaciones para hacer clustering de tipos de comportamientos
replications<-100
dataModel1 <- run_model1(replications,100)
dataModel1_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel1_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel1_ <- rbind(dataModel1_,dataModel1[dataModel1$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel1_ <- rbind(dataModel1_,dataModel1[dataModel1$replicationNumber==d,]$frequency)
}
colnames(dataModel1_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


#Utilizo los datos con nombre acabado en _ para clusterizar

dist_mat <- dist(dataModel1_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#observo los grupos en el dedrograma para obtener k

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 6)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 6)
table(clusters)


clusterdata <- dataModel1[dataModel1$replicationNumber %in% which(clusters==3),]

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
    scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("Examples_Model1_100_cluster3_hierarchicalClust_average.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

#### Resumen de 15 para el SI

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
    scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}
setEPS()
postscript("15Examples_Model1_100_cluster3_hierarchicalClust_average.eps", height=4*5, width=8.5*3, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,15), 5, 3, byrow=TRUE))
dev.off()


## Replicación que representa cada cluster
rep=12

setEPS()
postscript("Model1_100_cluster1_rep12.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel1[dataModel1$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel1[dataModel1$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel1[dataModel1$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel1[dataModel1$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
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
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()

############# N=1000 ###########
replications<-100
dataModel1_1000 <- run_model1(replications,1000)
dataModel1_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel1_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel1_1000_ <- rbind(dataModel1_1000_,dataModel1_1000[dataModel1_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel1_1000_ <- rbind(dataModel1_1000_,dataModel1_1000[dataModel1_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel1_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')





################# 1 de noviembre me quedo aquí

# Para cada modelo hay que hacer análisis clust de hierarchical, con bastantes grupos para que se agrupen bien los patrones
# Veo los patrones en una figura con todos, pero al supplementary subo una de tamaño 15 ejemplos o menos (las de menos, edito el pdf)
# Comento en el texto y para las representativas, muestro un ejemplo como figura sola
# Fisher test solo contra el cluster que se parezca, contra PGG_H, PGG_H2, y PGG100

#repito todo con 1000 agentes



###################################################



















distMatrix <- dist(dataModel1_1000_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel1_1000_)
plot(hc, labels=observedLabels, main="")

#num clusters selected
clustnum <- 2
clusterCut <- cutree(hc, clustnum)
table(clusterCut) #Percentaje of replications on each cluster.  99  1 

rep=1

setEPS()
postscript("Model1_1000_cluster1_99perc.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel1_1000[dataModel1_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


# Scatterplot all
setEPS()
postscript("Scatter_Bayesianos1000_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel1_1000, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


## Fisher tests

#Using data from real humans from the paper
counts_todas_rondas<-read.csv("/Users/mariapereda/Dropbox/UPM/investigacion/Mis_trabajos_en_curso/PGG_modelito_v2_bayesiano/paraMiPaper/datos/counts_todas_rondas.csv") #DATA FROM ZENODO https://zenodo.org/record/2590686
dist_agregada_PGG_H<-dcast(counts_todas_rondas[counts_todas_rondas$treatment=="PGG_H",], round ~ contribution, value.var = "COUNT", fun.aggregate=sum)
colnames(dist_agregada_PGG_H)<-c("round","c0","c2","c4","c6","c8","c10")


# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:final_round){
  for (r in 1:min(10,max(dataModel1$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel1[dataModel1$round==t & dataModel1$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # p-value>alpha means the data fit the model
    # We need at least one data in each category, so I approximate the p-value ensuring all categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001) #simulation data does not fit experimental results
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel1$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)

#PPGH_2
dist_agregada_PGG_H2<-dcast(counts_todas_rondas[counts_todas_rondas$treatment=="PGG_H2",], round ~ contribution, value.var = "COUNT", fun.aggregate=sum)
colnames(dist_agregada_PGG_H2)<-c("round","c0","c2","c4","c6","c8","c10")


# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:final_round){
  for (r in 1:min(10,max(dataModel1$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel1[dataModel1$round==t & dataModel1$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H2[dist_agregada_PGG_H2$round==t,-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel1$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)



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


run_model2 <- function(myreplications=1, N=100) {
  #N Number of agents
  
  ### PGG 14 rondas
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- 14             # Number of rounds
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

replications<-100
dataModel2 <- run_model2(replications,100)
dataModel2_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel2_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel2_ <- rbind(dataModel2_,dataModel2[dataModel2$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel2_ <- rbind(dataModel2_,dataModel2[dataModel2$replicationNumber==d,]$frequency)
}
colnames(dataModel2_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


library(dtw)
distMatrix <- dist(dataModel2_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel2_)
plot(hc, labels=observedLabels, main="")

#num clusters selected
clustnum <- 3
clusterCut <- cutree(hc, clustnum)
table(clusterCut) #Percentaje of replications on each cluster.  94  3  3 

#selecciono la replicacion a imprimir
#1, 8, 38 

rep=38

setEPS()
postscript("Model2_cluster3_4perc.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel2[dataModel2$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel2[dataModel2$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel2[dataModel2$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel2[dataModel2$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()





############# N=1000 ###########
replications<-100
dataModel2_1000 <- run_model2(replications,1000)
dataModel2_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel2_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel2_1000_ <- rbind(dataModel2_1000_,dataModel2_1000[dataModel2_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel2_1000_ <- rbind(dataModel2_1000_,dataModel2_1000[dataModel2_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel2_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


distMatrix <- dist(dataModel2_1000_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel2_1000_)
plot(hc, labels=observedLabels, main="")

#num clusters selected
clustnum <- 2
clusterCut <- cutree(hc, clustnum)
table(clusterCut) #Percentaje of replications on each cluster.  99  1 

rep=1

setEPS()
postscript("Model2_1000_cluster1_98perc.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel2_1000[dataModel2_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()




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
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
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
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()

# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round-1)){
  for (r in 1:min(10,max(dataModel2$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel2[dataModel2$round==t & dataModel2$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel2$replicationNumber))
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.94 SD0.22

#PPGH_2
# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round-1)){
  for (r in 1:min(10,max(dataModel2$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel2[dataModel2$round==t & dataModel2$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H2[dist_agregada_PGG_H2$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel2$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.90 SD0.25



### Comparison with PGG 100
dist_agregada_PGG100<-dcast(counts_todas_rondas[counts_todas_rondas$treatment=="PGG100",], round ~ contribution, value.var = "COUNT", fun.aggregate=sum)
colnames(dist_agregada_PGG100)<-c("round","c0","c2","c4","c6","c8","c10")


# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round-1)){
  for (r in 1:min(10,max(dataModel2$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel2[dataModel2$round==t & dataModel2$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG100[dist_agregada_PGG100$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel2$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.92 SD0.07


################################## ################################## ################################## ################################## 

set.seed(19)

########### MODEL C: Bayesian and non-bayesian agents with data-based informative prior ############################################

# Lets define free rider: person who donates 0 in at least a 80% times: 11 rounds
# Lets define full cooperator: person who donates 10 in at least a 80% times: 11 rounds


fulls <- dataset[dataset$player_contribution==0,] #full defectors
percentage_full_def <- sum(as.numeric(table(fulls$id))>11)/participants_number #Donate 0 more than 80% rounds

coops <- dataset[dataset$player_contribution==10,] #full cooperators
percentage_full_cop<-sum(as.numeric(table(coops$id))>11)/participants_number #Donate 10 more than 80% rounds


run_model3 <- function(myreplications=1, N=100) {
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- 14             # Number of rounds
  #N <- 100                      # Number of agents
  n_prime=N-round(percentage_full_def,0)-round(percentage_full_cop,0) # Bayesian agents
  round(percentage_full_def*N,0) #Full defector
  round(percentage_full_cop*N,0) #Full cooperators
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
      x<-c(x,rep(1,round(percentage_full_def,0)), rep(6, round(percentage_full_cop,0)))  # Category 1 y donate 0, category 6 is donate 10
      
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

replications<-100
dataModel3 <- run_model2(replications,100)
dataModel3_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel3_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel3_ <- rbind(dataModel3_,dataModel3[dataModel3$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel3_ <- rbind(dataModel3_,dataModel3[dataModel3$replicationNumber==d,]$frequency)
}
colnames(dataModel3_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

#Cluster analysis
library(dtw)
distMatrix <- dist(dataModel3_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel3_)
plot(hc, labels=observedLabels, main="")

#num clusters selected
clustnum <- 3
clusterCut <- cutree(hc, clustnum)
table(clusterCut) #Percentaje of replications on each cluster.  94  3  3 

#selecciono la replicacion a imprimir
#24, 7, 29
which(clusterCut==1)
# 1   2   3   4   5   6   8   9  10  11  12  16  17  18  19  20  21  22
#23  ·24  25  26  27  30  31  33  34  37  38  39  40  43  44  45  47  49
#52  54  55  56  57  59  60  61  62  63  66  68  69  ·70  71  72  ·73  75
#76  78  80  82  83  85  87  89  90  92  94  95  96  97  98  99 100
which(clusterCut==2)
# ·7 14 15 32 35 41 42 50 51 67 74 79 81 84 93
which(clusterCut==3)
# 13 28 ·29 36 46 48 53 58 64 65 77 86 88 91

rep=29

setEPS()
postscript("Model3_cluster3_14perc.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3[dataModel3$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel3[dataModel3$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel3[dataModel3$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel3[dataModel3$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()



############# N=1000 ###########
replications<-100
dataModel3_1000 <- run_model3(replications,1000)
dataModel3_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel3_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel3_1000_ <- rbind(dataModel3_1000_,dataModel3_1000[dataModel3_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel3_1000_ <- rbind(dataModel3_1000_,dataModel3_1000[dataModel3_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel3_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


distMatrix <- dist(dataModel3_1000_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel3_1000_)
plot(hc, labels=observedLabels, main="")

#num clusters selected
clustnum <- 4
clusterCut <- cutree(hc, clustnum)
table(clusterCut) #Percentaje of replications on each cluster.  99  1 

which(clusterCut==1) #75
which(clusterCut==2) #6

rep=6

setEPS()
postscript("Model3_1000_cluster2_32perc.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel3_1000[dataModel3_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


#Scatterplots
setEPS()
postscript("Scatter_Bayesianos_DataPrior_100_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()

setEPS()
postscript("Scatter_Bayesianos_DataPrior_1000_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel3_1000, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  #geom_line(data=distributions, aes(x=round, y=mean))+
  #geom_point(data=distributions, aes(x=round, y=mean))+
  #geom_text(data=distributions, aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round-1)){
  for (r in 1:min(10,max(dataModel3$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel3[dataModel3$round==t & dataModel3$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel3$replicationNumber))
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.93 SD0.22

#PPGH_2
# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round-1)){
  for (r in 1:min(10,max(dataModel3$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel3[dataModel3$round==t & dataModel3$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H2[dist_agregada_PGG_H2$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel3$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.88 SD0.24



### Comparison with PGG 100
dist_agregada_PGG100<-dcast(counts_todas_rondas[counts_todas_rondas$treatment=="PGG100",], round ~ contribution, value.var = "COUNT", fun.aggregate=sum)
colnames(dist_agregada_PGG100)<-c("round","c0","c2","c4","c6","c8","c10")


# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round-1)){
  for (r in 1:min(10,max(dataModel3$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel3[dataModel3$round==t & dataModel3$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG100[dist_agregada_PGG100$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel3$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.82 SD0.14

#####################################################
# Fisher tests for specific clusters data
distMatrix <- dist(dataModel3_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel3_)
clustnum <- 3
clusterCut <- cutree(hc, clustnum)

### We compare Cluster 3 with PGG_H data
clusterdata <- dataModel3[dataModel3$replicationNumber %in% which(clusterCut==3),]
# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14
for (t in 1:(final_round-1)){
  for (r in unique(clusterdata$replicationNumber)){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(clusterdata[clusterdata$round==t & clusterdata$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}
summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(length(unique(clusterdata$replicationNumber)),max(dataModel3$replicationNumber))
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.97 SD0.1


### We compare Cluster 1 with PGG100 data
clusterdata <- dataModel3[dataModel3$replicationNumber %in% which(clusterCut==1),]
# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14
for (t in 1:(final_round-1)){
  for (r in unique(clusterdata$replicationNumber)){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(clusterdata[clusterdata$round==t & clusterdata$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG100[dist_agregada_PGG100$round==(t+1),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}
summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(length(unique(clusterdata$replicationNumber)),max(dataModel3$replicationNumber))
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.82 SD0.11



################################## ################################## ################################## ################################## 

set.seed(21)

########### MODEL D: Bayesian and non-bayesian agents with data-based noninformative prior ############################################

run_model4 <- function(myreplications=1, N=100) {
  decisiones <- c(0,2,4,6,8,10) # For figure axis
  final_round <- 14             # Number of rounds
  #N <- 100                     # Number of agents
  n_prime=N-round(percentage_full_def,0)-round(percentage_full_cop,0) # Bayesian agents
  round(percentage_full_def*N,0) #Full defector
  round(percentage_full_cop*N,0) #Full cooperators
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
      x<-c(x,rep(1,round(percentage_full_def,0)), rep(6, round(percentage_full_cop,0)))  # Category 1 y donate 0, category 6 is donate 10
      
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

replications<-100
dataModel4 <- run_model4(replications,100)
dataModel4_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel4_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel4_ <- rbind(dataModel4_,dataModel4[dataModel4$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel4_ <- rbind(dataModel4_,dataModel4[dataModel4$replicationNumber==d,]$frequency)
}
colnames(dataModel4_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


distMatrix <- dist(dataModel4_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel4_)
plot(hc, labels=observedLabels, main="")

#num clusters selected
clustnum <- 4
clusterCut <- cutree(hc, clustnum)
table(clusterCut) #Percentaje of replications on each cluster.  94  4 2 1 

#selecciono la replicacion a imprimir
#1
which(clusterCut==1)
#   1   2   3   4   5   6   7   8   9  10  11  12  ·13  14  15  17  18  19  20  ·21  22  23  24  25  ·26  27  28
#  29  30  31  33  34  35  36  ·37  38  39  40  41  42  43  44  45  46  47  48  ·49  50  52  53  54  55  ·56  57
#  58  59  60  61  64  65  66  67  ·68  69  71  72  73  74  ·75  76  77  78  ·79  80  ·81  82  83  84  85  86  87
#  88  89  90  91  92  93  95  96  97  98  99 100

# Solo imprimo cluster 1 con 93% y después 3 clusters de tamaño < 5%
rep=87

#setEPS()
#postscript("Model4_cluster1_94perc.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4[dataModel4$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel4[dataModel4$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel4[dataModel4$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel4[dataModel4$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
#theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
#dev.off()



############# N=1000 ###########
replications<-100
dataModel4_1000 <- run_model4(replications,1000)
dataModel4_1000_ <- data.frame(matrix(ncol = 84, nrow = 0))
colnames(dataModel4_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')

dataModel4_1000_ <- rbind(dataModel4_1000_,dataModel4_1000[dataModel4_1000$replicationNumber==1,]$frequency)
for (d in 2:replications){
  dataModel4_1000_ <- rbind(dataModel4_1000_,dataModel4_1000[dataModel4_1000$replicationNumber==d,]$frequency)
}
colnames(dataModel4_1000_) <- c('1d0','1d2','1d4','1d6','1d8','1d10','2d0','2d2','2d4','2d6','2d8','2d10','3d0','3d2','3d4','3d6','3d8','3d10','4d0','4d2','4d4','4d6','4d8','4d10','5d0','5d2','5d4','5d6','5d8','5d10','6d0','6d2','6d4','6d6','6d8','6d10','7d0','7d2','7d4','7d6','7d8','7d10','8d0','8d2','8d4','8d6','8d8','8d10','9d0','9d2','9d4','9d6','9d8','9d10','10d0','10d2','10d4','10d6','10d8','10d10','11d0','11d2','11d4','11d6','11d8','11d10','12d0','12d2','12d4','12d6','12d8','12d10','13d0','13d2','13d4','13d6','13d8','13d10','14d0','14d2','14d4','14d6','14d8','14d10')


distMatrix <- dist(dataModel4_1000_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel4_1000_)
plot(hc, labels=observedLabels, main="")

#num clusters selected
clustnum <- 2
clusterCut <- cutree(hc, clustnum)
table(clusterCut) #Percentaje of replications on each cluster.  99  1 

which(clusterCut==1) #1
which(clusterCut==2) #66

rep=89

setEPS()
postscript("Model4_1000_cluster1_99perc.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(x=round, y=contribution)) + 
  geom_point(size=20, shape=15, aes(colour = frequency)) + 
  #scale_colour_gradientn(colours = heat.colors(10), trans = "reverse") +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  geom_line(data=dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_point(data=dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(x=round, y=mean))+
  geom_text(data=dataModel4_1000[dataModel4_1000$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)
theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()



#Scatterplots
setEPS()
postscript("Scatter_Bayesianos_JeffreysPrior_Fulls_100_100rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()


setEPS()
postscript("Scatter_Bayesianos_JeffreysPrior_Fulls_100_1000rep.eps", height=4, width=8.5, family="serif",horizontal=FALSE)
ggplot(dataModel4_1000, aes(x=round, y=contribution)) + 
  #geom_point(size=2, shape=16, aes(colour = frequency)) + 
  geom_jitter(size=0.8, width = 0.4, height = 0.7, aes(colour = frequency), alpha = 1) +
  scale_colour_gradient2(low = "yellow", mid= "red", high = "black",midpoint = 0.5, limits=c(0,1)) + 
  scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(limits=c(0, 15), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
  theme_bw() +theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=14), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab("Round")+ ylab("Contribution")
dev.off()



# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round)){
  for (r in 1:min(10,max(dataModel4$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel4[dataModel4$round==t & dataModel4$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel4$replicationNumber))
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.84 SD0.28

#PPGH_2
# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round)){
  for (r in 1:min(10,max(dataModel4$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel4[dataModel4$round==t & dataModel4$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H2[dist_agregada_PGG_H2$round==(t),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel4$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.92 SD0.12



### Comparison with PGG 100
dist_agregada_PGG100<-dcast(counts_todas_rondas[counts_todas_rondas$treatment=="PGG100",], round ~ contribution, value.var = "COUNT", fun.aggregate=sum)
colnames(dist_agregada_PGG100)<-c("round","c0","c2","c4","c6","c8","c10")


# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round)){
  for (r in 1:min(10,max(dataModel4$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel4[dataModel4$round==t & dataModel4$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG100[dist_agregada_PGG100$round==(t),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(10,max(dataModel4$replicationNumber))
summarytests
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.84 SD0.15

#####################################################
# Fisher tests for specific clusters data
distMatrix <- dist(dataModel4_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel3_)
clustnum <- 4
clusterCut <- cutree(hc, clustnum)

### We compare Cluster 1 with PGG_H data
clusterdata <- dataModel4[dataModel4$replicationNumber %in% which(clusterCut==1),]
# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14
for (t in 1:(final_round)){
  for (r in unique(clusterdata$replicationNumber)){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(clusterdata[clusterdata$round==t & clusterdata$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    # We do not run Chi Square tests because samples are small and so the calculation of p-values may be incorrect
    #chisq.test(x=unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),y=unlist(dist_agregada1[dist_agregada1$round==t,-1]))
    #We use Fisher exact test https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
    # We need at least one data in each category, so I approximate the p-value ensuring al categories have at least one data
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    # Cuidado que aquí hay que comparar simulaciones 1-13 con datos de experimento 2-14
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==(t),-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}
summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(length(unique(clusterdata$replicationNumber)),max(dataModel4$replicationNumber))
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)
# 0.82 SD0.27



##### prueba para ver cuánto tarda con una muestra de 50 

# Conducting r x t tests (r replications, t rounds)
tests <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tests)[1]<-"round"
colnames(tests)[2]<-"replication"
colnames(tests)[3]<-"pvalue" # We use LLR p-value https://cran.r-project.org/web/packages/XNomial/vignettes/XNomial.html
colnames(tests)[4]<-"different" #p-value < 0.001
indice<-1 #just to store the data
final_round <- 14

for (t in 1:(final_round)){
  for (r in 1:min(50,max(dataModel4$replicationNumber))){ #10 replications at most to save computational power,
    dist_agregada1 <- dcast(dataModel4[dataModel4$round==t & dataModel4$replicationNumber==r,1:3], round ~ contribution, value.var = "frequency", fun.aggregate=sum)
    colnames(dist_agregada1)<-c("round","c0","c2","c4","c6","c8","c10")
    dist_agregada1<-dist_agregada1*100 #N 100 agents
    dist_agregada1$round <- dist_agregada1$round/100
    dist_agregada1<-replace(dist_agregada1, dist_agregada1==0, 1)
    
    tests[indice,]$round <- t
    tests[indice,]$replication <- r
    invisible(capture.output(tests[indice,]$pvalue <- xmulti(unlist(dist_agregada_PGG_H[dist_agregada_PGG_H$round==t,-1]),unlist(dist_agregada1[dist_agregada1$round==t,-1]))[4])) #invisible(capture.output to avoid function printing
    tests[indice,]$different <- as.numeric(tests[indice,]$pvalue < 0.001)
    indice<- indice+1
  }
}

summarytests <- dcast(tests[,-3], round ~ different, value.var = "different", fun.aggregate=sum)
summarytests<-summarytests[,-2]
colnames(summarytests)<-c('round','different')
summarytests$percentagediff <- summarytests$different / min(50,max(dataModel4$replicationNumber))
round(mean(summarytests$percentagediff),2)
round(sd(summarytests$percentagediff),2)

#0.83 SD0.24

####### Plot of plots


distMatrix <- dist(dataModel4_, method='DTW')
hc <- hclust(distMatrix, method="average")
observedLabels <- rownames(dataModel4_)
clustnum <- 4
clusterCut <- cutree(hc, clustnum)
table(clusterCut)
clusterdata <- dataModel4[dataModel4$replicationNumber %in% which(clusterCut==1),]

sample_indexes <- sample(unique(clusterdata$replicationNumber), 25)
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
      scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
      geom_line(data=sample[sample$replicationNumber==rep,], aes(x=round, y=mean))+
      geom_point(data=sample[sample$replicationNumber==rep,], aes(x=round, y=mean))+
      geom_text(data=sample[sample$replicationNumber==rep,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
      ggtitle(paste("Replication",i))+
      theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
    indice<-indice+1
}
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,25), 5, 5, byrow=TRUE))


setEPS()
postscript("Examples_Model4_100_cluster1.eps", height=4*5, width=8.5*5, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,25), 5, 5, byrow=TRUE))
dev.off()



############ CRISIS MAXIMA, NO ME GUSTAN LOS CLUSTERS QUE SALEN, TODO MEZCLADO, VOY A PROBAR CLUSTERING JERARQUICO

dist_mat <- dist(dataModel4_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

suppressPackageStartupMessages(library(dendextend))

avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 11)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 11)
table(clusters)


clusterdata <- dataModel4[dataModel4$replicationNumber %in% which(clusters==6),]

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
    scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}


setEPS()
postscript("Examples_Model4_100_cluster6_hierarchicalClust_average.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

### PARECE QUE POR FIN ESTE MÉTODO TIENE SENTIDO


#### complete linkage


dist_mat <- dist(dataModel4_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'complete')
plot(hclust_avg)

suppressPackageStartupMessages(library(dendextend))

avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 5)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 5)
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
    scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}


setEPS()
postscript("Examples_Model4_100_cluster1_hierarchicalClust_complete.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()



#### ward.D2 linkage


dist_mat <- dist(dataModel4_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'ward.D2')
plot(hclust_avg)

suppressPackageStartupMessages(library(dendextend))

avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 4)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 2)
table(clusters)


clusterdata <- dataModel4[dataModel4$replicationNumber %in% which(clusters==2),]

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
    scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}


setEPS()
postscript("Examples_Model4_100_cluster2_hierarchicalClust_ward.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()





#### single linkage


dist_mat <- dist(dataModel4_, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'single')
plot(hclust_avg)

suppressPackageStartupMessages(library(dendextend))

avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 4)
plot(avg_col_dend)

clusters <- cutree(hclust_avg, k = 2)
table(clusters)


clusterdata <- dataModel4[dataModel4$replicationNumber %in% which(clusters==2),]

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
    scale_x_continuous(limits=c(1, 14), breaks=c(1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14)) +
    geom_line(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_point(data=sample[sample$replicationNumber==i,], aes(x=round, y=mean))+
    geom_text(data=sample[sample$replicationNumber==i,], aes(label=round(mean,2),y=0.5+mean,x=round), cex=3)+
    ggtitle(paste("Replication",i))+
    theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=22),legend.text=element_text(size=18), legend.title=element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title=element_text(family='', colour='black', size=14, margin=margin(t=0,b=0, l=50))) + xlab("Round")+ ylab("Contribution")+ theme(legend.position = "none")
  indice<-indice+1
}


setEPS()
postscript("Examples_Model4_100_cluster2_hierarchicalClust_single.eps", height=4*8, width=8.5*4, family="serif",horizontal=FALSE)
grid.arrange(grobs=lista_plots, layout_matrix= matrix(seq(1,32), 8, 4, byrow=TRUE))
dev.off()

