### Deprecated 

## Fisher tests. Reasons not to include:
# - Low statistical power since it is high computationally expensive, we need to run 14 tests per replication (one per round)


########### MODEL A: Bayesian agents with non-informative prior ############################################

#Old text for model A:
#For each of a total of ten replications of the model and for each round, we run Fisher's exact tests (one per round) to compare these simulation results (14000 data points) to those decisions of people in the experiment with the same number of people, 100 in this case. 80\% of the tests found empirical differences in the comparison of the simulated data to the PGG\_H data and 84\% to the PGG\_H2 data.

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


########### MODEL B: Bayesian agents with data-based informative prior ############################################

