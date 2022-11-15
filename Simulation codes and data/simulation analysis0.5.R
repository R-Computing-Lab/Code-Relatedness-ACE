setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("07121000single0.5.RData")
### 0608 single simulations negtive models

df_neg_perc_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                              names(l.finaldf_0712[["modelset1"]])))

df_neg_perc_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                              names(l.finaldf_0712[["modelset1"]])))

df_neg_perc_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                              names(l.finaldf_0712[["modelset1"]])))
df_neg_perc_one <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                                names(l.finaldf_0712[["modelset1"]])))

for (k in 1:10){
     for (j in 1:10){
          all_a <- numeric()
          all_c <- numeric()
          all_e <- numeric()
          ifneg <- numeric()
          for (i in 1: 1000){
               all_a[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VA"]][j]
               all_c[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VC"]][j]
               all_e[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VE"]][j]
               ifneg[i] <- ifelse(all_a[i] < 0 | all_c[i] < 0 | all_e[i] < 0, .001, 0 )
               
          }
          neg_ace <- c(sum(all_a <= 0)/1000, sum(all_c <= 0)/1000, sum(all_e <= 0)/1000)
          df_neg_perc_a[j,k] <- neg_ace[1]
          df_neg_perc_c[j,k] <- neg_ace[2]
          df_neg_perc_e[j,k] <- neg_ace[3]
          df_neg_perc_one[j,k] <- sum(ifneg)
     }
     
}

### 0529 single simulations negtive models


### correct the function-Done

### Average each cells 0529 single simulations

df_mean_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0712[["modelset1"]])))
df_mean_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0712[["modelset1"]])))
df_mean_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0712[["modelset1"]])))


l.finaldf_0712[["modelset1"]][["N =30"]][["df_ACE"]][1,4] / l.finaldf_0712[["modelset1"]][["N =30"]][["df_ACE"]][1,3]

for (k in 1:10) {
     for (j in 1:10){
          all_a <- numeric()
          all_c <- numeric()
          all_e <- numeric()
          for (i in 1: 1000){
               all_a[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VA"]][j] / l.finaldf_0712[[i]][[k]][["df_ACE"]][["V"]][j]
               all_c[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VC"]][j] / l.finaldf_0712[[i]][[k]][["df_ACE"]][["V"]][j]
               all_e[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VE"]][j] / l.finaldf_0712[[i]][[k]][["df_ACE"]][["V"]][j]
          }
          ace_mean <- c(mean(all_a), mean(all_c), mean(all_e))
          df_mean_a[j,k] <- ace_mean[1]
          df_mean_c[j,k] <- ace_mean[2]
          df_mean_e[j,k] <- ace_mean[3]
     }
}

### 5th percentile and 95th percentile of the 1000 models

df_5th_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                         names(l.finaldf_0712[["modelset1"]])))

df_95th_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                          names(l.finaldf_0712[["modelset1"]])))

df_5th_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                         names(l.finaldf_0712[["modelset1"]])))

df_95th_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                          names(l.finaldf_0712[["modelset1"]])))

df_5th_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                         names(l.finaldf_0712[["modelset1"]])))

df_95th_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                          names(l.finaldf_0712[["modelset1"]])))



for (k in 1:10) {
     for (j in 1:10){
          all_a <- numeric()
          all_c <- numeric()
          all_e <- numeric()
          for (i in 1: 1000){
               all_a[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VA"]][j] 
               all_c[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VC"]][j] 
               all_e[i] <- l.finaldf_0712[[i]][[k]][["df_ACE"]][["VE"]][j]
          }
          ace_5th <- c(quantile(all_a, c(.05)), quantile(all_c, c(.05)), quantile(all_e, c(.05)))
          df_5th_a[j,k] <- ace_5th[1]
          df_5th_c[j,k] <- ace_5th[2]
          df_5th_e[j,k] <- ace_5th[3]
          ace_95th <- c(quantile(all_a, c(.95)), quantile(all_c, c(.95)), quantile(all_e, c(.95)))
          df_95th_a[j,k] <- ace_95th[1]
          df_95th_c[j,k] <- ace_95th[2]
          df_95th_e[j,k] <- ace_95th[3]
     }
}


### If ace model has the lowest AIC compared to ae model and ce model

for(k in 1:10) {
     for (j in 1:10){
          for (i in 1:1000){
               AICs <- l.finaldf_0712[[i]][[k]][["Nest"]][[j]]$AIC[1:3]
               if(which.min(AICs)==1){ch <- "ace"}
               if(which.min(AICs)==2){ch <- "ae"}
               if(which.min(AICs)==3){ch <- "ce"}
               l.finaldf_0712[[i]][[k]][["Nest"]][[j]]$lowest <- c(ch,ch,ch,ch)
          }
     }
}

df_ace_lowest <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                              names(l.finaldf_0712[["modelset1"]])))

for (k in 1:10){
     for (j in 1:10){
          all_aic <- as.character()
          for (i in 1: 1000){
               all_aic[i] <- l.finaldf_0712[[i]][[k]][["Nest"]][[j]]$lowest[1]
          }
          lowest_ace <- c(sum(all_aic == "ace")/1000)
          df_ace_lowest[j,k] <- lowest_ace
     }
     
}


## Power calculation
df_power <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0712[["modelset1"]][["N =60"]][["Nest"]]),
                                                         names(l.finaldf_0712[["modelset1"]])))

for(k in 1:10){
     for(j in 1:10){
          N <- substr(names(l.finaldf_0712[["modelset1"]])[j],4,nchar(names(l.finaldf_0712[["modelset1"]]))[j]) |> as.numeric()
          #print(N)
          DiffLL <- numeric()
          for(i in 1:1000){
               DiffLL[i] <- l.finaldf_0712[[i]][[k]][["Nest"]][[j]]$diffLL[3]
          }
          meanDiffLL <- mean(DiffLL)
          #print(meanDiffLL)
          df_power[j,k] <- 1- pchisq(qchisq(1-.05, 1), 1, meanDiffLL)
     }
}

#Visualize the matrixs
library(plot.matrix)
library("viridis")  
library(Cairo)
par(mar=c(5.1, 4.1, 4.1, 4.1))

Cairo( 840, 630,file="0.5neg.png", type="png", bg="white")
plot(df_neg_perc_one, 
     breaks =c(0,.05,.10,.15,.25,.40,.60,.80,.96,.98,.99,1), 
     fmt.cell='%.3f', 
     digits=3, 
     text.cell=list(cex=1.5), 
     col = magma(11, direction = 1, alpha = 1, begin = 1, end = .5),
     border=NA,
     key = NULL,
     axis.col=list(cex.axis=1.4), axis.row=list(cex.axis=1.08),
     xlab="Sample Size", ylab="R of SS twins",
     cex = 3,
     main= NA,
     font.lab = 2
)
dev.off() 

Cairo( 840, 630,file="0.5AIC.png", type="png", bg="white")
plot(df_ace_lowest, 
     breaks =c(0,.05,.10,.15,.25,.35,.45,.55,.65,.80,.90,1), 
     fmt.cell='%.3f', 
     digits=3, 
     text.cell=list(cex=1.5), 
     col = magma(11, direction = -1, alpha = 1, begin = 1, end = .5),
     border=NA,
     key = NULL,
     axis.col=list(cex.axis=1.4), axis.row=list(cex.axis=1.08),
     xlab="Sample Size", ylab="R of SS twins",
     cex = 3,
     main= NA,
     font.lab = 2
)
dev.off()

Cairo( 840, 630,file="0.5power.png", type="png", bg="white")
plot(df_power, 
     breaks =c(0,.05,.10,.15,.25,.40,.60,.80,.96,.98,.99,1), 
     fmt.cell='%.3f', 
     digits=3, 
     text.cell=list(cex=1.5), 
     col = magma(11, direction = -1, alpha = 1, begin = 1, end = .5),
     border=NA,
     key = NULL,
     axis.col=list(cex.axis=1.4), axis.row=list(cex.axis=1.08),
     xlab="Sample Size", ylab="R of SS twins",
     cex = 3,
     main= NA,
     font.lab = 2
)
dev.off()
