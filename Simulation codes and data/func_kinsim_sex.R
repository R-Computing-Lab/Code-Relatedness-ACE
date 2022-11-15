setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


rmvn <- function(n,sigma) {
     Sh <- with(svd(sigma),
                v%*%diag(sqrt(d))%*%t(u))
     matrix(stats::rnorm(ncol(sigma)*n),
            ncol = ncol(sigma))%*%Sh
}

kinsim_single <- function(
     r=1,
     r_c = 1,
     npg=100,
     mu=0,
     ace=c(1,1,1),
     ...){
     
     sA <- ace[1]^0.5
     sC <- ace[2]^0.5
     sE <- ace[3]^0.5
     
     S2 <- matrix(c(0,1,
                    1,0),2)
     datalist <- list()
     
     
     id <- 1:sum(npg)
     
     n = npg
     
     A.r <- sA*rmvn(n,
                    sigma = diag(2) + S2*r)
     C.r <- sC*rmvn(n,
                    sigma = diag(2) + S2*r_c)
     #C.r <- cbind(C.r,
     #             C.r)
     E.r <- cbind(stats::rnorm(n,
                               sd = sE),
                  stats::rnorm(n,
                               sd = sE))
     
     y.r <- mu + A.r + C.r + E.r
     
     
     r_ <- rep(r,n)
     
     data.r <- data.frame(A.r,C.r,E.r,y.r,r_)
     names(data.r) <- c("A1","A2","C1","C2","E1","E2","y1","y2","r")
     
     
     merged.data.frame <- data.r
     merged.data.frame$id <- id
     
     return(merged.data.frame)
}

diff_size <- function(Ngroup1, Ngroup2, rel=c(1,.5), ace=c(1,1,1), mu = 0){
     df_N1 <- kinsim_single(r = rel[1], npg = Ngroup1, mu = mu, ace = ace)
     df_N2 <- kinsim_single(r = rel[2], npg = Ngroup2, mu = mu, ace = ace)
     df_final <- rbind(df_N1, df_N2)
     return(df_final)
     
     
}


kinsim_comb <- function(N_AM, N_DZ, rel = .75, ace=c(1,1,1), mu = 0, r_c = 1){
     if(rel != 1){
          N_MZ <- round((rel-.5)*2*N_AM)
          N_AM_DZ <- N_AM - N_MZ
          
          df_OSDZ <- kinsim_single(r = .5, npg = N_DZ, mu = mu, ace = ace, r_c = r_c)
          df_SSDZ <- kinsim_single(r = .5, npg = N_AM_DZ, mu = mu, ace = ace, r_c = 1)

          df_MZ <- kinsim_single(r = 1, npg = N_MZ, mu = mu, ace = ace)
          
          df_DZ <- df_OSDZ
          df_DZ$tag <- "DZ"
          
          df_AM <- rbind(df_SSDZ,
                         df_MZ)
          df_AM$tag <- "AM"
          
     }
     if(rel == 1){
          df_DZ <- kinsim_single(r = .5, npg = N_DZ, mu = mu, ace = ace, r_c = r_c)
          df_DZ$tag <- "DZ"
          
          df_AM <- kinsim_single(r = 1, npg = N_AM, mu = mu, ace = ace)
          df_AM$tag <- "AM"
          
     }
     df_final <- rbind(df_AM,df_DZ) 
     return(df_final)
     
     
}

sep_kin <- function(data, basis = "rel"){
     if(basis == "rel"){
          two_r <- sort(unique(data[,9]),decreasing = TRUE)
          data_one <- data[which(data[,9]==two_r[1]),]
          data_two <- data[which(data[,9]==two_r[2]),]
          
          data_one <- data_one[,7:8]
          data_two <- data_two[,7:8]
          
          l_twodata <- list(data_one,data_two)
          return(l_twodata)
     }
     if(basis == "tag"){
          two_t <- sort(unique(data[,11]))
          data_onet <- data[which(data[,11]==two_t[1]),]
          data_twot <- data[which(data[,11]==two_t[2]),]
          
          data_onet <- data_onet[,7:8]
          data_twot <- data_twot[,7:8]
          
          l_twodata <- list(data_onet,data_twot)
          return(l_twodata)
     }
     
}

fit_siACE <- function(data_one, data_two, coe_am, coe_OS, elbound = FALSE){
     # Load Libraries & Options
     require(OpenMx)
     #require(psych)
     #require(polycor)
     # source("miFunctions.R")
     # # Create Output
     # filename <- "oneACEc"
     # sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)
     
     # ----------------------------------------------------------------------------------------------------------------------
     # PREPARE DATA
     # Load Data
     mzData    <- data_one
     dzData    <- data_two
     xxx <- mxMatrix(type = "Full", nrow = 1, ncol = 1,free = FALSE, values = coe_am, name = "coeAM")
     OSOS <- mxMatrix(type = "Full", nrow = 1, ncol = 1,free = FALSE, values = coe_OS, name = "coeOS")
     #coeAM <- coe_am
     
     # covMZ <- cov(mzData, use = "pairwise")
     # covDZ <- cov(dzData, use = "pairwise")
     # # 
     # mean(rbind(mzData,dzData)[,1], na.rm = TRUE)
     
     nv <- 1
     ntv <- 2
     selVars   <- c('y1','y2')
     
     #start values
     svBe <- .01
     svMu <- 0
     svVa <- .2
     svVe <- .5
     
     #variance matrix
     covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVa, label = "VA11", name = "VA")
     covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVa, label = "VC11", name = "VC")
     
     if(elbound == TRUE){
          covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVe, lbound = .001, label = "VE11", name = "VE")
     }else{
          covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv,free = TRUE, values = svVe, label = "VE11", name = "VE")
     }
     
     
     #expected variance matrix
     covP <- mxAlgebra(expression = VA+VC+VE, name = "V")
     covMZ <- mxAlgebra(expression = coeAM*VA+VC, name = "cMZ")
     covDZ <- mxAlgebra(expression = 0.5%x%VA+coeOS*VC, name = "cDZ")
     expCovMz <- mxAlgebra(expression = rbind(cbind(V,cMZ), cbind(t(cMZ),V)), name = "expCovMz")
     expCovDz <- mxAlgebra(expression = rbind(cbind(V,cDZ), cbind(t(cDZ),V)), name = "expCovDz")
     
     #create data
     dataMZ       <- mxData( observed=mzData, type="raw" )
     dataDZ       <- mxData( observed=dzData, type="raw" )
     
     # Mean Matrix
     intercept <- mxMatrix(type = "Full", nrow= 1 , ncol = ntv, free = TRUE, values = 0, labels = "interC", name = "intercept")
     expMean <- mxAlgebra(expression = 1*intercept , name = "expMean")
     
     # Create expectation objects
     expMZ <- mxExpectationNormal(covariance = "expCovMz", means ="expMean", dimnames = selVars)
     expDZ <- mxExpectationNormal(covariance = "expCovDz", means ="expMean", dimnames = selVars)
     funML <- mxFitFunctionML()
     
     #Create models
     pars <- list(intercept, covA, covC, covE, covP)
     modelMZ <- mxModel(pars, expMean,covMZ,expCovMz,dataMZ,expMZ,funML,xxx,name = "MZ")
     #MZfit <- mxRun(modelMZ, intervals = TRUE)
     #summary(MZfit)
     modelDZ <- mxModel(pars, expMean,covDZ,expCovDz,dataDZ,expDZ,funML,OSOS,name = "DZ")
     #DZfit <- mxRun(modelDZ, intervals = TRUE)
     #summary(DZfit)
     
     multi <- mxFitFunctionMultigroup(c("MZ","DZ"))
     
     
     #Algebra for Variance components
     rowUS <- rep("US",nv)
     colUS <- rep(c("VA","VC","VE","SA","SC","SE"),each = nv)
     estUS <- mxAlgebra(expression = cbind(VA,VC,VE,VA/V,VC/V,VE/V), name = "US", dimnames = list(rowUS,colUS))
     
     #CI
     ciACE <- mxCI("US[1,1:6]")
     modelACE <- mxModel("oneACEvc_1cov", pars, modelMZ, modelDZ, multi,estUS,ciACE)
     fitACE <- mxRun(modelACE, intervals = TRUE, silent = TRUE)
     sumACE <-summary(fitACE)
     #sumACE
     
     # ----------------------------------------------------------------------------------------------------------------------
     # RUN SUBMODELS
     # Run AE model
     modelAE <- mxModel( modelACE, name="oneAEvc" )
     modelAE <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
     fitAE <- mxRun( modelAE, intervals=T, silent = TRUE )
     #fitGofs(fitAE); fitEstCis(fitAE)
     # Run CE model
     modelCE <- mxModel( modelACE, name="oneCEvc" )
     modelCE <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
     modelCE <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
     fitCE <- mxRun( modelCE, intervals=T, silent = TRUE )
     #fitGofs(fitCE); fitEstCis(fitCE)
     # Run E model
     modelE <- mxModel( modelACE, name="oneEvc" )
     modelE <- omxSetParameters( modelE, labels=c("VA11","VC11"), free=FALSE, values=0 )
     fitE <- mxRun( modelE, intervals=T, silent = TRUE )
     #fitGofs(fitE); fitEstCis(fitE)
     # Print Comparative Fit Statistics
     df_nested <- mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) )
     #(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result),4)
     l.modeloutput <- list(df_nested,fitACE)
     return(l.modeloutput)
}




#### Functions to do the simulation:
#### Alert: The nmodel is the parameter to determine the steps from .50 to 1.00. Not the numbers of groups of models being run. 
simulation_whole <- function(nmodel = 50, npg = c(500,500), V_ace = c(1,1,1), mode="comb", r_c = 1){
     if (mode == "direct"){
          l.sim <- list()
          
          for(i in 1: nmodel){
               r_ambi <- .5 + i/(nmodel*2)
               l.sim[[i]] <- diff_size(
                    Ngroup1 = npg[1],
                    Ngroup2 = npg[2],
                    rel = c(r_ambi, 0.5),
                    ace = V_ace
               )
               
          }   
          # put simulated data into the ACE model
          l.outputACE <- list()
          l.outputNest <- list()
          for(i in 1:nmodel){
               r_ambi <- .5 + i/(nmodel*2)
               l_doubledata <- sep_kin(l.sim[[i]])
               fitModel <- fit_siACE(l_doubledata[[1]], l_doubledata[[2]],r_ambi)
               l.outputACE[[i]] <- fitModel[[2]]
               l.outputNest[[i]]<- as.data.frame(fitModel[[1]])[,c(2,4:9)]
               names(l.outputNest)[i] <- paste("r=",r_ambi)
          }
     }
     if(mode == "comb"){
          l.sim <- list()
          
          for(i in 1: nmodel){
               r_ambi <- .5 + i/(nmodel*2)
               l.sim[[i]] <- kinsim_comb(
                    N_AM = npg[1],
                    N_DZ = npg[2],
                    rel = r_ambi,
                    ace = V_ace,
                    r_c = r_c 
               )
               
          }
          # put simulated data into the ACE model
          l.outputACE <- list()
          l.outputNest <- list()
          for(i in 1:nmodel){
               r_ambi <- .5 + i/(nmodel*2)
               l_doubledata <- sep_kin(l.sim[[i]], basis = "tag")
               fitModel <- fit_siACE(l_doubledata[[1]], l_doubledata[[2]],r_ambi, coe_OS = r_c )
               l.outputACE[[i]] <- fitModel[[2]]
               l.outputNest[[i]]<- as.data.frame(fitModel[[1]])[,c(2,4:9)]
               names(l.outputNest)[i] <- paste("r=",r_ambi)
          }
     }
     
     
     # Get the estimated variance components for each model
     df_output_ACE <- data.frame("Model"= 1:nmodel, 
                                 "coeAM"= as.numeric(NA),
                                 "V" = as.numeric(NA),
                                 "VA" = as.numeric(NA),
                                 "VC" = as.numeric(NA),
                                 "VE" = as.numeric(NA)
     )
     for(i in 1:nmodel){
          r_ambi <- .5 + i/(nmodel*2)
          df_output_ACE[i,2] <- r_ambi
          df_output_ACE[i,3] <- l.outputACE[[i]]@algebras$V$result[1]
          df_output_ACE[i,4:6] <- l.outputACE[[i]]@algebras$US$result[1:3]
     }
     # Get more information like se, ci, variance and covariance model status
     l.info <- list()
     for(i in 1:nmodel){
          modelInfo <- list(
               SE = l.outputACE[[i]]@output[["standardErrors"]],
               CI = l.outputACE[[i]]@output[["confidenceIntervals"]],
               cov = l.outputACE[[i]]@output[["vcov"]],
               status = l.outputACE[[i]]@output[["status"]],
               iterations = l.outputACE[[i]]@output[["iterations"]],
               evaluations = l.outputACE[[i]]@output[["evaluations"]]
          )
          l.info[[i]] <- modelInfo
          
     }
     ## Get the nested model comparisons
     names(l.info) <- names(l.outputNest)
     df_output <- list(Nest = l.outputNest, df_ACE = df_output_ACE, info = l.info)
     return(df_output)
}


#07.07.2022 10*8*100 0.55start: test the effects of OSDZ twins, ace=1.5,.6,.9
l.finaldf_0707 <- list()
N_level0707 <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:1000){
     l.finaldf_k <- list()
     set.seed(k)
     for(j in 1: length(N_level0707)){
          n1 = N_level0707[j]
          n2 = N_level0707[j]
          l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(1.5,.6,.9), mode = "comb", r_c = .95)
          names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
     }  
     l.finaldf_0707[[k]] <- l.finaldf_k
     names(l.finaldf_0707)[[k]] <- paste("modelset",k, sep = "")
     print(c(k, format(Sys.time(), "%c")))
}


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("0707single_gender95.RData")


### 0529 single simulations negtive models

df_neg_perc_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =60"]][["Nest"]]),
                                                              names(l.finaldf_0707[["modelset1"]])))

df_neg_perc_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =60"]][["Nest"]]),
                                                              names(l.finaldf_0707[["modelset1"]])))

df_neg_perc_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =60"]][["Nest"]]),
                                                              names(l.finaldf_0707[["modelset1"]])))
df_neg_perc_one <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =60"]][["Nest"]]),
                                                                names(l.finaldf_0707[["modelset1"]])))

for (k in 1:10){
        for (j in 1:10){
                all_a <- numeric()
                all_c <- numeric()
                all_e <- numeric()
                ifneg <- numeric()
                for (i in 1: 1000){
                        all_a[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VA"]][j]
                        all_c[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VC"]][j]
                        all_e[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VE"]][j]
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

df_mean_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0707[["modelset1"]])))
df_mean_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0707[["modelset1"]])))
df_mean_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0707[["modelset1"]])))


l.finaldf_0707[["modelset1"]][["N =30"]][["df_ACE"]][1,4] / l.finaldf_0707[["modelset1"]][["N =30"]][["df_ACE"]][1,3]

for (k in 1:10) {
     for (j in 1:10){
          all_a <- numeric()
          all_c <- numeric()
          all_e <- numeric()
          for (i in 1: 1000){
               all_a[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VA"]][j] / l.finaldf_0707[[i]][[k]][["df_ACE"]][["V"]][j]
               all_c[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VC"]][j] / l.finaldf_0707[[i]][[k]][["df_ACE"]][["V"]][j]
               all_e[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VE"]][j] / l.finaldf_0707[[i]][[k]][["df_ACE"]][["V"]][j]
          }
          ace_mean <- c(mean(all_a), mean(all_c), mean(all_e))
          df_mean_a[j,k] <- ace_mean[1]
          df_mean_c[j,k] <- ace_mean[2]
          df_mean_e[j,k] <- ace_mean[3]
     }
}

### 5th percentile and 95th percentile of the 1000 models

df_5th_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                         names(l.finaldf_0707[["modelset1"]])))

df_95th_a <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0707[["modelset1"]])))

df_5th_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                         names(l.finaldf_0707[["modelset1"]])))

df_95th_c <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0707[["modelset1"]])))

df_5th_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                         names(l.finaldf_0707[["modelset1"]])))

df_95th_e <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                          names(l.finaldf_0707[["modelset1"]])))



for (k in 1:10) {
     for (j in 1:10){
          all_a <- numeric()
          all_c <- numeric()
          all_e <- numeric()
          for (i in 1: 1000){
               all_a[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VA"]][j] 
               all_c[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VC"]][j] 
               all_e[i] <- l.finaldf_0707[[i]][[k]][["df_ACE"]][["VE"]][j]
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

### compare results with natural relatedness 

### diffLL table

### If ace model has the lowest AIC compared to ae model and ce model

for(k in 1:10) {
     for (j in 1:10){
          for (i in 1:1000){
               AICs <- l.finaldf_0707[[i]][[k]][["Nest"]][[j]]$AIC[1:3]
               if(which.min(AICs)==1){ch <- "ace"}
               if(which.min(AICs)==2){ch <- "ae"}
               if(which.min(AICs)==3){ch <- "ce"}
               l.finaldf_0707[[i]][[k]][["Nest"]][[j]]$lowest <- c(ch,ch,ch,ch)
          }
     }
}

df_ace_lowest <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =30"]][["Nest"]]),
                                                              names(l.finaldf_0707[["modelset1"]])))

for (k in 1:10){
     for (j in 1:10){
          all_aic <- as.character()
          for (i in 1: 1000){
               all_aic[i] <- l.finaldf_0707[[i]][[k]][["Nest"]][[j]]$lowest[1]
          }
          lowest_ace <- c(sum(all_aic == "ace")/1000)
          df_ace_lowest[j,k] <- lowest_ace
     }
     
}

df_ae_lowest <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =60"]][["Nest"]]),
                                                             names(l.finaldf_0707[["modelset1"]])))
for (k in 1:10){
        for (j in 1:10){
                all_aic <- as.character()
                for (i in 1: 1000){
                        all_aic[i] <- l.finaldf_0707[[i]][[k]][["Nest"]][[j]]$lowest[1]
                }
                lowest_ace <- c(sum(all_aic == "ae")/1000)
                df_ae_lowest[j,k] <- lowest_ace
        }
        
}

df_ce_lowest <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =60"]][["Nest"]]),
                                                             names(l.finaldf_0707[["modelset1"]])))
for (k in 1:10){
        for (j in 1:10){
                all_aic <- as.character()
                for (i in 1: 1000){
                        all_aic[i] <- l.finaldf_0707[[i]][[k]][["Nest"]][[j]]$lowest[1]
                }
                lowest_ace <- c(sum(all_aic == "ce")/1000)
                df_ce_lowest[j,k] <- lowest_ace
        }
        
}


## Power calculation
df_power <- matrix(nrow = 10, ncol = 10, dimnames = list(names(l.finaldf_0707[["modelset1"]][["N =60"]][["Nest"]]),
                                                         names(l.finaldf_0707[["modelset1"]])))

for(k in 1:10){
        for(j in 1:10){
                N <- substr(names(l.finaldf_0707[["modelset1"]])[j],4,nchar(names(l.finaldf_0707[["modelset1"]]))[j]) |> as.numeric()
                #print(N)
                DiffLL <- numeric()
                for(i in 1:1000){
                        DiffLL[i] <- l.finaldf_0707[[i]][[k]][["Nest"]][[j]]$diffLL[3]
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

Cairo( 840, 630,file="sex1.5neg.png", type="png", bg="white")
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

Cairo( 840, 630,file="sex1.5AIC.png", type="png", bg="white")
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

Cairo( 840, 630,file="sex1.5power.png", type="png", bg="white")
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
