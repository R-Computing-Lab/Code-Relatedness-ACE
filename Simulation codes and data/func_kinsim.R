#' Simulate Biometrically informed Multivariate Data



kinsim <- function(
     r_all=c(1,.5),
     npg_all=500,
     npergroup_all=rep(npg_all,length(r_all)),
     mu_all=0,
     variables=2,
     mu_list=rep(mu_all,variables),
     r_vector=NULL, # alternative specification, give vector of rs
     ace_all=c(1,1,1), # variance default
     ace_list=matrix(rep(ace_all,variables),byrow=TRUE,nrow=variables),
     cov_a=0, #default shared covariance for genetics across variables
     cov_c=0, #default shared variance for c across variables
     cov_e=0, #default shared variance for e across variables
     ...){
     mu <- NULL
     sA <- ace_list[,1]^0.5
     sC <- ace_list[,2]^0.5
     sE <- ace_list[,3]^0.5
     S2 <- diag(4)*-1+1
     
     datalist <- list()
     if(variables==1){
          data_v<-kinsim_internal(r = r_all,
                                  npergroup = npergroup_all,	#
                                  mu = mu_list[1],			#intercept
                                  ace = ace_list[[1]],
                                  r_vector = r_vector
          )
          data_v$A1_u<-data_v$A1
          data_v$A2_u<-data_v$A2
          data_v$C1_u<-data_v$C1
          data_v$C2_u<-data_v$C2
          data_v$E1_u<-data_v$E1
          data_v$E2_u<-data_v$E2
          data_v$y1_u<-data_v$y1
          data_v$y2_u<-data_v$y2
          
          merged.data.frame = data_v
          names(merged.data.frame)[c(1,10)]<-c("id","r")
     }
     if(variables > 2){
          stop("You have tried to generate data beyond the current limitations of this program. Maximum variables 2.")
     }
     if(is.null(r_vector)){
          id <- 1:sum(npergroup_all)
          for(i in 1:length(r_all)){
               n = npergroup_all[i]
               
               # Genetic Covariance
               sigma_a <- diag(4) + S2*r_all[i]
               sigma_a[1,3] <- cov_a
               sigma_a[3,1] <- cov_a
               sigma_a[2,4] <- cov_a
               sigma_a[4,2] <- cov_a
               sigma_a[1,4] <- cov_a*r_all[i]
               sigma_a[4,1] <- cov_a*r_all[i]
               sigma_a[3,2] <- cov_a*r_all[i]
               sigma_a[2,3] <- cov_a*r_all[i]
               A.r <- rmvn(n,
                           sigma = sigma_a)
               
               A.r[,1:2]<- A.r[,1:2]*sA[1]
               A.r[,3:4]<- A.r[,3:4]*sA[2]
               
               # Shared C Covariance
               sigma_c<-diag(4)+S2*1
               sigma_c[1,3]<-cov_c
               sigma_c[3,1]<-cov_c
               sigma_c[2,4]<-cov_c
               sigma_c[4,2]<-cov_c
               sigma_c[1,4]<-cov_c*1
               sigma_c[4,1]<-cov_c*1
               sigma_c[3,2]<-cov_c*1
               sigma_c[2,3]<-cov_c*1
               C.r <- rmvn(n,
                           sigma = sigma_c)
               C.r[,1:2]<- C.r[,1:2]*sC[1]
               C.r[,3:4]<- C.r[,3:4]*sC[2]
               
               # Shared E Covariance
               sigma_e <- diag(4)+S2*0
               sigma_e[1,3] <- cov_e
               sigma_e[3,1] <- cov_e
               sigma_e[2,4] <- cov_e
               sigma_e[4,2] <- cov_e
               E.r <- rmvn(n,
                           sigma = sigma_e)
               E.r[,1:2]<- E.r[,1:2]*sE[1]
               E.r[,3:4]<- E.r[,3:4]*sE[2]
               
               # total score
               y.r <-  A.r + C.r + E.r
               
               
               y.r[,1:2] <- y.r[,1:2] + mu_list[1]
               y.r[,3:4] <- y.r[,3:4] + mu_list[2]
               r_ <- rep(r_all[i],
                         n)
               
               data.r<-data.frame(A.r,C.r,E.r,y.r,r_)
               names(data.r)<-c("A1_1","A1_2",
                                "A2_1","A2_2",
                                "C1_1","C1_2",
                                "C2_1","C2_2",
                                "E1_1","E1_2",
                                "E2_1","E2_2",
                                "y1_1","y1_2",
                                "y2_1","y2_2",
                                "r")
               
               datalist[[i]] <- data.r
               names(datalist)[i] <- paste0("datar",r_all[i])
          }
          merged.data.frame <- Reduce(function(...) merge(..., all=T), datalist)
          merged.data.frame$id <- id
     }else{
          id=1:length(r_vector)
          data_vector <- data.frame(id,
                                    r_vector,
                                    matrix(rep(as.numeric(NA),
                                               length(id)*4),
                                           nrow=length(id),
                                           ncol=4))
          
          names(data_vector) <- c("id","r",
                                  "A1_1","A1_2",
                                  "A2_1","A2_2")
          
          unique_r= matrix(unique(r_vector))
          
          for(i in 1:length(unique_r)){
               n <- length(r_vector[r_vector==unique_r[i]])
               
               # Genetic Covariance
               sigma_a<-diag(4)+S2*unique_r[i]
               sigma_a[1,3]<-cov_a
               sigma_a[3,1]<-cov_a
               sigma_a[2,4]<-cov_a
               sigma_a[4,2]<-cov_a
               sigma_a[1,4]<-cov_a*unique_r[i]
               sigma_a[4,1]<-cov_a*unique_r[i]
               sigma_a[3,2]<-cov_a*unique_r[i]
               sigma_a[2,3]<-cov_a*unique_r[i]
               A.r <- rmvn(n,
                           sigma = sigma_a)
               data_vector$A1_1[data_vector$r_vector==unique_r[i]] <- A.r[,1]*sA[1]
               data_vector$A1_2[data_vector$r_vector==unique_r[i]] <- A.r[,2]*sA[1]
               data_vector$A2_1[data_vector$r_vector==unique_r[i]] <- A.r[,3]*sA[2]
               data_vector$A2_2[data_vector$r_vector==unique_r[i]] <- A.r[,4]*sA[2]
               A.r[,1:2]<- A.r[,1:2]
               A.r[,3:4]<- A.r[,3:4]*sA[2]
          }
          n <- length(r_vector)
          A.r <- matrix(c(data_vector$A1_1,
                          data_vector$A1_2,
                          data_vector$A2_1,
                          data_vector$A2_2),
                        ncol=4,
                        nrow=n)
          # Shared C Covariance
          sigma_c<-diag(4)+S2*1
          sigma_c[1,3]<-cov_c;sigma_c[3,1]<-cov_c;sigma_c[2,4]<-cov_c;sigma_c[4,2]<-cov_c
          sigma_c[1,4]<-cov_c*1;sigma_c[4,1]<-cov_c*1;sigma_c[3,2]<-cov_c*1;sigma_c[2,3]<-cov_c*1
          C.r <- rmvn(n,sigma=sigma_c)
          C.r[,1:2]<- C.r[,1:2]*sC[1]; C.r[,3:4]<- C.r[,3:4]*sC[2]
          
          # Shared E Covariance
          sigma_e<-diag(4)+S2*0
          sigma_e[1,3]<-cov_e;sigma_e[3,1]<-cov_e;sigma_e[2,4]<-cov_e;sigma_e[4,2]<-cov_e
          E.r <- rmvn(n,sigma=sigma_e)
          E.r[,1:2]<- E.r[,1:2]*sE[1]; E.r[,3:4]<- E.r[,3:4]*sE[2]
          
          
          y.r <- A.r
          y.r[,1:2]<-A.r[,1:2]*ace_list[1,1] + C.r[,1:2]*ace_list[1,2] + E.r[,1:2]*ace_list[1,3]
          y.r[,3:4]<-A.r[,3:4]*ace_list[2,1] + C.r[,3:4]*ace_list[2,2] + E.r[,3:4]*ace_list[2,3]
          y.r[,1:2]<-y.r[,1:2]+mu_list[1]
          y.r[,3:4]<-y.r[,3:4]+mu_list[2]
          y.r <- mu + A.r + C.r + E.r
          data.r<-data.frame(A.r,C.r,E.r,y.r,r_vector,id)
          names(data.r)<-c("A1_1","A1_2","A2_1","A2_2","C1_1","C1_2","C2_1","C2_2","E1_1","E1_2","E2_1","E2_2","y1_1","y1_2","y2_1","y2_2","r","id")
          
          
          datalist[[i]] <- data.r
          names(datalist)[i]<-paste0("datar",r_all[i])
          merged.data.frame = data.r
     }
     return(merged.data.frame)
}


rmvn <- function(n,sigma) {
     Sh <- with(svd(sigma),
                v%*%diag(sqrt(d))%*%t(u))
     matrix(stats::rnorm(ncol(sigma)*n),
            ncol = ncol(sigma))%*%Sh
}

kinsim_internal <- function(
     r=c(1,.5),
     npg=100,
     npergroup=rep(npg,length(r)),
     mu=0,
     ace=c(1,1,1),
     r_vector=NULL,
     ...){
     
     sA <- ace[1]^0.5
     sC <- ace[2]^0.5
     sE <- ace[3]^0.5
     
     S2 <- matrix(c(0,1,
                    1,0),2)
     datalist <- list()
     
     if(is.null(r_vector)){
          
          id <- 1:sum(npergroup)
          
          for(i in 1:length(r)){
               
               n = npergroup[i]
               
               A.r <- sA*rmvn(n,
                              sigma = diag(2) + S2*r[i])
               C.r <- stats::rnorm(n,
                                   sd = sC)
               C.r <- cbind(C.r,
                            C.r)
               E.r <- cbind(stats::rnorm(n,
                                         sd = sE),
                            stats::rnorm(n,
                                         sd = sE))
               
               y.r <- mu + A.r + C.r + E.r
               
               
               r_ <- rep(r[i],n)
               
               data.r <- data.frame(A.r,C.r,E.r,y.r,r_)
               names(data.r) <- c("A1","A2","C1","C2","E1","E2","y1","y2","r")
               datalist[[i]] <- data.r
               names(datalist)[i] <- paste0("datar",r[i])
          }
          merged.data.frame <- Reduce(function(...) merge(..., all=T), datalist)
          merged.data.frame$id <- id
     }else{
          id=1:length(r_vector)
          data_vector=data.frame(id,r_vector)
          data_vector$A.r1<-as.numeric(NA)
          data_vector$A.r2<-as.numeric(NA)
          unique_r= matrix(unique(r_vector))
          for(i in 1:length(unique_r)){
               n <- length(r_vector[r_vector==unique_r[i]])
               A.rz <- sA*rmvn(n,
                               sigma=diag(2)+S2*unique_r[i])
               data_vector$A.r1[data_vector$r_vector==unique_r[i]] <- A.rz[,1]
               data_vector$A.r2[data_vector$r_vector==unique_r[i]] <- A.rz[,2]
          }
          n=length(r_vector)
          A.r <- matrix(c(data_vector$A.r1,
                          data_vector$A.r2),ncol=2)
          C.r <- stats::rnorm(n,sd=sC)
          C.r <- cbind(C.r,C.r)
          E.r <- cbind(stats::rnorm(n,
                                    sd = sE),
                       stats::rnorm(n,
                                    sd = sE))
          
          y.r <- mu + A.r + C.r + E.r
          
          data.r <- data.frame(id,A.r,C.r,E.r,y.r,r_vector)
          names(data.r) <- c("id","A1","A2","C1","C2","E1","E2","y1","y2","r")
          datalist[[i]] <- data.r
          names(datalist)[i]<- paste0("datar",r[i])
          
          merged.data.frame <- data.r
     }
     
     return(merged.data.frame)
}

# Create simulated kin data to achieve the unusual relatedness by combining 1.0 and .5 relatedness twins. 
kinsim_comb <- function(N_AM, N_DZ, rel = .75, ace=c(1,1,1), mu = 0){
        #DZ data generation
        df_DZ <- kinsim(r_all = c(1,.5),
                          npg_all = N_DZ,
                          ace_all = ace,
                          variable = 1,
                          mu_all = mu)
        df_DZ <- df_DZ[which(df_DZ[,17]==sort(unique(df_DZ[,17]),decreasing = TRUE)[2]),]
        df_DZ$tag <- "DZ"
        
        #AM data generation
        if(rel != 1){
                N1 <- round((1-rel)*2*N_AM)
                N2 <- N_AM - N1
                df_AM <- diff_size(N1,N2,rel=c(1,.5), ace = ace, mu = mu)
                df_AM$tag <- "AM"
                
                 
        }
        if(rel == 1){
                df_AM <- kinsim(r_all = c(1,.5),
                                npg_all = N_AM,
                                ace_all = ace,
                                variable = 1,
                                mu_all = mu)
                df_AM <- df_AM[which(df_AM[,17]==sort(unique(df_AM[,17]),decreasing = TRUE)[1]),]
                df_AM$tag <- "AM"
        }
        df_final <- rbind(df_AM,df_DZ) 
        return(df_final)
        
        
}

### A wrapper to create different sample sizes for each group
diff_size <- function(Ngroup1, Ngroup2, rel=c(1,.5), ace=c(1,1,1), mu = 0){
        if(Ngroup1 == Ngroup2){
                df_final <- kinsim(r_all = rel,
                                   npg_all = Ngroup1,
                                   ace_all = ace,
                                   variable = 1,
                                   mu_all = mu)
        }
        if(Ngroup1 > Ngroup2){
                N1 <- Ngroup2
                N2 <- Ngroup1-Ngroup2
                df_core <- kinsim(r_all = rel,
                                  npg_all = N1,
                                  ace_all = ace,
                                  variable = 1,
                                  mu_all = mu)
                df_supp <- kinsim(r_all = rel,
                                  npg_all = N2,
                                  ace_all = ace,
                                  variable = 1,
                                  mu_all = mu)
                two_r <- sort(unique(df_supp[,17]),decreasing = TRUE)
                df_att <- df_supp[which(df_supp[,17]==two_r[1]),]
                df_final <- rbind(df_core, df_att)
        }
        if(Ngroup2 > Ngroup1){
                N1 <- Ngroup1
                N2 <- Ngroup2-Ngroup1
                df_core <- kinsim(r_all = rel,
                                  npg_all = N1,
                                  ace_all = ace,
                                  variable = 1,
                                  mu_all = mu)
                df_supp <- kinsim(r_all = rel,
                                  npg_all = N2,
                                  ace_all = ace,
                                  variable = 1,
                                  mu_all = mu)
                two_r <- sort(unique(df_supp[,17]),decreasing = TRUE)
                df_att <- df_supp[which(df_supp[,17]==two_r[2]),]
                df_final <- rbind(df_core, df_att)
        }
        return(df_final)
        
        
}


sep_kin <- function(data, basis = "rel"){
     if(basis == "rel"){
             two_r <- sort(unique(data[,17]),decreasing = TRUE)
             data_one <- data[which(data[,17]==two_r[1]),]
             data_two <- data[which(data[,17]==two_r[2]),]
             
             data_one <- data_one[,13:14]
             data_two <- data_two[,13:14]
             
             l_twodata <- list(data_one,data_two)
             return(l_twodata)
     }
     if(basis == "tag"){
             two_t <- sort(unique(data[,19]))
             data_onet <- data[which(data[,19]==two_t[1]),]
             data_twot <- data[which(data[,19]==two_t[2]),]
             
             data_onet <- data_onet[,13:14]
             data_twot <- data_twot[,13:14]
             
             l_twodata <- list(data_onet,data_twot)
             return(l_twodata)
     }

}



fit_siACE <- function(data_one, data_two, coe_am, elbound = FALSE){
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
     #coeAM <- coe_am
     
     # covMZ <- cov(mzData, use = "pairwise")
     # covDZ <- cov(dzData, use = "pairwise")
     # # 
     # mean(rbind(mzData,dzData)[,1], na.rm = TRUE)

     nv <- 1
     ntv <- 2
     selVars   <- c('y1_1','y1_2')
     
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
     covDZ <- mxAlgebra(expression = 0.5%x%VA+VC, name = "cDZ")
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
     modelDZ <- mxModel(pars, expMean,covDZ,expCovDz,dataDZ,expDZ,funML,name = "DZ")
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
     #sumACE <-summary(fitACE)
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
simulation_whole <- function(nmodel = 50, npg = c(500,500), V_ace = c(1,1,1), mode="direct"){
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
                                ace = V_ace
                        )

                }
                # put simulated data into the ACE model
                l.outputACE <- list()
                l.outputNest <- list()
                for(i in 1:nmodel){
                        r_ambi <- .5 + i/(nmodel*2)
                        l_doubledata <- sep_kin(l.sim[[i]], basis = "tag")
                        fitModel <- fit_siACE(l_doubledata[[1]], l_doubledata[[2]],r_ambi)
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


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        require(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}
