setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## Functions for simulation
source("func_kinsim_single.R")

## testing trials

# l.finaldf_test <- list()
# for(k in 1:5){
#      l.finaldf_k <- list()
#      for(j in 1: 5){
#           seed = as.numeric(paste(k,j,sep = ""))
#           set.seed(seed)
#           n1 = 30*j
#           n2 = 30*j
#           l.finaldf_k[[j]] <- simulation_whole(5,c(n1,n2),c(10,10,10), mode = "comb")
#           names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
#      }
#      l.finaldf_test[[k]] <- l.finaldf_k
#      names(l.finaldf_test)[[k]] <- paste("modelset",k, sep = "")
#      print(c(k, format(Sys.time(), "%c")))
# }


### 
#05.29.2022 10*8*100 5.46start: simulate use new functions
l.finaldf_0529 <- list()
N_level0529 <- c( 30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:100){
     l.finaldf_k <- list()
     set.seed(k)
     for(j in 1: length(N_level0529)){
          n1 = N_level0529[j]
          n2 = N_level0529[j]
          l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(10,10,10), mode = "comb")
          names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
     }  
     l.finaldf_0529[[k]] <- l.finaldf_k
     names(l.finaldf_0529)[[k]] <- paste("modelset",k, sep = "")
     print(c(k, format(Sys.time(), "%c")))
}

### 
#05.30.2022 10*8*100 2.50start: simulate use new functions and ace=(8,1,1)
l.finaldf_0530 <- list()
N_level0530 <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:100){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0530)){
                n1 = N_level0530[j]
                n2 = N_level0530[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(8,1,1), mode = "comb")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0530[[k]] <- l.finaldf_k
        names(l.finaldf_0530)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}


#05.31.2022 10*8*100 0.55start: simulate use new functions and the mode of "direct"
l.finaldf_0531 <- list()
N_level0531 <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:100){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0531)){
                n1 = N_level0531[j]
                n2 = N_level0531[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(10,10,10), mode = "direct")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0531[[k]] <- l.finaldf_k
        names(l.finaldf_0531)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}


#06.01.2022 10*8*100 0.55start: simulate use new functions and the mode of "comb", ace=1,1,1
l.finaldf_0601 <- list()
N_level0601 <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:100){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0601)){
                n1 = N_level0601[j]
                n2 = N_level0601[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(1,1,1), mode = "comb")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0601[[k]] <- l.finaldf_k
        names(l.finaldf_0601)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}

#06.01.2022b 10*8*1000 0.55start: simulate use new functions and the mode of "comb", ace=1,1,1
l.finaldf_0601b <- list()
N_level0601b <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:1000){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0601b)){
                n1 = N_level0601b[j]
                n2 = N_level0601b[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(1,1,1), mode = "comb")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0601b[[k]] <- l.finaldf_k
        names(l.finaldf_0601b)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}


#06.03.2022 10*8*1000 0.55start: simulate use new functions and the mode of "comb", ace=2.4,.3,.3
l.finaldf_0603 <- list()
N_level0603 <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:1000){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0603)){
                n1 = N_level0603[j]
                n2 = N_level0603[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(2.4,.3,.3), mode = "comb")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0603[[k]] <- l.finaldf_k
        names(l.finaldf_0603)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}

#06.08.2022 10*10*1000 0.55start: simulate use new functions and the mode of "comb", ace=1.5,.6,.9
l.finaldf_0608 <- list()
N_level0608 <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)
# left N_levels c()
for(k in 1:1000){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0608)){
                n1 = N_level0608[j]
                n2 = N_level0608[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(1.5,.6,.9), mode = "comb")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0608[[k]] <- l.finaldf_k
        names(l.finaldf_0608)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}


#06.14.2022 10*8*1000 for CFPS BGAstart: simulate use new functions and the mode of "comb", ace=0.22*V, 0.74*V, 0.04*V
l.finaldf_0614 <- list()
V = 2.128281
N_level0614a <- c(30*2,30*4,30*6,30*8,30*10,30*16,30*24,30*30)
N_level0614b <- c(30*1,30*2,30*3,30*4,30*5,30*8,30*12,30*15)
# left N_levels c()
for(k in 1:1000){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0614a)){
                n1 = N_level0614a[j]
                n2 = N_level0614b[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(0.22*V, 0.74*V, 0.04*V), mode = "comb")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0614[[k]] <- l.finaldf_k
        names(l.finaldf_0614)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}

#07.12.2022 10*10*1000 forstart: simulate use new functions and the mode of "comb", ace=0.5, 2, 0.5
l.finaldf_0712 <- list()
V = 3
N_level0712 <- c(30*1,30*2,30*3,30*5,30*7,30*10,30*15,30*25,30*40,30*65)

# left N_levels c()
for(k in 1:1000){
        l.finaldf_k <- list()
        set.seed(k)
        for(j in 1: length(N_level0712)){
                n1 = N_level0712[j]
                n2 = N_level0712[j]
                l.finaldf_k[[j]] <- simulation_whole(10,c(n1,n2),c(0.5*V, 2*V, 0.5*V), mode = "comb")
                names(l.finaldf_k)[[j]] <- paste("N =",n1, sep = "")
        }  
        l.finaldf_0712[[k]] <- l.finaldf_k
        names(l.finaldf_0712)[[k]] <- paste("modelset",k, sep = "")
        print(c(k, format(Sys.time(), "%c")))
}


