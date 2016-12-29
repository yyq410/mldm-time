setwd("~/LDM-time/kmLDM/")
#setwd("F://R code/west english channel time/")

library(lbfgs)
library(QUIC)
library(mclust)
library(huge)

source("./program/generateSimulation.R")
args <- commandArgs()
cluster <- as.numeric(args[6])

testTrue <- FALSE
Toy <- TRUE
testOther <- FALSE
addRandom <- FALSE
showFlag <- FALSE

dir <- 'simulation-test-new-best-varyCluster-AGP/'
dir.create(dir)

graph_list <- list()

if(Toy){
  #load("./toyData.RData")
  #load("./ag_noise_mclust.RData")
  load("./ag_otu_meta_547_7_5305.RData")
  X <- x_final
  M <- m_final
  p <- ncol(X)
  q <- ncol(M)
  K <- cluster
  n <- nrow(X)
  t <- 1
  #X <- X[(n+1):(2*n),]
  #print(cor(X))
  #M <- as.matrix(M[(n+1):(2*n),])
  #print(cor(X,M))
  #K <- 1
  #trueParameters <- list(trueParameters[[2]])
} else {
X <- matrix(0, sum(n_list), p)
M <- matrix(0, sum(n_list), q)

basis <- 0
trueParameters <- list()
  for(i in 1:K) {
    n <- n_list[i]
    graph <- sample(c(1,2,3,4), 1)
    res <- generateSimulation(n, p, q, graph, i)
    graph_list[[i]] <- res[[4]]
    X[(basis + 1):(basis + n), ] <- res[[1]]
    M[(basis + 1):(basis + n), ] <- res[[2]]
    print("colMeans X")
    print(colMeans(res[[1]]))
    print("colMeans M")
    print(colMeans(res[[2]]))
    trueParameters[[i]] <- res[[3]]
  if(showFlag) {
    print("B:")
    print(trueParameters[[i]][[1]])
    print("Theta:")
    print(trueParameters[[i]][[2]])
    print("Sigma:")
    print(solve(trueParameters[[i]][[2]]))
    print("B0:")
    print(trueParameters[[i]][[3]])
    print("muM:")
    print(trueParameters[[i]][[4]])
    print("covM:")
    print(trueParameters[[i]][[5]])
  }
    basis <- basis + n
  }
}

graph <- "AGP"
subDir <- paste(dir, 'p-', as.character(p), '-q-', as.character(q), '-n-', as.character(n),'-K-', as.character(K), '-', graph, "/",sep="")
dir.create(subDir)

Mscale <- scale(M, center = TRUE, scale = TRUE)

otherResults <- list()
if(testOther) {
  source("./program/testOtherMethods.R")
  otherResults <- testOtherMethods(X, M, subDir, t)
}


# Cluster datasets via meta data, and initialization
if(K == 2) {
  load("./m_final_cluater2.RData")
  model <- mod
} else {
  model <- Mclust(M, G = K)
}
classification <- model$classification

print("mclust classification:")
print(table(classification))

Z_mean <- 1
Z_init <- log(X+1) + Z_mean
print('Z init max min mean:')
print(max(Z_init))
print(min(Z_init))
print(mean(Z_init))

# record parameters for every cluster
initParameters <- list()
# record the weight of every cluster
initPi <- table(classification) / n
print("init Pi:")
print(initPi)
# Record lambda1 and lambda2
corx_list <- c()
corm_list <- c()

kmLDM_result <- list()
kmLDM_result_all <- list()
source("./program/mLDM-B-Theta-noScaleM.R")
source("./program/Lognormal-Dirichlet-Multinomial-lbfgs-proximal-split-q-active-set-quic-B-Theta-speed.R")

print("##########################################")
print("######### mLDM initialization ############")
print("##########################################")

Zmldm <- matrix(0, n, p)
res <- list()
for(i in 1:K) {
  print(paste("init Cluster : ", as.character(i)))
  perIndex <- which(classification == i)
  Xs <- as.matrix(X[perIndex,])
  Ms <- as.matrix(Mscale[perIndex,])
  res[[i]] <- mLDM(X = Xs, M = Ms, model_selection_num = 1, verbose = TRUE, max_iteration = 1, threshold = 1e-2, max_iteration_B = 50)  
  resi <- res[[i]][[1]]
  initParameters[[i]] <- list(resi[[1]], resi[[3]], resi[[2]], colMeans(Ms), cov(Ms))
  Zmldm[perIndex,] <- resi[[4]]
}
initParametersmLDM <- initParameters 
Z_init <- Zmldm

ZK_init <- matrix(0, K*n, p)
for(i in 1:n) {
  for(j in 1:K) {
    ZK_init[(i-1)*K + j,] <- Z_init[i,]
  }
}

# Choose the lambda1 and lambda2 for k-mLDM
ratio1 <- 0.6
ratio2 <- 0.9
modelSelectionNumber <- 1

if(abs(ratio1 - ratio2) < 1e-3) {
  modelSelectionNumber <- 1
}

lambda1_list <- seq(ratio1, ratio2, length = modelSelectionNumber)
lambda2_list <- seq(ratio1, ratio2, length = modelSelectionNumber)

# Record the minimum of EBIC
EBIC_min <- 1e+20

count <- 0

length1 <- length(lambda1_list)
length2 <- length(lambda2_list)

exist_best <- FALSE
exist_per <- rep(0,length1)
best_i <- 0

loopFlag = FALSE
max_iteration <- 500
threshold <- 1e-4
approx_num_Z <- 10
max_linesearch_Z <- 30
approx_num_B <- 10
max_linesearch_B <- 30
max_iteration_B <- 100
threshold_B <- 1e-3
delta1_threshold_B <- 1e-4
delta2_threshold_B <- 0.9
sy_threshold_B <- 1e-6
max_iteration_B_coor <- 20
threshold_B_coor <- 1e-6
verbose <- TRUE
source("./program/k-Lognormal-Dirichlet-Multinomial-lfbgs-proximal-QUIC-lambda-testZ-AGP-speed.R")
source("./program/proximal-qusi-newton-coor.R")

print("#####################################################")
print("##############      k-mLDM Begin       ##############")
print("#####################################################")

time0 <- proc.time()
#length1 <- 1
#length2 <- 1
for(j in 1:length2){
  if(j > 1){
    loopFlag = TRUE;
  }
  for(i in 1:length1){
    solution <- list()
    
    if(!exist_best || exist_per[i]){
      lambda1 <- lambda1_list[i]
      lambda2 <- lambda2_list[j]
      
      solution_EBIC <- 1e+30
      
      solution <- kmLDM(X, Mscale, K, initParameters, initPi, ZK_init, lambda1, lambda2, max_iteration, threshold, approx_num_Z, max_linesearch_Z, approx_num_B, max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, delta2_threshold_B, sy_threshold_B, max_iteration_B_coor, threshold_B_coor, threshold_lbfgs = 1e-5,loopFlag, verbose)
      print("lambda1:")
      print(solution[[7]])
      print("lambda2:")
      print(solution[[8]])
      if(length(solution)!=0){
        solution_EBIC <- solution[[5]]
      }
      
      if(solution_EBIC > 0 && solution_EBIC < EBIC_min){
        EBIC_min <- solution[[5]]
        kmLDM_result <- solution
        best_i <- i
        if(j == 1 && i == 1){
          initParameters <- solution[[1]]
          initPi <- solution[[2]]
          ZK_init <- solution[[3]]
        }
      }
      #print('memory used:')
      #print(memory.profile())
      
    }
    
    count <- count + 1
    kmLDM_result_all[[count]] <- solution
  }
  
  if(best_i!=0){
    exist_per[best_i] <- TRUE
    exist_best <- TRUE
  }
}
time1 <- proc.time()

print("################################################")
print("###########        mLDM Tuning        ##########")
print("################################################")
source("./program/mLDM-B-Theta-ScaleM.R")
source("./program/Lognormal-Dirichlet-Multinomial-lbfgs-proximal-split-q-active-set-quic-B-Theta-speed.R")

classification_final <- kmLDM_result[[9]]
print("classification_final:")
print(table(classification_final))

Zmldm_final <- matrix(0, n, p)
lambda1_final <- rep(0, K)
lambda2_final <- rep(0, K)
Pi_final <- kmLDM_result[[2]]
ZK_final <- kmLDM_result[[3]]
finalParameters <- list()
for(i in 1:K) {
  perIndex <- which(classification_final == i)
  Xi <- X[perIndex, ]
  Mi <- Mscale[perIndex,]
  res <- mLDM(X = Xi, M = Mi, model_selection_num = 1, verbose = FALSE)  
  resi <- res[[1]]
  finalParameters[[i]] <- list(resi[[1]], resi[[3]], resi[[2]], colMeans(Mi), cov(Mi))
  Zmldm_final[perIndex,] <- resi[[4]]
  lambda1_final[i] <- resi[[7]][i]
  lambda2_final[i] <- resi[[8]][i]
}

print("mLDM tuning finish:")
kmLDM_after <- list(finalParameters, Pi_final, ZK_final, Zmldm_final, lambda1_final, lambda2_final, classification_final)

# record all results for LDM model
kmLDM_record <- list(kmLDM_result, kmLDM_result_all, lambda1_list, lambda2_list)
simulation <- list(kmLDM_record, otherResults, kmLDM_after)

resultFile <- paste(subDir, 'result-', as.character(p), '-', as.character(q), '-', as.character(n), "-", as.character(K), "-", as.character(t), '-', graph , ".RData", sep="")

save(file = resultFile, simulation)

warnings()
