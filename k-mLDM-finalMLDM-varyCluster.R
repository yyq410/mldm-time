setwd("~/LDM-time/kmLDM/")
#setwd("F://R code/west english channel time/")

library(lbfgs)
library(QUIC)
library(mclust)
library(huge)

source("./program/generateSimulation.R")

testTrue <- TRUE
Toy <- FALSE
testOther <- FALSE
addRandom <- FALSE
showFlag <- FALSE

args <- commandArgs()
p <- as.numeric(args[6])
q <- as.numeric(args[7])
# the number of clusters
K <- as.numeric(args[8])
t <- as.numeric(args[9])

# generate n for evary cluster
n_list <- rep(0, K)
nLeft <- 5*p 
nRight <- 10*p

for(i in 1:K){
    n_list[i] <- floor(runif(1, nLeft, nRight))
}

dir <- 'simulation-test-new-best-varyCluster-new/'
dir.create(dir)

graph_list <- list()

if(Toy){
  #load("./toyData.RData")
  #load("./test-nonlocal-4-2-100-2-5.RData")
  load("./simulation-test-new-lambda-afterTuning/p-4-q-2-n-100-K-2-scale-free-random/graph-4-2-100-2-1-scale-free-random.RData")
  X <- graph_parameters[[1]]
  M <- graph_parameters[[2]]
  trueParameters <- graph_parameters[[9]]
  p <- ncol(X)
  q <- ncol(M)
  K <- 2
  n <- nrow(X) / K
  t <- 200
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

showTrueParameters(trueParameters)

graph <- ""
sampleNum <- ""
for(i in 1:(K-1)) {
  graph <- paste(graph, graph_list[[i]], '-', sep="")
  sampleNum <- paste(sampleNum, as.character(n_list[i]), '-', sep="")
}
graph <- paste(graph, graph_list[[K]], sep="")
sampleNum <- paste(sampleNum, as.character(n_list[K]), sep="")

subDir <- paste(dir, 'p-', as.character(p), '-q-', as.character(q), '-n-', sampleNum,'-K-', as.character(K), '-', graph, "/",sep="")
dir.create(subDir)

Mscale <- scale(M, center = TRUE, scale = TRUE)
basis <- 0
for(i in 1:K) {
  n <- n_list[i]
  print("cluster M:")
  print(i)
  print("range M:")
  print(range(Mscale[(basis+1):(basis + n),]))
  basis <- basis + n
}

otherResults <- list()
if(testOther) {
  source("./program/testOtherMethods.R")
  otherResults <- testOtherMethods(X, M, subDir, t)
}


# Cluster datasets via meta data, and initialization
BIC <- mclustBIC(Mscale, G = c(K))
model <- Mclust(Mscale, x = BIC)
print("mclust cluster parameters:")
print("Mean:")
print(model$parameters$mean)
print("Sigma:")
print(model$parameters$variance$sigma)

print("cluster true number:")
print(n_list)
classification <- model$classification
print("classification mclust:")
print(table(classification))
basis <- 0
for(i in 1:K) {
    print("clsuter i:")
    print(i)
    print(table(classification[(basis+1):(basis + n_list[i])]))
    basis <- basis + n_list[i]
}


Z_mean <- 1
Z_init <- log(X+1) + Z_mean
print('Z init max min mean:')
print(max(Z_init))
print(min(Z_init))
print(mean(Z_init))

# record parameters for every cluster
initParameters <- list()
# record the weight of every cluster
initPi <- table(classification) / nrow(X)
print("init Pi:")
print(initPi)
# Record lambda1 and lambda2
corx_list <- c()
corm_list <- c()

kmLDM_result <- list()
kmLDM_result_all <- list()
source("./program/mLDM-B-Theta-noScaleM.R")
source("./program/Lognormal-Dirichlet-Multinomial-lbfgs-proximal-split-q-active-set-quic-B-Theta.R")

print("##########################################")
print("######### mLDM initialization ############")
print("##########################################")

Zmldm <- matrix(0, sum(n_list), p)
res <- list()
basis <- 0
for(i in 1:K) {
  perIndex <- which(classification == i)
  Xs <- as.matrix(X[perIndex,])
  Ms <- as.matrix(Mscale[perIndex,])
  res[[i]] <- mLDM(X = Xs, M = Ms, model_selection_num = 1, verbose = FALSE, max_iteration = 5, threshold = 1e-2)  
  resi <- res[[i]][[1]]
  initParameters[[i]] <- list(resi[[1]], resi[[3]], resi[[2]], colMeans(Ms), cov(Ms))
  #Zmldm[(basis + 1):(basis + length(perIndex)),] <- resi[[4]]
  Zmldm[perIndex,] <- resi[[4]]
  basis <- basis + length(perIndex)
}
initParametersmLDM <- initParameters 
if(showFlag) {
print("init Parameters:")
print(initParameters)
}
Z_init <- Zmldm

ZK_init <- matrix(0, K*(sum(n_list)), p)
for(i in 1:(sum(n_list))) {
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
threshold_B <- 1e-5
delta1_threshold_B <- 1e-4
delta2_threshold_B <- 0.9
sy_threshold_B <- 1e-6
max_iteration_B_coor <- 20
threshold_B_coor <- 1e-6
verbose <- TRUE
source("./program/k-Lognormal-Dirichlet-Multinomial-lfbgs-proximal-QUIC-lambda-testZ-random.R")
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
      
      solution <- kmLDM(X, Mscale, K, initParameters, initPi, ZK_init, lambda1, lambda2, max_iteration, threshold, approx_num_Z, max_linesearch_Z, approx_num_B, max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, delta2_threshold_B, sy_threshold_B, max_iteration_B_coor, threshold_B_coor, threshold_lbfgs = 1e-5,loopFlag, verbose, n_list)
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
#print('used time for LDM:')
#print(time1 - time0)

print("cluster true number:")
print(n_list)

print("################################################")
print("###########        mLDM Tuning        ##########")
print("################################################")
source("./program/mLDM-B-Theta-ScaleM.R")
source("./program/Lognormal-Dirichlet-Multinomial-lbfgs-proximal-split-q-active-set-quic-B-Theta.R")

classification_final <- kmLDM_result[[9]]
print("classification_final:")
basis <- 0
for(i in 1:K) {
    print("cluster i:")
    print(i)
    print(table(classification_final[(basis + 1) : (basis + n_list[i])]))
    basis <- basis + n_list[i]
}
Zmldm_final <- matrix(0, sum(n_list), p)
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

print("************************************************")
print("after mLDM tuning:")
if(showFlag) {
print(finalParameters)
}
kmLDM_after <- list(finalParameters, Pi_final, ZK_final, Zmldm_final, lambda1_final, lambda2_final, classification_final)
print("************************************************")

if(testTrue) {
################################################
# test give true label 
print("#####################################################")
print("##############      TEST TRUE       #################")
print("#####################################################")

basis <- 0
for(i in 1:K) {
  n <- n_list[i]
  classification[(basis+1):(basis+n)] <- i
  basis <- basis + n_list[i]
}

Z_mean <- 1
Z_init <- log(X+1) + Z_mean
print('Z init max min mean:')
print(max(Z_init))
print(min(Z_init))
print(mean(Z_init))

# record parameters for every cluster
initParameters <- list()
# record the weight of every cluster
initPi <- table(classification) / nrow(X)
print("init Pi:")
print(initPi)

# Record lambda1 and lambda2
corx_list <- c()
corm_list <- c()

source("./program/mLDM-B-Theta-noScaleM.R")
source("./program/Lognormal-Dirichlet-Multinomial-lbfgs-proximal-split-q-active-set-quic-B-Theta.R")

Zmldm <- matrix(0, sum(n_list), p)
res <- list()
basis <- 0
for(i in 1:K) {
  perIndex <- which(classification == i)
  Xs <- as.matrix(X[perIndex,])
  Ms <- as.matrix(Mscale[perIndex,])
 # res[[i]] <- mLDM(X = Xs, M = Ms, model_selection_num = 1, max_iteration = 5, threshold <- 1e-2)  
  res[[i]] <- mLDM(X = Xs, M = Ms, model_selection_num = 1)  
  resi <- res[[i]][[1]]
  initParameters[[i]] <- list(resi[[1]], resi[[3]], resi[[2]], colMeans(Ms), cov(Ms))
  Zmldm[(basis + 1):(basis + length(perIndex)),] <- resi[[4]]
  basis <- basis + length(perIndex)
}
#initParametersmLDM <- initParameters
if(showFlag) {
print(initParameters)
}

Z_init <- Zmldm

ZK_init <- matrix(0, K*(sum(n_list)), p)
for(i in 1:(sum(n_list))) {
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

lambda1_list_true <- seq(ratio1, ratio2, length = modelSelectionNumber)
lambda2_list_true <- seq(ratio1, ratio2, length = modelSelectionNumber)

# Record the minimum of EBIC
EBIC_min <- 1e+20

#kmLDM_result <- list()
#kmLDM_result_all <- list()
#count <- 0

length1 <- length(lambda1_list_true)
length2 <- length(lambda2_list_true)

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
threshold_B <- 1e-5
delta1_threshold_B <- 1e-4
delta2_threshold_B <- 0.9
sy_threshold_B <- 1e-6
max_iteration_B_coor <- 20
threshold_B_coor <- 1e-6
verbose <- TRUE
source("./program/k-Lognormal-Dirichlet-Multinomial-lfbgs-proximal-QUIC-lambda-testZ-random.R")
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
      lambda1 <- lambda1_list_true[i]
      lambda2 <- lambda2_list_true[j]
      
      solution_EBIC <- 1e+30
      
      solution <- kmLDM(X, Mscale, K, initParameters, initPi, ZK_init, lambda1, lambda2, max_iteration, threshold, approx_num_Z, max_linesearch_Z, approx_num_B, max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, delta2_threshold_B, sy_threshold_B, max_iteration_B_coor, threshold_B_coor, threshold_lbfgs = 1e-5, loopFlag, verbose, n_list)
      print("lambda1:")
      print(solution[[7]])
      print("lambda2:")
      print(solution[[8]])
      if(length(solution)!=0){
        solution_EBIC <- solution[[5]]
      }
      
      if(solution_EBIC > 0 && solution_EBIC < EBIC_min){
        EBIC_min <- solution[[5]]
        #kmLDM_result <- solution
        best_i <- i
        if(j == 1 && i == 1){
          initParameters <- solution[[1]]
          initPi <- solution[[2]]
          ZK_init <- solution[[3]]
        }
      }
    }
  }
  
  if(best_i!=0){
    exist_per[best_i] <- TRUE
    exist_best <- TRUE
  }
}
}

# record all results for LDM model
kmLDM_record <- list(kmLDM_result, kmLDM_result_all, lambda1_list, lambda2_list)
if(showFlag) {
print("#############################################")
print("##########     TRUE PARAMETERS     ##########")
print("#############################################")

for(i in 1:K) {
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
}
simulation <- list(kmLDM_record, otherResults, kmLDM_after)

resultFile <- paste(subDir, 'result-', as.character(p), '-', as.character(q), '-', sampleNum, "-", as.character(K), "-", as.character(t), '-', graph , ".RData", sep="")
graphFile <- paste(subDir, 'graph-', as.character(p), '-', as.character(q), '-', sampleNum, '-', as.character(K), "-", as.character(t), '-', graph , ".RData", sep="")

graph_parameters <- list(X, M, p, q, n_list, K, initParametersmLDM, Zmldm, trueParameters)

save(file = resultFile, simulation, graph_parameters, graph_list)
save(file = graphFile, graph_parameters, graph_list)

warnings()
