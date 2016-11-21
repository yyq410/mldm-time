setwd("~/LDM-time/kmLDM/")

source("./simulation_needs.R")
if(FALSE) {
args <- commandArgs()
p <- as.numeric(args[6])
q <- as.numeric(args[7])
n <- as.numeric(args[8])
K <- as.numeric(args[9])
t <- as.numeric(args[10]) 
}

load("./simulation-test-new-lambda-afterTuning-scaleMTuning-shortB-filter-exchange/p-50-q-5-n-400-K-2-cluster-cluster/graph-50-5-400-2-3-cluster-cluster.RData")
load("./simulation-test-new-lambda-afterTuning-scaleMTuning-shortB-filter-exchange/p-50-q-5-n-400-K-2-cluster-cluster/result-50-5-400-2-3-cluster-cluster.RData")

trueParameters <- graph_parameters[[9]]
kmLDM_after <- simulation[[3]]
finalParameters <- kmLDM_after[[1]]
classificationFinal <- kmLDM_after[[7]]

p <- graph_parameters[[3]]
q <- graph_parameters[[4]]
n <- graph_parameters[[5]]
K <- graph_parameters[[6]]

X <- graph_parameters[[1]]
M <- graph_parameters[[2]]

filterEdges <- function(x) {
  x[abs(x) < 1e-3] <- 0
  return(sign(x))
}

ignoreEdges <- function(x) {
  x[abs(x) < 1e-3] <- 0
  x[x != 0] <- 1
  return(x)
}

ThetaGraphTrue <- list()
ThetaGraphTrueIgnore <- list()
BGraphTrue <- list()
BGraphTrueIgnore <- list()
ThetaGraphLearn <- list()
ThetaGraphLearnIgnore <- list()
BGraphLearn <- list()
BGraphLearnIgnore <- list()

test_num <- 20
kmldm_tf <- matrix(0, test_num + 1, K*4)
kmldm_auc <- matrix(0, K, 2) 
kmldm_classPrecision <- matrix(0, K, 2)

for(i in 1:K) { 
  perTrue <- trueParameters[[i]]
  perLearn <- finalParameters[[i]]
  ThetaGraphTrue[[i]] <- filterEdges(perTrue[[2]])
  ThetaGraphTrueIgnore[[i]] <- ignoreEdges(perTrue[[2]])
  ThetaGraphLearn[[i]] <- filterEdges(perLearn[[2]])
  ThetaGraphLearnIgnore[[i]] <- ignoreEdges(perLearn[[2]])
  BGraphTrue[[i]] <- filterEdges(perTrue[[1]])
  BGraphTrueIgnore[[i]] <- ignoreEdges(perTrue[[1]])
  BGraphLearn[[i]] <- filterEdges(perLearn[[1]])
  BGraphLearnIgnore[[i]] <- ignoreEdges(perLearn[[1]])

  Theta_true <- perTrue[[2]]
  Sigma_true <- solve(Theta_true)
  B_true <- perTrue[[1]]
# real conditional independent graph
  Edges <- ignoreEdges(Theta_true)

  Edges1_Theta_list <- recover_Edges1(Theta_true, p, 1e-8)
  Edges1_Theta_true <- Edges1_Theta_list[[1]]
  Theta_filter <- Edges1_Theta_list[[2]]
    
  Edges1_Sigma_list <- recover_Edges1(Sigma_true, p, 1e-1)
  Edges1_Sigma_true <- Edges1_Sigma_list[[1]]
  Sigma_filter <- Edges1_Sigma_list[[2]] 
    
  # if Theta[i,j] != 0 Edges1_ignore[i,j] <- 1
  Edges1_Theta_ignore <- Edges
  Edges1_Sigma_ignore <- Edges1_Sigma_true
  Edges1_Sigma_ignore[Edges1_Sigma_ignore != 0] <- 1
  # recover edges set according to the B
  # if B[i,j] > 0 edges2[i,j] = 1; 
  # if B[i,j] < 0 edges2[i,j] = -1;
  # if B[i,j] = 0 edges2[i,j] = 0
  Edges2_true <- recover_Edges2(B_true, q, p, 1e-3) 
  Edges2_ignore <- Edges2_true
  Edges2_ignore[Edges2_ignore != 0] <- 1

  # k-mldm
  kmldm_Theta <- perLearn[[2]]
  kmldm_B <- perLearn[[1]]

  LDM_th <- kmldm_Theta
  diag(LDM_th) <- 0
  threshold_LDM <- c(rep(LDM_th),rep(kmldm_B))
  threshold_list <- selectThreshold(threshold_LDM, test_num)
    
  LDM_tf <- computeGLTF(kmldm_Theta, kmldm_B, Edges1_Theta_ignore, Edges2_ignore, threshold_list) 
  #print(LDM_tf) 
  kmldm_auc[i,] <- c(computeAUC(LDM_tf[,1], LDM_tf[,2]), computeAUC(LDM_tf[,3],LDM_tf[,4]))
  kmldm_tf[,((i-1)*4+1) : (i*4)] <- LDM_tf 

  perIndex <- which(classificationFinal == i)
  allNum <- length(perIndex)
  errorNum <- sum(perIndex > (i*n)) + sum(perIndex <= ((i-1)*n))
  kmldm_classPrecision[i,] <- c((allNum - errorNum) / allNum, errorNum / allNum)
}
print("kmldm TF list:")
print(kmldm_tf)
print("kmldm AUC scores:")
print(kmldm_auc) 
print("kmldm Classification precision:") 
print(kmldm_classPrecision)
