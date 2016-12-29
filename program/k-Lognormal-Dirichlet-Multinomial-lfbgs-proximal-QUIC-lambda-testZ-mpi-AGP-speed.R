# compute the det via cholesky decomposition
computeLogDet <- function(X){
  eigen_value <- eigen(X, symmetric=TRUE, only.values=TRUE)$values
  
  return(sum(log(eigen_value)))
}

# compute the log probability for every point
getLogProbPer <- function(x, m, z, init, piProb) {
  B <- init[[1]]
  B0 <- init[[3]]
  Theta <- init[[2]]
  muM <- init[[4]]
  covM <- init[[5]]
  alpha <- exp(t(B)%*%m + z)
  part1 <- sum(lgamma(alpha + x) - lgamma(alpha)) + (lgamma(sum(alpha)) - lgamma(sum(alpha + x)))
  part2 <- 0.5*(computeLogDet(Theta) - computeLogDet(covM))
  part3 <- log(piProb)
  part4 <- - 0.5*(t(z - B0)%*%Theta%*%(z - B0) + t(m - muM)%*%solve(covM)%*%(m - muM))
  
  return(part1 + part2 + part3 + part4)
}

computeL2 <- function(x){
  return(sqrt(sum(x^2)))
}

objZi <- function(zi, mi, xi, parameters, n) {
  B <- parameters[[1]]
  Theta <- parameters[[2]]
  B0 <- parameters[[3]]
  alphai <- exp(t(B) %*% mi + zi)
  alphaix <- alphai + xi
  
  part1 <- sum(lgamma(alphaix) - lgamma(alphai)) + (lgamma(sum(alphai)) - lgamma(sum(alphaix)))
  part2 <- t(zi - B0) %*% Theta %*% (zi - B0) / 2

  obj <- - part1 + part2
  #print("zi obj:")
  #print(obj / n)
  
 # if(sum(is.na(lgamma(alphai))) || sum(is.infinite(lgamma(alphai)))) {
 #   print("lgamma NaN zi:")
 #   print(zi)
 # }

  return(as.numeric(obj / n))
}

derObjZi <- function(zi, mi, xi, parameters, n) {
  B <- parameters[[1]]
  Theta <- parameters[[2]]
  B0 <- parameters[[3]]
  alphai <- exp(t(B) %*% mi + zi)
  alphaix <- alphai + xi

  part1 <- (digamma(alphaix) - digamma(alphai)) + (digamma(sum(alphai)) - digamma(sum(alphaix)))
  part2 <- Theta %*% (zi - B0)

  der <- - part1 * alphai + part2 
  
  #if(sum(is.na(digamma(alphai))) || sum(is.infinite(digamma(alphai)))) {
  #  print("digamma NaN zi:")
  #  print(zi)
  #}

  #print("der zi:")
  #print(mean(abs(obj / n)))
  return(as.numeric(der) / n)
}

# compute lgamma matrix part of objective funciton for matrix B
computeAlphaMatrix <- function(AlphaMX, AlphaM, prob){
  obj <- -sum(t(lgamma(AlphaMX) - lgamma(AlphaM))*prob) + sum((lgamma(colSums(AlphaMX)) - lgamma(colSums(AlphaM)))*prob)
  return(obj)
}

# compute digamma matrix part of objective funciton for the matrix B
computeAlphaDerMatrix <- function(AlphaMX, AlphaM, prob){
  der <- -t(digamma(AlphaMX) - digamma(AlphaM)) + (digamma(colSums(AlphaMX)) - digamma(colSums(AlphaM)))
  der <- t(der*prob)
  return(der)
}

# compute objective funciton for index-th row of the matrix B
objbp <- function(xv, index, BZ, X, M, q, p, n, prob){
  bv <- xv
  AlphaM <- exp(BZ + matrix(bv,p,1)%*%matrix(M[,index],1,n))
  AlphaMX <- AlphaM + t(X)
  obj <- computeAlphaMatrix(AlphaMX, AlphaM, prob)

  return(obj/n)
}

# compute derivate function for index-th row of the matrix B
derObjbp <- function(xv, index, BZ, X, M, q, p, n, prob){
  bv <- xv
  AlphaM <- exp(BZ + matrix(bv,p,1)%*%matrix(M[,index],1,n))
  AlphaMX <- AlphaM + t(X)
  der <- (computeAlphaDerMatrix(AlphaMX, AlphaM, prob)*AlphaM)%*%M[,index]
  der <- der / n
  
  #print("derBi")
  #print(der)
  
  return(as.vector(der))
}
# compute objective funciton for index-th row of the matrix B
objbp_owl <- function(xv, index, BZ, X, M, q, p, n, prob){
  bv <- xv
  AlphaM <- exp(BZ + matrix(bv,p,1)%*%matrix(M[,index],1,n))
  AlphaMX <- AlphaM + t(X)
  obj <- computeAlphaMatrix(AlphaMX, AlphaM, prob)

  return(obj)
}

# compute derivate function for index-th row of the matrix B
derObjbp_owl <- function(xv, index, BZ, X, M, q, p, n, prob){
  bv <- xv
  AlphaM <- exp(BZ + matrix(bv,p,1)%*%matrix(M[,index],1,n))
  AlphaMX <- AlphaM + t(X)
  der <- (computeAlphaDerMatrix(AlphaMX, AlphaM, prob)*AlphaM)%*%M[,index]
  
  #print("derBi")
  #print(der)
  
  return(as.vector(der))
}

# compute  the objective function for the matrix B
objB <- function(B, ZB0, X, M, lambda2, q, p, n, prob){
  AlphaM <- exp(t(M%*%B) + ZB0)
  AlphaMX <- AlphaM + t(X)
  
  obj <- computeAlphaMatrix(AlphaMX, AlphaM, prob)
  
  return(obj/n)
}

# compute the derivate of the objective function for the  matrix B
derObjB <- function(B, ZB0, X, M, lambda2, q, p, n, prob){
  AlphaM <- exp(t(M%*%B) + ZB0)
  AlphaMX <- AlphaM + t(X)
  der <- (computeAlphaDerMatrix(AlphaMX, AlphaM, prob)*AlphaM)%*%M
  der <- t(der)/n
  
  return(rep(der))
}

showClassification <- function(PiProb, n, K, round, nList) {
  print("Round:")
  print(round)
  classification <- rep(0, n)
  classiIndex <- rep(0, K)
  for(i in 1:n) {
    perIndex <- which(PiProb[i,] == max(PiProb[i,]))
    classification[i] <- perIndex
    classiIndex[perIndex] <- classiIndex[perIndex] + i
  }
  
  print("classification result:")
  print(table(classification))
  basis <- 0
  for(i in 1:K) {
    print("clsuter i:")
    print(i)
    print(table(classification[(basis + 1) : (basis + nList[i])]))
    basis <- basis + nList[i]
  }

}

getPiProb <- function(PiProb) {
    PiLog <- PiProb 
    PiProbMax <- apply(PiProb, 1, max)
    PiProbSum <- PiProbMax +  log(rowSums(exp(PiProb - PiProbMax)))
    PiProb <- exp(PiProb - PiProbSum)
    # filter near 0 value
    PiProb[PiProb < 1e-10] <- 0
    PiProb <- PiProb / rowSums(PiProb)

    return(list(PiLog, PiProbSum, PiProb))
}

# estimate the OTU-OTU associations and EF-OTU associations within K clusters
kmLDM <- function(X, M, K, initParameters, initPi, Z_init, lambda1, lambda2, max_iteration,
                  threshold, approx_num_Z, max_linesearch_Z, approx_num_B, 
                  max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, 
                  delta2_threshold_B, sy_threshold_B, max_iteration_B_coor, 
                  threshold_B_coor, threshold_lbfgs = 1e-5, loopFlag, verbose) {
  p <- ncol(X)
  minCluster <- 0.5*p +  q
  q <- ncol(M)
  n <- nrow(X)
  # record the value of objective function
  objOld <- 0
  objNew <- 0
  objOldB <- 0
  objOldTheta <- 0
  objNewB <- 0
  objNewTheta <- 0

  # test if the exp(Z) is infinite or NaN
  isZ_abnormal <- FALSE
  maxZ <- 700

  ZK <- Z_init
  
  parameters <- initParameters
  parametersOld <- parameters
  Pi <- initPi
  derB <- matrix(0, q, p)
  # record the round of the iterations
  round <- 0
  roundB <- 0
  setRound <- 20
  roundExchange <- 0
  initRoundTheta <- 50
  initRoundB <- 100
  roundTheta <- 0
  roundB <- 0

  # Rik
  PiProb <- matrix(0, n, K)
  PiProbSum <- rep(0, n)
  PiNeed <- matrix(0, n, 4)

  # record the best lambda
  lambda1_list <- rep(0, K)
  lambda2_list <- rep(0, K)
  Xr <- X / rowSums(X)
  classIndex <- rep(0, n)
  
  # init MPI
  maxCores <- 10
  cl <- startMPIcluster(count = min(K, maxCores))
  registerDoMPI(cl)

  while(TRUE) {
    round <- round + 1
    print(paste("round: ", as.character(round)))


    timeZ0 <- proc.time()

    # estimate the ZK
    ZK <- foreach(j = 1:K, .combine = rbind, .packages=c("lbfgs"), .export=c("objZi", "derObjZi")) %dopar% {
      Zkp <- matrix(0, n, p)
      for(i in 1:n) {
        xi <- X[i,]
        mi <- M[i,]
        zik <- ZK[(j-1)*n + i,]
        zik_result <- lbfgs(call_eval=objZi, call_grad=derObjZi, vars=zik, mi=mi, xi=xi, 
                         parameters=parameters[[j]], n=n, invisible=1, m=approx_num_Z, max_linesearch=max_linesearch_Z, epsilon = threshold_lbfgs, max_iterations = 500)
        
        Zkp[i,] <- zik_result$par
      }

      Zkp
    }
    timeZ1 <- proc.time()
    print("ZK time:")
    print(timeZ1 - timeZ0)

    # test exp(Zk) -> infinite or exp(ZK) == 0
    if(sum(abs(ZK) >= maxZ) || sum(is.na(ZK)) || sum(is.infinite(ZK))){
      isZ_abnormal <- TRUE
      print("iteration stop!")
      print("values of Z is abnormal!!!!!!!!!!") 
      print("return the latest result~")

      break
    }
    
    timeE0 <- proc.time()
    # E-step
    for(j in 1:K) {
      timeA0 <- proc.time()
      B <- parameters[[j]][[1]]
      B0 <- parameters[[j]][[3]]
      Theta <- parameters[[j]][[2]]
      muM <- parameters[[j]][[4]]
      covM <- parameters[[j]][[5]]
      piProb <- Pi[j]
      timeA1 <- proc.time()
      print("Assigan time:")
      print(timeA1 - timeA0)

      timeC0 <- proc.time()
      for(i in 1:n) {
        z <- ZK[(j-1)*n + i,]
        m <- M[i,]
        x <- X[i,]
        alpha <- exp(t(B)%*%m + z)
        part1 <- sum(lgamma(alpha + x) - lgamma(alpha)) + (lgamma(sum(alpha)) - lgamma(sum(alpha + x)))
        part2 <- 0.5*(computeLogDet(Theta) - computeLogDet(covM))
        part3 <- log(piProb)
        part4 <- - 0.5*(t(z - B0)%*%Theta%*%(z - B0) + t(m - muM)%*%solve(covM)%*%(m - muM))
        
        PiProb[i,j] <- part1 + part2 + part3 + part4
      }
      timeC1 <- proc.time()
      print("time ComPi:")
      print(timeC1 - timeC0)
    }
    
    PiLog <- PiProb 
    PiProbMax <- apply(PiProb, 1, max)
    PiProbSum <- PiProbMax +  log(rowSums(exp(PiProb - PiProbMax)))
    PiProb <- exp(PiProb - PiProbSum)
    # filter near 0 value
    PiProb[PiProb < 1e-10] <- 0
    PiProb <- PiProb / rowSums(PiProb)
    
    timeE1 <- proc.time()
    print("E time:")
    print(timeE1 - timeE0)

    timeS0 <- proc.time()
    # calculate the lambda1 and lambda2 for every cluster
    for(i in 1:n) {
        classIndex[i] <- which(PiProb[i,] == max(PiProb[i,]))
    }

    # test not zero ratio on every cluster 
    if(FALSE) {
    for(i in 1:K) {
        print(paste("cluster ", as.character(i)))
        perIndex <- which(classIndex == i)
        Xk <- X[perIndex,]
        testNotZero(Xk)
    }
    }

    # remove the samll cluster
    removeList <- c()
    needRemove <- FALSE
    removeCluster <- c()
    for(i in 1:K) {
        perIndex <- which(classIndex == i)
        perNum <- length(perIndex)
        if(perNum < minCluster) {
            removeList <- c(removeList, perIndex)
            removeCluster <- c(removeCluster, i)
            needRemove <- TRUE
        }
    }

    if(needRemove) {
        print("some small cluster need Remove!")
        print(removeCluster)
        print("removed clusters' size:")
        print(length(removeList))
        X <- X[-removeList,]
        Xr <- Xr[-removeList,]
        M <- M[-removeList,]
        parameters <- parameters[-removeCluster]
        PiLog <- PiLog[-removeList, -removeCluster]
        removeNum <- length(removeCluster)
        ZK <- ZK[-rep(vapply(removeCluster, FUN = function(x) {return(as.integer(((x-1)*n + c(1:n))))}, FUN.VALUE = integer(removeNum*n))),]
        K <- K - removeNum
        ZK <- ZK[-rep(vapply(removeList, FUN = function(x) {return(as.integer((c(0:(K-1))*n + x)))}, FUN.VALUE = integer(K))),]
        
        PiProb <- PiLog 
        PiProbMax <- apply(PiProb, 1, max)
        PiProbSum <- PiProbMax +  log(rowSums(exp(PiProb - PiProbMax)))
        PiProb <- exp(PiProb - PiProbSum)
        # filter near 0 value
        PiProb[PiProb < 1e-10] <- 0
        PiProb <- PiProb / rowSums(PiProb)
        
        n <- nrow(X)
        lambda1_list <- rep(0, K)
        lambda2_list <- rep(0, K)
        classIndex <- rep(0, n)
        # calculate the lambda1 and lambda2 for every cluster
        for(i in 1:n) {
            classIndex[i] <- which(PiProb[i,] == max(PiProb[i,]))
        }
    }
    timeS1 <- proc.time()
    print("time remove Cluster:")
    print(timeS1 - timeS0)

    timeL0 <- proc.time()
    for(i in 1:K) {
        perIndex <- which(classIndex == i)
        Xrk <- Xr[perIndex,]
        Mk <- M[perIndex,]
        corx <- cor(X[perIndex,], method = "spearman")
        corx[is.na(corx)] <- 0
        lambda1_list[i] <- quantile(unique(rep(abs(corx))), lambda1)
        corxm <- cor(x = Xrk, y = Mk, method = "spearman")
        corxm[is.na(corxm)] <- 0
        lambda2_list[i] <- quantile(unique(rep(abs(corxm))), lambda2)
    }
    timeL1 <- proc.time()
    print("time lambda:")
    print(timeL1 - timeL0)
   
    #M-step
    timeMO0 <- proc.time()
    # estimate the B0 for every cluster
    for(i in 1:K) {
      basis <- (i-1)*n
      parameters[[i]][[3]] <- colSums(ZK[basis + c(1:n),] * PiProb[,i]) / sum(PiProb[,i])
    }
    
    # estimate the muM for every cluster
    for(i in 1:K) {
      parameters[[i]][[4]] <- colSums(M * PiProb[,i]) / sum(PiProb[,i])
    }
    
    # estimate the covM for every cluster
    for(i in 1:K) {
      temp <- t(M) - parameters[[i]][[4]]
      temp2 <- t(temp) * PiProb[,i]
      parameters[[i]][[5]] <- temp %*% temp2 / sum(PiProb[,i])
    }
    
    # estiimate the Pi for every cluster
    Pi <- colSums(PiProb) / sum(PiProb)
    timeMO1 <- proc.time()
    print("time M step other:")
    print(timeMO1 - timeMO0)
    
    if(loopFlag) {
      print("in B loop:")
      roundB <- roundB + 1 
      objOldB <- objNewB 
        
      # estimate the B for every cluster
      timeB0 <- proc.time()
      B_list <- foreach(i = 1:K, .export=c("proximalQusiNewtonB", "objbp", "derObjbp", "computeAlphaMatrix", "computeAlphaDerMatrix", "getCond", "computeBtQ", "softthrehold", "filter_dir", "scale_dir", "getActiveSet", "L2_norm", "coordinateDescentD", "L1_norm", "linesearch")) %dopar% {
        basis <- (i-1)*n
        Z <- ZK[basis + c(1:n), ]
        ZB0 <- t(Z)
        Bi <- parameters[[i]][[1]]
        for(j in 1:q) {
          bv <- Bi[j,]
          Bp <- Bi
          Bp[j,] <- 0
          BMp <- M%*%Bp
          BZ <- t(BMp) + ZB0
            
          b_result_per <- proximalQusiNewtonB(objFunc=objbp, derObjFunc=derObjbp, w=bv, lambda=lambda2_list[i], approx_num = approx_num_B, max_linesearch = max_linesearch_B, max_iteration = max_iteration_B, threshold = threshold_B, delta1_threshold = delta1_threshold_B, delta2_threshold = delta2_threshold_B, sy_threshold = sy_threshold_B, max_iteration_coor = max_iteration_B_coor, threshold_coor = threshold_B_coor, BZ=BZ, index=j, X=X, M=M, q=q, p=p, n=n, prob = PiProb[,i])
          Bi[j,] <- b_result_per$par
        }
          
        Bi
      }
      timeB1 <- proc.time()
      print("B time:")
      print(timeB1 - timeB0)


      for(i in 1:K) {
        parameters[[i]][[1]] <- B_list[[i]]
      }

      objNewB <- computeObjf(PiProbSum, n, K, parameters, lambda1_list, lambda2_list)
      deltaB <- abs(objNewB - objOldB)

      if(deltaB < threshold || round >= max_iteration || (roundExchange < setRound && roundB > initRoundB)) {
          loopFlag <- FALSE
          roundExchange <- roundExchange + 1
          objNew <- objNewB
          if(verbose) {
            print("Round B is:")
            print(roundB)
            print("round:")
            print(round)
            print("obj New:")
            print(objNewB) 
            print("delta B:")
            print(objNewB - objOldB)
            print('norm2 ZK:')
            print(computeL2(rep(ZK)))
            print("max, min, mean ZK:")
            print(c(max(ZK), min(ZK), mean(abs(ZK))))
            for(i in 1:K) {
              basis <- (i-1)*n
              Z <- ZK[basis + c(1:n), ]
              ZB0 <- t(Z)
              Bi <- parameters[[i]][[1]]
              derBi <- derObjB(Bi, ZB0, X, M, lambda2_list[i], q, p, n, PiProb[,i])  
              print('norm2 B:')
              print(computeL2(rep(Bi)))
              print('norm2 derB:')
              print(computeL2(derBi))
            }
            print("classification:")
            print(table(classIndex))
            warnings()
          }
          roundB <- 0
      } else {
        next
      }

    } else {
      print("in Theta loop:")
      roundTheta <- roundTheta + 1
      objOldTheta <- objNewTheta
      
      timeT0 <- proc.time()
      # estimate the Theta for every cluster
      for(i in 1:K) {
        basis <- (i-1)*n
        Z <- ZK[basis + c(1:n), ]
        B0per <- parameters[[i]][[3]]
        temp <- t(Z) - B0per
        temp2 <- t(temp) * PiProb[,i]
        S <- temp %*% temp2 / sum(PiProb[,i])
        timeQ0 <- proc.time() 
        quic_res <- QUIC(S, rho=lambda1_list[i], msg = 0, maxIter = 500, tol = 1e-3)
        timeQ1 <- proc.time()
        print("QUIC time:")
        print(timeQ1 - timeQ0)

        parameters[[i]][[2]] <- quic_res$X
      }
      timeT1 <- proc.time()
      print("Theta time:")
      print(timeT1 - timeT0)

      objNewTheta <- computeObjf(PiProbSum, n, K, parameters, lambda1_list, lambda2_list)
      deltaTheta <- abs(objNewTheta - objOldTheta)
      if(deltaTheta < threshold || round >= max_iteration || (roundExchange < setRound && roundTheta >= initRoundTheta)) {
          loopFlag <- TRUE
          roundExchange <- roundExchange + 1
          objNew <- objNewTheta
          if(verbose) {
              print("Round Theta is:")
              print(roundTheta)
              print("round:")
              print(round)
              print("obj New:")
              print(objNewTheta)
              print("delta:")
              print(objNewTheta - objOldTheta)
              print('norm2 ZK:')
              print(computeL2(rep(ZK)))
              print("max, min, mean ZK:")
              print(c(max(ZK), min(ZK), mean(abs(ZK))))
              print("classification")
              print(table(classIndex))
              warnings()
          }

          roundTheta <- 0
      } else {
        next
      }
    }
    
    # compute the value of objective function
    #objNew <- computeObjf(PiProbSum, n, K, parameters, lambda1, lambda2)
    if(loopFlag) {
        objOld <- objNewB
    } else {
        objOld <- objNewTheta
    }
    
    # test if ternimate
    delta <- abs(objNew - objOld)
    if(delta < threshold || round > max_iteration) {
        if(verbose) {
        print("All iterations:")
        print(round)
        print("obj New:")
        print(objNew)
        print("delta : ")
        print(objNew - objOld)
        print("Pi:")
        print(Pi)
        for(i in 1:K) {
          perTheta <- parameters[[i]][[2]]
          perB <- parameters[[i]][[1]]
          print("Cluster i = ")
          print(i)
          edgesTheta <- computeEdgesTheta(perTheta)
          print(edgesTheta)
          oldTheta <- parametersOld[[i]][[2]]
          edgesOld <- computeEdgesTheta(oldTheta)
          print("edges vary number:")
          print(edgesTheta - edgesOld)
        }
        }

        print("iteration stop~")
        break
    }
    
    parametersOld <- parameters
    objOld <- objNew
  }

  closeCluster(cl)

  g <- 0.5
  # PiProbSum don't have the value of combination and constant
  EBIC <- - 2*sum(PiProbSum)
  E1 <- 0
  E2 <- 0
  for(i in 1:K) {
    E1 <- E1 + (sum(parameters[[i]][[2]]!=0) - p)/2
    E2 <- E2 + sum(parameters[[i]][[1]]!=0)
  }
  EBIC <- EBIC + (E1 + E2)*log(n) + 4*g*E1*log(p) + 2*g*E2*log(p*q)
    
  print("EBIC:")
  print(EBIC)
  
  optimalSolution <- list(parameters, Pi, ZK, objNew, EBIC, objNew, lambda1_list, lambda2_list, classIndex, isZ_abnormal)
  return(optimalSolution)
}

computeEdgesTheta <- function(Theta) {
  diag(Theta) <- 0
  return(sum(Theta != 0) / 2)
}

computeObjf <- function(PiProbSum, n, K, parameters, lambda1, lambda2) {
    objNew <- - sum(PiProbSum) / n
    for(i in 1:K) {
      objNew <- objNew + lambda1[i]*sum(abs(parameters[[i]][[2]]))/2 + lambda2[i]*sum(abs(parameters[[i]][[1]]))
    }

    return(objNew)
}

testNotZero <- function(x) {
    not0 <- colSums(x != 0) / nrow(x)
    p <- ncol(x)
    sampleSize <- rowSums(x)
    sortRes <- sort(not0)
    sortRes2 <- sort(sampleSize)
    print("min 20, top 5 otu couts:")
    print(sortRes2[1:20])
}
