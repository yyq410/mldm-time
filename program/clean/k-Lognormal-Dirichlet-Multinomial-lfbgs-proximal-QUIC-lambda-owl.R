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

showClassification <- function(PiProb, n, K, round) {
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
  for(i in 1:K) {
    pn <- n / K
    print("clsuter i:")
    print(i)
    print(table(classification[((i-1)*pn+1):(i*pn)]))
  }
  print("index mean:")
  print(classiIndex / table(classification))
#      print("PiProb: ")
#      print(PiProb)
#      print("PiLog:")
#      print(PiLog)

}

# estimate the OTU-OTU associations and EF-OTU associations within K clusters
kmLDM <- function(X, M, K, initParameters, initPi, Z_init, lambda1, lambda2, max_iteration,
                  threshold, approx_num_Z, max_linesearch_Z, approx_num_B, 
                  max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, 
                  delta2_threshold_B, sy_threshold_B, max_iteration_B_coor, 
                  threshold_B_coor, threshold_lbfgs = 1e-5, loopFlag, verbose) {
  p <- ncol(X)
  q <- ncol(M)
  n <- nrow(X)
  # record the value of objective function
  objOld <- 0
  objNew <- 0
  objOldB <- 0
  objOldTheta <- 0
  objNewB <- 0
  objNewTheta <- 0

  Z <- Z_init
  # record every z for every clusters
  ZK <- matrix(0, n*K, p)
  for(i in 1:n) {
      for(j in 1:K) {
          ZK[(i-1)*K + j,] <- Z[i,]
      }
  }
  derZK <- matrix(0, n*K, p)
  
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
  initRoundB <- 200
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
  
  while(TRUE) {
    round <- round + 1
    
    # estimate the ZK
    for(i in 1:n) {
      xi <- X[i,]
      mi <- M[i,]
      
      for(j in 1:K) {
        zik <- ZK[(i-1)*K + j,]
        zik_result <- lbfgs(call_eval=objZi, call_grad=derObjZi, vars=zik, mi=mi, xi=xi, 
                         parameters=parameters[[j]], n=n, invisible=1, m=approx_num_Z, max_linesearch=max_linesearch_Z, epsilon = threshold_lbfgs)
        
        ZK[(i-1)*K + j,] <- zik_result$par
        derZK[(i-1)*K + j,] <- derObjZi(ZK[(i-1)*K + j,], mi, xi, parameters[[j]], n)
      }
      
    }
    
    # E-step
    for(i in 1:n) {
      for(j in 1:K) {
        PiProb[i,j] <- getLogProbPer(X[i,], M[i,], ZK[(i-1)*K + j,], parameters[[j]], Pi[j])
      }
    }
    PiLog <- PiProb 
    PiProbMax <- apply(PiProb, 1, max)
    PiProbSum <- PiProbMax +  log(rowSums(exp(PiProb - PiProbMax)))
    PiProb <- exp(PiProb - PiProbSum)
    # filter near 0 value
    PiProb[PiProb < 1e-10] <- 0
    PiProb <- PiProb / rowSums(PiProb)

    # calculate the lambda1 and lambda2 for every cluster
    for(i in 1:n) {
        classIndex[i] <- which(PiProb[i,] == max(PiProb[i,]))
    }
    for(i in 1:K) {
        perIndex <- which(classIndex == i)
        Xrk <- Xr[perIndex,]
        Mk <- M[perIndex,]
        
        lambda1_list[i] <- quantile(unique(rep(abs(cor(X[perIndex,], method = "spearman")))), lambda1)
        lambda2_list[i] <- quantile(unique(rep(abs(cor(x = Xrk, y = Mk, method = "spearman")))), lambda2)
    }

    # M-step

    # estimate the B0 for every cluster
    for(i in 1:K) {
      index <- i %% K
      parameters[[i]][[3]] <- colSums(ZK[c(1:(K*n)) %% K == index ,] * PiProb[,i]) / sum(PiProb[,i])
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
    
    if(loopFlag) {
      roundB <- roundB + 1 
      objOldB <- objNewB 
        
      # estimate the B for every cluster
      for(i in 1:K) {
        index <- i %% K
        Z <- ZK[c(1:(K*n)) %% K == index, ]
        ZB0 <- t(Z)
        Bi <- parameters[[i]][[1]]
        for(j in 1:q) {
          bv <- Bi[j,]
          Bp <- Bi
          Bp[j,] <- 0
          BMp <- M%*%Bp
          BZ <- t(BMp) + ZB0
            
          #b_result_per <- proximalQusiNewtonB(objFunc=objbp, derObjFunc=derObjbp, w=bv, lambda=lambda2_list[i], approx_num = approx_num_B, max_linesearch = max_linesearch_B, max_iteration = max_iteration_B, threshold = threshold_B, delta1_threshold = delta1_threshold_B, delta2_threshold = delta2_threshold_B, sy_threshold = sy_threshold_B, max_iteration_coor = max_iteration_B_coor, threshold_coor = threshold_B_coor, BZ=BZ, index=j, X=X, M=M, q=q, p=p, n=n, prob = PiProb[,i])
          b_result_per <- lbfgs(call_eval = objbp, call_grad = derObjbp, vars = bv, BZ=BZ, index=j, X=X, M=M, q=q, p=p, n=n, prob = PiProb[,i], orthantwise_c = lambda2_list[i], invisible = 1, max_iterations = 500)
          Bi[j,] <- b_result_per$par
        }
          
        parameters[[i]][[1]] <- Bi
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
            showClassification(PiProb, n, K, round)
            print('norm2 ZK:')
            print(computeL2(rep(ZK)))
            print('norm2 der ZK:')
            print(computeL2(rep(derZK)))
            print("Pi:")
            print(Pi)
            for(i in 1:K) {
              index <- i %% K
              Z <- ZK[c(1:(K*n)) %% K == index, ]
              ZB0 <- t(Z)
              Bi <- parameters[[i]][[1]]
              objBivalue <- objB(Bi, ZB0, X, M, lambda2_list[i], q, p, n, PiProb[,i])
              derBi <- derObjB(Bi, ZB0, X, M, lambda2_list[i], q, p, n, PiProb[,i])  
              #print("per B:")
              #print(i)
              #print(Bi)
              print('norm2 B:')
              print(computeL2(rep(Bi)))
              print('norm2 derB:')
              print(computeL2(derBi))
            }
          }
          roundB <- 0
      } else {
          if(0 && verbose && (round %% setRound == 1)){
            print("Round B is:")
            print(roundB)
            print("round:")
            print(round)
            print("obj New:")
            print(objNewB)
            print('norm2 ZK:')
            print(computeL2(rep(ZK)))
            print('norm2 der ZK:')
            print(computeL2(rep(derZK)))
            print("Pi:")
            print(Pi)
            for(i in 1:K) {
              index <- i %% K
              Z <- ZK[c(1:(K*n)) %% K == index, ]
              ZB0 <- t(Z)
              Bi <- parameters[[i]][[1]]
              objBivalue <- objB(Bi, ZB0, X, M, lambda2_list[i], q, p, n, PiProb[,i])
              derBi <- derObjB(Bi, ZB0, X, M, lambda2_list[i], q, p, n, PiProb[,i])  
              print('norm2 B:')
              print(computeL2(rep(Bi)))
              print('norm2 derB:')
              print(computeL2(derBi))
            }
            
          }
        next
      }

    } else {
      roundTheta <- roundTheta + 1
      objOldTheta <- objNewTheta

      # estimate the Theta for every cluster
      for(i in 1:K) {
        index <- i %% K
        Z <- ZK[c(1:(K*n)) %% K == index, ]
        B0per <- parameters[[i]][[3]]
        temp <- t(Z) - B0per
        temp2 <- t(temp) * PiProb[,i]
        S <- temp %*% temp2 / sum(PiProb[,i])
      
        quic_res <- QUIC(S, rho=lambda1_list[i], msg = 0)
        parameters[[i]][[2]] <- quic_res$X
      }

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
             # for(i in 1:K) {
             #     print("Per Theta:")
             #     print(i)
             #     print(parameters[[i]][[2]])
             # }
              print("obj New:")
              print(objNewTheta)
              print("delta:")
              print(objNewTheta - objOldTheta)
              showClassification(PiProb, n, K, round)
              print('norm2 ZK:')
              print(computeL2(rep(ZK)))
              print('norm2 der ZK:')
              print(computeL2(rep(derZK)))
              print("Pi:")
              print(Pi)
          }

          roundTheta <- 0
      } else {
          if(0 && verbose && (round %% setRound == 1)) {
              print("Round Theta is:")
              print(roundTheta)
              print("round:")
              print(round)
              print("obj New:")
              print(objNewTheta)
              print('norm2 ZK:')
              print(computeL2(rep(ZK)))
              print('norm2 der ZK:')
              print(computeL2(rep(derZK)))
              print("Pi:")
              print(Pi)
            
          }
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
          print("Theta: ")
          print(perTheta)
          print("B: ")
          print(perB)
          print("muM:")
          print(parameters[[i]][[4]])
          print("Edges Theta")
          edgesTheta <- computeEdgesTheta(perTheta)
          print(edgesTheta)
          oldTheta <- parametersOld[[i]][[2]]
          edgesOld <- computeEdgesTheta(oldTheta)
          print("edges vary number:")
          print(edgesTheta - edgesOld)
        }
        print("der ZK:")
        print("mean der ZK:")
        print(mean(abs(derZK)))
        print("max min")
        print(max(derZK))
        print(min(derZK))
        }
        #print("PiProb:")
        #print(PiProb)

        for(i in 1:K) {
            print("cluster i Prob:")
            print(i)
            print("Prob range:")
            print(range(PiProb[which(classIndex == i),i]))
            print("Prob error:")
            errorIndex <- Filter(function(x) { (x <= ((i-1)*n / K)) || (x > (i*n / K)) }, which(classIndex == i))
            print(errorIndex)
            print(PiProb[errorIndex,])
        }

        print("iteration stop~")
        break
    }
    
    parametersOld <- parameters
    objOld <- objNew
  }
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
  
  optimalSolution <- list(parameters, Pi, ZK, objNew, EBIC, objNew, lambda1_list, lambda2_list, classIndex)
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
