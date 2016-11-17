

# compute Bt for proximal qusi newton method
computeBtQ <- function(gammat, St, Yt, approx_num, approx_count, approx_index){
  d <- nrow(St)
  Q <- matrix(0,d,2)
  Qh <- matrix(0,2,d)
  BtQ <- list()
  
  if(approx_count != 0){
    real_index <- rep(0, approx_count)
    for(i in 1:approx_count){
      delta <- approx_count - i
      real_index[i] <- (approx_index - delta + approx_num - 1)%%approx_num + 1
    }
    
    Stp <- St[,real_index]
    Ytp <- Yt[,real_index]
    
    #   print("Stp")
    #    print(Stp)
    #   print("Ytp")
    #    print(Ytp)
    
    Q <- cbind(gammat*Stp, Ytp)
    # print("Qin ")
    # print(Q)
    SY <- t(Stp)%*%Ytp
    
    if(approx_count!=1){
      Dt <- diag(diag(SY))
    }else{
      Dt <- diag(SY)
    }
    Dt <- as.matrix(Dt)
    BtQ$Dt <- Dt
    
    Lt <- SY
    Lt[upper.tri(Lt, diag=TRUE)] <- 0
    #print("gammat")
    #print(gammat)
    #print("SY")
    # print(SY)
    # print("Dt")
    # print(Dt)
    Rt <- rbind(cbind(gammat*t(Stp)%*%Stp, Lt), cbind(t(Lt),-Dt))
    Rt <- solve(Rt)
    
    Qh <- Rt%*%t(Q)
    #print("Qh")
    #print(Qh)
  }
  
  BtQ$Q <- Q
  BtQ$Qh <- Qh
  
  return(BtQ)
}

softthrehold <- function(a, b){
  return(sign(a)*max(abs(a) - b, 0))
}

# set near zero value to zero~
filter_dir <- function(direction){
  threshold <- 1e-9
  direction[abs(direction) < threshold] <- 0
  return(direction)
}

scale_dir <- function(direction){
  direction <- direction / sqrt(sum(direction*direction))
  return(direction)
}

# select the not zero elements in w and gradient of w
# as the active set
getActiveSet <- function(w, gt, lambda, d){
  active <- list()
  threshold <- 1e-6
  act <- c()
  count <- 0
  
  for(i in 1:d){
    a <- w[i]
    b <- gt[i]
    add_flag <- TRUE
    if(abs(a) < threshold){
      subgt <- max(0, b - lambda, -(b + lambda))
      if(abs(subgt) < threshold){
        add_flag <- FALSE
      }
    }
    
    if(add_flag){
      count <- count + 1
      act[count] <- i
    }
  }
  
  active$act <- act
  active$num <- count
  
  return(active)
}

# get the direciton via coordinate descent
coordinateDescentD <- function(w, gammat, Q, Qh, gt, lambda, max_iteration, threshold){
  d <- nrow(Q)
  wt <- w
  wt_old <- w
  activeSet <- getActiveSet(w, gt, lambda, d)
  act <- activeSet$act
  act_size <- activeSet$num
  #act <- c(1:d)
  #act_size <- d
  #print('active set size:')
  #print(act_size)
  
  # set the number of iteration to the half of size of active set
  #max_iteration <- round(activeSet$num*0.5)
  round <- 0
  Bt <- diag(gammat*rep(1,d)) - Q%*%Qh
  delta <- 1 
  isNan <- FALSE
  
  while(delta > threshold && round < max_iteration){
    round <- round + 1
    #print('round coor:')
    #print(round)
    
    wt_old <- wt
    for(index in 1:act_size){
      i <- act[index]
      a <- Bt[i,i]
      b <- gt[i] + Bt[i,]%*%(wt-w) - a*wt[i]
      
      wt[i] <- softthrehold(-b, lambda)/a
    }
    
    delta <- sum(abs(wt - wt_old)) / d
    
    if(sum(is.infinite(delta)) || sum(is.na(delta))){
      print('coordinate error!')
      print(delta)
      isNan <- TRUE
      break
    }
    #print('delta:')
    #print(delta)
  }
  if(!isNan) {
    direction <- wt - w
  } else {
    direction <- wt_old - w 
    print("coordinate error direction:")
    print(direction)
  }

  direction <- scale_dir(direction)
  direction <- filter_dir(direction)
  # print('direction scale not zero number:')
  # print(sum(abs(direction)>1e-10))
  
  return(direction)
}


L1_norm <- function(x){
  return(sum(abs(x)))
}


# find new w via linesearch based on strong wolfe condition
linesearch <- function(objFunc, derObjFunc, w, direction, lambda, max_linesearch, f0, g0, delta1, delta2, ...){
  beta <- 0.5
  alpha <- 1
  k <- 8
  exist <- FALSE
  wt <- w
  ret <- list()
  f1 <- f0
  g1 <- g0
  if(sum(is.na(direction))!=0){
    print('linesearch direction NaN!!!!')
    print('NaN direction:')
    print(direction)
    ret$exist <- FALSE
    ret$wt <- w
    ret$value <- f0
    ret$grad <- g0
    ret$k <- 0
    return(ret)
  }
  delta1 <- delta1
  delta2 <- delta2
  d1 <- as.numeric(t(g0)%*%direction)
  
  if(is.infinite(d1) || is.na(d1)){
    print('linesearch direction d1 error!')
    print('NaN d1:')
    print(d1)
    ret$exist <- FALSE
    ret$wt <- w
    ret$value <- f0
    ret$grad <- g0
    ret$k <- 0
    return(ret)
  }
  # if(d1 > 0){
  #     print('linesearch is not decrease !!!!')
  #     print('increase d1:')
  #     print(d1)
  #     ret$exist <- FALSE
  #     ret$wt <- w
  #     ret$value <- f0
  #     ret$grad <- g0
  #     ret$k <- 0
  #     return(ret)
  # }
  dg0 <- d1 + lambda*L1_norm(w)
  #print('dg0:')
  #print(dg0)
  #print('lambda*L1_norm(w):')
  #print(lambda*L1_norm(w))
  d1 <- d1 + lambda*(L1_norm(w+direction) - L1_norm(w))
  d1 <- d1*delta1
  
  while(k <= max_linesearch){
    alpha <- beta^k
    f1 <- objFunc(xv=w+alpha*direction, ...)
    f1 <- f1 + lambda*L1_norm(w+alpha*direction)
    if(is.infinite(f1) || is.na(f1)){
      print('line search error!')
      print(f1)
      print('alpha:')
      print(alpha)
    }
    part <- alpha*d1
    
    k <- k + 1
    # print('f0:')
    # print(f0)
    # print('f1:')
    # print(f1)
    # print('part:')
    # print(part)
    if(f1 <= f0 + part){
      g1 <- derObjFunc(xv=w+alpha*direction, ...)
      dg1 <- as.numeric(t(g1)%*%direction) + alpha*lambda*L1_norm(w+direction)
      #print('dg1:')
      #print(dg1)
      #print('lambda*L1_norm(w+alpha*direction):')
      #print(alpha*lambda*L1_norm(w+direction))
      if(abs(dg1/dg0) <= delta2){
        # print('line search ok~')
        exist <- TRUE
        wt <- w+alpha*direction
        break
      }
    }
  }
  
  #wt <- filter_dir(wt)
  
  ret$wt <- wt
  ret$exist <- exist
  ret$value <- f1
  ret$grad <- g1
  ret$k <- k
  
  return(ret)
}


# estimate w via proximal quasi-newton method
proximalQusiNewtonB <- function(objFunc, derObjFunc, w, lambda, approx_num, max_linesearch, max_iteration, threshold, delta1_threshold, delta2_threshold, sy_threshold, max_iteration_coor, threshold_coor, ...){
  approx_count <- 0
  approx_index <- 0
  gammat <- 1
  
  d <- length(w)
  St <- matrix(0,d,approx_num)
  Yt <- matrix(0,d,approx_num)
  
  round <- 0
  # when SYdot < sy_threshold, change to steepest descent
  steepest_active <- FALSE
  steepest_active_height <- 1
  steepest_count <- 0
  
  objNew <- 0
  objOld <- objFunc(xv=w, ...) + lambda*L1_norm(w)
  gt0 <- derObjFunc(xv=w, ...)
  
  # print("gt0")
  #print(gt0)
  
  ret <- list()
  
  while(TRUE){
    round <- round + 1
    
    # print('Round B:')
    # print(round)
    if(!steepest_active){
      #print("St")
      #print(St)
      #print("Yt")
      #print(Yt)
      BtQ <- computeBtQ(gammat, St, Yt, approx_num, approx_count, approx_index)
      Q <- BtQ$Q
      Qh <- BtQ$Qh
      
      #print("Q")
      #print(Q)
      #print("Qh")
      #print(Qh)
      # adjust gammat and guarantee Bt positive definite
      #gammat <- checkBtPD(gammat, BtQ, round, d)
      #print('gammat adjust:')
      # print(gammat)
      
      #Bt <- diag(gammat*rep(1,d)) - Q%*%Qh
      
      direction <- coordinateDescentD(w, gammat, Q, Qh, gt0, lambda, max_iteration = max_iteration_coor, threshold = threshold_coor)
    }else{
      direction <- coordinateDescentD(w, 1, matrix(0,d,2), matrix(0,2,d), gt0, lambda, max_iteration = max_iteration_coor, threshold = threshold_coor)
      steepest_count <- steepest_count + 1
      if(steepest_count >= steepest_active_height){
        steepest_active <- FALSE
      }
    }
    #print('direciton norm1:')
    #print(sum(abs(direction)))
    
    # print("direction")
    # print(direction)
    
    line <- linesearch(objFunc=objFunc, derObjFunc=derObjFunc, w=w, direction=direction, lambda=lambda, max_linesearch=max_linesearch, f0=objOld, g0=gt0, delta1=delta1_threshold, delta2=delta2_threshold, ...)
    
    exist <- line$exist
    if(!exist){
      # print('line search failed!!!')
      break
    }
    wt <- line$wt
    gt1 <- line$grad
    objNew <- line$value
    
    # print('line search k:')
    # print(line$k)
    
    # print('obj new:')
    # print(objNew)
    
    delta <- objNew - objOld
    # print('delta:')
    # print(delta)
    
    # print("wt")
    # print(wt)
    
    # update BFGS
    St_per <- wt - w
    Yt_per <- gt1 - gt0
    SYdot <- as.numeric(t(St_per)%*%Yt_per)
    
    if(SYdot > sy_threshold){
      approx_count <- approx_count + 1
      if(approx_count > approx_num){
        approx_count <- approx_num
      }
      #  print('approx_count:')
      #  print(approx_count)
      approx_index <- approx_index %% approx_num + 1
      St[,approx_index] <- St_per
      Yt[,approx_index] <- Yt_per
      gammat <- SYdot/(t(St_per)%*%St_per)
      gammat <- as.numeric(gammat)
    }else{
      #print(paste('bi sydot <= 0 when round:', as.character(round)))
      steepest_active <- TRUE
      steepest_count <- 0
      #print('steepest descent active!')
    }
    # print('gammat:')
    # print(gammat)
    
    objOld <- objNew
    gt0 <- gt1
    w <- wt
    
    if(abs(delta) < threshold || round > max_iteration){
      # print('proximal B stop!')
      break
    }
  }
  
  ret$par <- w
  ret$value <- objOld
  ret$gradient <- gt0
  return(ret)
}
