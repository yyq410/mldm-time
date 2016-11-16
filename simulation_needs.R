computeAUC <- function(TP, FP){
  num <- length(TP)
  AUC <- 0
  tpi_1 <- 0
  fpi_1 <- 0
  for(i in 1:num){
    tp <- TP[i]
    fp <- FP[i]
    AUC <- AUC + (2 - fp - fpi_1)*(tp - tpi_1)/2
    tpi_1 <- tp
    fpi_1 <- fp
  }
  
  return(AUC)
}

# select suitable p_value as thresholds
selectThreshold <- function(pvalue, num){
    pvalue <- unique(sort(abs(pvalue)))
    percent <- seq(0,1,length = num)
    return(quantile(pvalue,percent))
}

getSymmetry <- function(X){
  p <- ncol(X)
  for(i in 1:p){
    for(j in i:p){
      X[j,i] <- X[i,j]
    }
  }
  return(X)
}


filter_edge <- function(e, threshold){
  if(abs(e) < threshold)
    return(0)
  if(e > 0){
    return(1)
  }else if(e < 0){
    return(-1)
  }
}

# set the diag -> 0
recover_Edges1 <- function(Edges, p, threshold){
  Edges_origin <- Edges
  for(i in 1:(p-1)){
    #print(i)
    Edges[i,i] <- 0
    for(j in (i+1):p){
      #print(j)
      Edges[i,j] <- filter_edge(Edges[i,j], threshold)
      Edges[j,i] <- Edges[i,j]
      if(Edges[i,j]==0){
        Edges_origin[i,j] <- 0
        Edges_origin[j,i] <- 0
      }
    }
  }
  
  Edges[p,p] <- 0
  return(list(Edges,Edges_origin))
}

recover_Edges2 <- function(Edges, q, p, threshold){
  for(i in 1:q)
    for(j in 1:p){
      Edges[i,j] <- filter_edge(Edges[i,j], threshold)
    }
  return(Edges)
}

filter_p_value <- function(pm, a, b, cm, threshold){
  Edges <- matrix(0,a,b)
  for(i in 1:a){
    for(j in 1:b){
      if(is.na(pm[i,j]) == TRUE)
        Edges[i,j] <- 0
      else{
        if(pm[i,j] != 0){
          Edges[i,j] <- filter_edge(cm[i,j], threshold)
        }
      }
    }
  }
  return(Edges)
}

# ignore the sign of edges
ignoreNP <- function(x){
  r <- x
  r[r!=0] <- 1
  return(r)
}

# compute the true positive and false positive for Edges1
computeTPFN1 <- function(real, origin){
  p <- nrow(real)
  TP <- sum(real*origin) / sum(origin)
  FP <- (sum(real) - sum(real*origin)) / (p*p - p - sum(origin))
  
  if(is.na(TP))
    TP <- 1
  if(is.na(FP))
    FP <- 1
  
  return(c(TP,FP))
}

# compute the true positive and false positive for Edges2
computeTPFN2 <- function(real, origin){
  q <- nrow(real)
  p <- ncol(real)
  TP <- sum(real*origin) / sum(origin)
  FP <- (sum(real) - sum(real*origin)) / (p*q - sum(origin))
  
  if(is.na(TP))
    TP <- 1
  if(is.na(FP))
    FP <- 1
  
  return(c(TP,FP))
}

computeD1 <- function(corr, Sigma){
  diag(Sigma) <- 0
  p <- ncol(corr)
  diag(corr) <- 0
  return(sum(abs(corr-Sigma))/(p*(p-1)))
}
computeD1_env <- function(env_corr, B){
  env_corr <- t(env_corr)
  p <- ncol(B)
  q <- nrow(B)
  return(sum(abs(env_corr-B))/(q*p))
}
computeD2 <- function(corr, Sigma){
  d <- corr - Sigma
  return(sqrt(sum(d*d)))
}
computeD2_env <- function(env_corr, B){
  env_corr <- t(env_corr)
  d <- env_corr - B
  return(sqrt(sum(d*d)))
}

computeDirTF <- function(Dir_B, Edges2, threshold_list){
  num <- length(threshold_list)
  Dir_tf <- matrix(0,num,2)
  for(i in 1:num){
    th <- threshold_list[i]
    Dir_cut <- Dir_B
    Dir_cut[abs(Dir_cut)< th] <- 0
    Dir_cut[abs(Dir_cut)>= th] <- 1
    Edges2_dir <- t(ignoreNP(Dir_cut))
    tf_2 <- computeTPFN2(Edges2_dir, Edges2)
    Dir_tf[i,] <- tf_2
  }
  return(Dir_tf)
}

filter_Dir <- function(Dir_B, Edges2, B){
  Edges_Dir <- t(ignoreNP(Dir_B))
  tf_2 <- computeTPFN2(Edges_Dir, Edges2)
  d1_2 <- computeD1_env(Dir_B, B)
  return(c(tf_2, d1_2))
}

# compute the true positive and false positive for pcc result on different threshold
computePCCTF <- function(threshold, pcc_pvalue, pcc_env_pvalue, pcc_corr, pcc_env_corr, Edges1, Edges2){
  # record the true positive and false positive
  num <- length(threshold)
  pcc_tf <- matrix(0,num,4)
  
  for(i in 1:num){
    th <- threshold[i]
    pcc_cut <- pcc_pvalue + 1e-20
    pcc_cut[pcc_cut > th] <- 0
    pcc_env_cut <- pcc_env_pvalue + 1e-20
    pcc_env_cut[pcc_env_cut > th] <- 0
    p <- nrow(pcc_env_pvalue)
    q <- ncol(pcc_env_pvalue)
    
    Edges1_pcc <- filter_p_value(pcc_cut,p,p,pcc_corr,1e-20)
    Edges2_pcc <- filter_p_value(t(pcc_env_cut),q,p,t(pcc_env_corr),1e-20)
    diag(Edges1_pcc) <- 0
    Edges1_pcc_ignore <- ignoreNP(Edges1_pcc)
    Edges2_pcc_ignore <- ignoreNP(Edges2_pcc)
    
    tf_1 <- computeTPFN1(Edges1_pcc_ignore, Edges1)
    tf_2 <- computeTPFN2(Edges2_pcc_ignore, Edges2)
    pcc_tf[i,] <- c(tf_1, tf_2)
  }
  
  return(pcc_tf)
}



filter_corr <- function(pvalue, env_pvalue, corr, env_corr, threshold, Edges1, Edges2, Sigma, B){
  pvalue <- pvalue + 1e-20
  pvalue[pvalue > threshold] <- 0
  env_pvalue[env_pvalue > threshold] <- 0
  
  corr[pvalue == 0] <- 0
  env_corr[env_pvalue == 0] <- 0
  Edges1_real <- ignoreNP(corr)
  diag(Edges1_real) <- 0
  Edges2_real <- ignoreNP(env_corr)
  
  tf_1 <- computeTPFN1(Edges1_real, Edges1)
  tf_2 <- computeTPFN2(t(Edges2_real), Edges2)
  
  d1_1 <- computeD1(corr,Sigma)
  d1_2 <- computeD1_env(env_corr, B)
  
  return(c(tf_1,tf_2,d1_1,d1_2))
}

# compute the true positive and false positive for spearman result on different threshold
computeSpearTF <- function(threshold, spearman_pvalue, spearman_env_pvalue, spearman_corr, spearman_env_corr, Edges1, Edges2){
  spearman_tf <- computePCCTF(threshold, spearman_pvalue, spearman_env_pvalue, spearman_corr, spearman_env_corr, Edges1, Edges2)
  return(spearman_tf)
}

computeGLTF <- function(icov_origin, icov_env_origin, Edges1, Edges2, threshold){
  lengtht <- length(threshold)
  gl_tf <- matrix(0,lengtht+1,4)
  
  for(i in 1:lengtht){
    th <- threshold[lengtht - i + 1]
    icov <- icov_origin
    icov_env <- icov_env_origin
    
    icov[abs(icov) < th] <- 0
    icov[abs(icov) >= th] <- 1
    icov_env[abs(icov_env) < th] <- 0
    icov_env[abs(icov_env) >= th] <- 1
    Edges1_gl <- ignoreNP(icov)
    diag(Edges1_gl) <- 0
    Edges2_gl <- ignoreNP(icov_env)
    
    tf_1 <- computeTPFN1(Edges1_gl, Edges1)
    tf_2 <- computeTPFN2(Edges2_gl, Edges2)
    
    gl_tf[i+1,] <- c(tf_1, tf_2)
  }
  
  return(gl_tf)
}

filter_gl_cond <- function(gl_cond, gl_cond_env, gl_icov, gl_icov_env, Edges1, Edges2, Theta, B){
  tf_1 <- computeTPFN1(gl_cond, Edges1)
  tf_2 <- computeTPFN2(gl_cond_env, Edges2)
  
  d1_1 <- computeD1(gl_icov, Theta)
  d1_2 <- computeD1_env(t(gl_icov_cond_env), B)
  
  return(c(tf_1,tf_2,d1_1,d1_2))
}

# compute the TF and FN for mlasso on different lambda
computeMLTF <- function(ml, Edges1, Edges2){
  ml_tf <- computeGLTF(ml, Edges1, Edges2)
  return(ml_tf)
}

filter_ml_cond <- function(cond, cond_env, Edges1, Edges2){
  tf_1 <- computeTPFN1(cond, Edges1)
  tf_2 <- computeTPFN2(cond_env, Edges2)
  
  return(c(tf_1,tf_2))
}

# compute the TF and FN for sparcc on different lambda
computeSpccTF <- function(threshold, sparcc_cov, Edges1){
  num <- length(threshold)
  sp_tf <- matrix(0,num,2)
  
  for(i in 1:num){
    th <- threshold[num - i + 1]
    sparcc_cut <- sparcc_cov
    sparcc_cut[abs(sparcc_cut) <= th] <- 0
    sparcc_cut <- ignoreNP(sparcc_cut)
    diag(sparcc_cut) <- 0
    
    tf_1 <- computeTPFN1(sparcc_cut, Edges1)
    #print(tf_1)
    sp_tf[i,] <- tf_1
  }
  
  return(sp_tf)
}

filter_sparcc_corr <- function(sparcc_cor, threshold, Edges1, Sigma){
  sparcc_cor[abs(sparcc_cor) < threshold] <- 0
  diag(sparcc_cor) <- 0
  Edges1_sparcc <- ignoreNP(sparcc_cor)
  
  tf_1 <- computeTPFN1(Edges1_sparcc, Edges1)
  d1_1 <- computeD1(sparcc_cor, Sigma)
  
  return(c(tf_1, d1_1))
}

filter_lsa <- function(lsa_otu, lsa_meta, Edges1, Edges2, threshold){
  lsa_otu[abs(lsa_otu) < threshold] <- 0
  lsa_meta[abs(lsa_meta) < threshold] <- 0
  Edges1_lsa <- ignoreNP(lsa_otu)
  Edges2_lsa <- ignoreNP(lsa_meta)

  tf_1 <- computeTPFN1(Edges1_lsa, Edges1)
  tf_2 <- computeTPFN2(Edges2_lsa, Edges2)

  return(c(tf_1,tf_2))
}

# compute the TF and FN for ccrepe on different threshold
computeCCreTF <- function(threshold, ccrepe_p_value, ccrepe_corr, Edges1){
  num <- length(threshold)
  cc_tf <- matrix(0,num,2)
  p <- ncol(ccrepe_p_value)
  
  for(i in 1:num){
    th <- threshold[i]
    ccrepe_p_cut <- ccrepe_p_value + 1e-20
    diag(ccrepe_p_cut) <- 0
    ccrepe_p_cut[ccrepe_p_cut > th] <- 0
    
    Edges1_ccrepe <- filter_p_value(ccrepe_p_cut, p, p, ccrepe_corr, 1e-20)
    Edges1_ccrepe <- ignoreNP(Edges1_ccrepe)
    
    tf_1 <- computeTPFN1(Edges1_ccrepe, Edges1)
    cc_tf[i,] <- tf_1
  }
  
  return(cc_tf)
}

filter_ccrepe_corr <- function(corr, p_value, threshold, Edges1, Sigma){
  diag(corr) <- 1e-20
  corr[p_value>threshold] <- 0
  Edges1_ccrepe <- ignoreNP(corr)
  
  tf_1 <- computeTPFN1(Edges1_ccrepe, Edges1)
  d1_1 <- computeD1(corr, Sigma)
  
  return(c(tf_1,d1_1))
}

computeSpiecTF <- function(spiec_icov, Edges1, threshold, flag){
  num <- length(threshold)
  p <- ncol(spiec_icov)
  spiec_tf <- matrix(0,num+1,2)
  for(i in 1:num){
    th <- threshold[num - i + 1]
    icov <- spiec_icov
    icov[abs(icov) < th] <- 0
    icov[abs(icov) >= th] <- 1
    
    if(flag){
     for(a in 1:p){
       for(b in 1:p){
         if(icov[a,b]!=0){
           icov[b,a] <- 1
         }
         else if(icov[b,a]!=0){
           icov[a,b] <- 1
         }
       }
     }
    }
    
    Edges1_gl <- ignoreNP(icov)
    diag(Edges1_gl) <- 0
    
    tf_1 <- computeTPFN1(Edges1_gl, Edges1)
    
    spiec_tf[i+1,] <- tf_1
  }
  
  return(spiec_tf)
}

filter_spiec_gl <- function(icov, cond, Edges1, Theta){
  tf_1 <- computeTPFN1(cond, Edges1)
  d1_1 <- computeD1(icov, Theta) 
  
  return(c(tf_1, d1_1))
}
filter_spiec_ml <- function(cond, Edges1){
  tf_1 <- computeTPFN1(cond, Edges1)
  
  return(c(tf_1))
}

filter_cclasso <- function(corr, Edges1_cclasso, Edges1, Sigma){
  tf_1 <- computeTPFN1(Edges1_cclasso, Edges1)
  d1_1 <- computeD1(corr, Sigma)
  
  return(c(tf_1, d1_1))
}

# compute the TF and FN for LDM model on different lambda
computeLDMTF <- function(LDM_result_all, Edges1, Edges2){
  num <- length(LDM_result_all)
  LDM_tf <- matrix(0,num+1,4)
  p <- ncol(Edges1)
  q <- nrow(Edges2)
  
  for(i in 1:num){
    per <- LDM_result_all[[i]]
    if(length(per) != 0){
      B <- per[[1]]
      Theta <- per[[3]]
      
      #Edges1_LDM <- recover_Edges1(Theta, p)
      Edges1_LDM <- ignoreNP(Theta)
      diag(Edges1_LDM) <- 0
      
      #Edges2_LDM <- recover_Edges2(B, q, p)
      Edges2_LDM <- ignoreNP(B)
      
      tf_1 <- computeTPFN1(Edges1_LDM, Edges1)
      tf_2 <- computeTPFN2(Edges2_LDM, Edges2)
      LDM_tf[i+1,] <- c(tf_1, tf_2)
    }
    
  }
  
  return(LDM_tf)
}

filter_LDM <- function(LDM_B, LDM_Theta, Edges1, Edges2, Theta, B){
  Edges1_LDM <- ignoreNP(LDM_Theta)
  diag(Edges1_LDM) <- 0
  Edges2_LDM <- ignoreNP(LDM_B)
  
  tf_1 <- computeTPFN1(Edges1_LDM, Edges1)
  tf_2 <- computeTPFN2(Edges2_LDM, Edges2)
  
  d1_1 <- computeD1(LDM_Theta, Theta)
  d1_2 <- computeD1_env(t(LDM_B), B)
  
  return(c(tf_1,tf_2,d1_1,d1_2))
}

computeF1 <- function(X, X_true){
  return(max(abs((X - diag(diag(X)) - (X_true - diag(diag(X_true)))))))
}
