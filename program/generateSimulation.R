getInitParameters <- function(X, M, Z_mean) {
  # set initial values for every variable: Z, B, Theta, B0
  Xr <- X / rowSums(X)
  cor_x <- cor(Xr, method = "spearman")
  cor_m <- cor(x = Xr, y = M, method = "spearman")
  corx_list <- unique(rep(abs(cor_x)))
  corm_list <- unique(rep(abs(cor_m)))
  
  # for B
  B_init <- t(cor_m)
  # for Theta
  Theta_init <- cor_x
  diag_value <- 1
  p <- ncol(X)
  n <- nrow(X)
  if(n < p) {
    print("The number of samples is less than the number of OTUs")
  }
  while(det(Theta_init) < p){
    diag(Theta_init) <- diag(Theta_init) + rep(diag_value, p)
  }
  Theta_init <- solve(Theta_init)
  # for B0
  B0_init <- colMeans(Z_init)
  
  muM <- colMeans(M)
  covM <- cov(M)
  if(det(covM) < 1e-8) {
    print("The number of samples is less than the number of meta factors")
  }
  
  result <- list()
  result[[1]] <- list(B_init, Theta_init, B0_init, muM, covM)
  result[[2]] <- list(corx_list, corm_list)
  return(result)
}

generateSimulation <- function(n, p, q, graph_number, index) {
  # the probability of B[i,j] != 0
  probB <- 0.85
  bound <- 0.5
  # the probability is a big OTU
  probB0 <- 0.2
  max_B0l <- 6
  max_B0r <- 8
  min_B0l <- 2
  min_B0r <- 4.5
  
  graph_set <- c('random','cluster','scale-free','hub','band')
  graph <- graph_set[graph_number]
  
  # generate theta with specific graph
  g <- huge.generator(n, p, graph=graph, v = 0.4)
  # the covariance matrix 
  Sigma <- g$sigma
  # the precision matrix
  Theta <- g$omega
  
  # the adjacent matrix for the real graph
  Edges <- g$theta
  
  mu <- rep(0,p)
  Z <- mvrnorm(n, mu, Sigma)
  
  # generate meta data
  #M <- mvrnorm(n, rep(0,q), diag(rep(1,q)))
  gm <- huge.generator(n, q, graph=graph, v = 0.5)
  M <- gm$data
  #muMLeft <- index * 2
  #muMRight <- (index + 0.5) * 2
  covM <- gm$sigma
  # Mscale <- scale(M, center=TRUE, scale=TRUE)
  
  # generate B
  B_whole <- runif(q*p, min=bound - 0.5, max=bound)
  
  filterZero <- function(x) {
    for (i in 1:length(x)) {
      r <- runif(1)
      
      if (r < probB) {
        x[i] <- 0
      } else {
          r2 <- runif(1)
          if(r2 < 0.5) {
              x[i] <- -x[i]
          }
      }
    }
    return(x)
  }
  
  B <- matrix(filterZero(B_whole), q, p)
  
  #generate B0
  B0 <- rep(0,p)
  
  for(i in 1:p){
    r <- runif(1)
    if(r < probB0){
      B0[i] <- runif(1,max_B0l,max_B0r)
    }else{
      B0[i] <- runif(1,min_B0l,min_B0r)
    }
  }
  
  
  # compute the alpha p*n
  part1 <- t(B)%*%t(M) + B0
  part2 <- t(Z)
  Alpha <- exp(part1 + part2)
  
  muMLeft <- index
  muMRight <- (index + 0.5) * index
  muM <- runif(q, muMLeft, muMRight)
  M <- M + muM
  
  library('HMP')
  # sampling from the dirichlet.multinomial distribution
  data <- matrix(0, n, p)
  # the minimization number of per sample
  start <- 5000
  # the maximization number of per sample
  end <- 10000
  perData <- rep(0, p)
  for (i in 1:n){
    perSum <- floor(runif(1, start, end))
    perData <- Dirichlet.multinomial(c(perSum), Alpha[,i])
    data[i,] <- perData
  }

  res <- list()
  res[[1]] <- data
  res[[2]] <- M
  res[[3]] <- list(B, Theta, B0, muM, covM)
  res[[4]] <- graph
  return(res)
}
