

# tab is the actual data set to be analysed containing the individuals 
# and all the variables coming from each each block. 
# group contains the number of variables per block 
# ndimension  is the number of common components to print

ComDim <- function(tab, group, ndimension) {
  
  ntab  <- length(group)                # Number of blocks 
  
  nind  <- nrow(tab)                    # Number of individuals
  
  p     <- ncol(tab)                    # Total number of variables 
  
  ndim  <- min(nind, p) - 1             # Number of dimensions
  
  W     <- array(data = 0, dim = c(nind, nind, ntab+1)) # Association matrix
  
  LAMBDA <- matrix(data = 0, nrow = ntab, ncol = ndim)   # matrix of saliences
  
  Q     <- matrix(data = 0, nrow = nind, ncol = ndim)   # matrix of components
  
  
  Y     <- scale(tab, center = TRUE, scale = FALSE) # to center each variable
  
  J     <- rep(1:ntab, times = group)    # it indicates the block a variable 
                                         # belongs to 
  
if(!require("multigroup")){install.packages("multigroup")}
  
  for(j in 1:ntab) { Y[ ,J==j] = normM(Y[ ,J==j]) }  # Normalize each block
  
  for(j in 1:ntab) {                                 # Scalar Product Matrix
    
    Yj=as.matrix(Y[,J==j]) 
    
    W[, ,j] <- Yj%*%t(Yj)
  } 
  
  Itot = 0                                      # Initialize the total variance 
  
  for(j in 1:ntab) {                            # Total Variance Computation
    
    Itot = Itot + sum(as.matrix(W[, ,j])^2)
    
  }    
  
  explained <- matrix(0,ndim)
  
  Res       <- NULL
  
  for( dimension in 1:ndim) { 
    
    previousfit  = 100000
    
    deltafit     = 1000000
    
    threshold    = 1e-10
    
    lambda       = matrix(1, nrow = ntab)
    
    while(deltafit > threshold){ 
      
      W[, ,ntab+1] = matrix(0,nrow = nind, ncol = nind)
      
      for( ta in 1:ntab) {
        
        W[ , ,ntab+1] = W[ , ,ntab+1]+lambda[ta]*W[ , ,ta]
        
      }
    Ws = as.matrix(W[ , ,ntab+1]) 
    
    
  Svdw=svd(Ws)                                 # Perform PCA
  
  q=as.matrix(Svdw$u[ ,1])                     # Extract the first eigen vector
  
  fit = 0  
  
    for (ta in 1:ntab) {                      # Estimation of the residuals
    
      lambda[ta]=t(q)%*%as.matrix(W[ , ,ta])%*%q
      
      aux=as.matrix(W[,,ta])-lambda[ta]*q%*%t(q)
      
      fit=fit+sum(aux^2);
    }  
  
  deltafit=previousfit-fit
  
  previousfit=fit
  
  }
  
  explained[dimension,1]  = 100*sum(lambda^2)/Itot
  
  LAMBDA[ ,dimension]     = lambda 
  
  Q[ ,dimension]          = q
  
  aux=diag(1,nind)-q%*%t(q)    # deflation procedure to update the centered 
                               # and normalized data set (Y)
  
  for (j in 1:ntab) { 
    
    Y                     =     as.matrix(Y)
    
    Y[ ,J==j]             =     aux%*%Y[,J==j]; 
    
    W[ , ,j]              =     Y[,J==j]%*%t(Y[,J==j])
    
  }

}

  Res$Q                   =     Q[ ,1:ndimension]
  
  expl                    <-    matrix(0,nrow=ndimension,ncol=2)
  
  expl[ ,1]               =     explained[1:ndimension]
  
  expl[ ,2]               =     cumsum(expl[ ,1])
  
  Res$saliences           =     LAMBDA[ ,1:ndimension]  
  
  rownames(Res$Q)         <-    rownames(tab)
  
  colnames(Res$Q)         <-    paste("Dim.",1:ndimension,sep="")
  
  colnames(Res$saliences) <-    paste("Dim.",1:ndimension,sep="")
  
  rownames(Res$saliences) <-    paste("Dataset ",1:ntab,sep="")
  
  Res$expl                <-    expl
  
  rownames(Res$expl)      <-    paste("Dim.",1:ndimension,sep="")
  
  colnames(Res$expl)      <-    c("%Total Var expl", "Cumul % total Var")
  
  
  # Computation of the compromise (overall agreement)
  
if(!require("FactoMineR")){install.packages("FactoMineR")}
  
  LambdaMoyen             <-    apply(LAMBDA,2,mean)
  
  D                       =     diag(LambdaMoyen)
  
  C                       =     Q%*%sqrt(D) # Compromise
  
  # Computation of the Escouffier RV coefficient
  
  Rv                     <-     matrix(0,nrow=ntab,ncol=1)
  
  rownames(Rv)           <-     paste("Dataset ",1:ntab,sep="")
  
  for(j in 1:ntab) { 
    
    Yj       =     as.matrix(tab[,J==j]) 
    
    resRV    <-    coeffRV(Yj,C)
  
  Rv[j,1]    <-    resRV$rv }
  
  Res$RV     <-    Rv

  
  return(Res)

}














