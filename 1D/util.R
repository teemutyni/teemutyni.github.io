
##############################################################################
# This file (and funcs.R) define all the utility functions to be used in main code (main.R)
#
# T. Tyni, M. Harju, Univ. Oulu 2016
##############################################################################
source("funcs.R")
library(parallel)

# Weights and abscissas of the 2n+1 point composite Simpson rule over general interval [a,b]
simpson <- function(a,b,n){
  h = (b-a)/(2*n)
  x = seq(a,b,by=h)
  w = h*c(1,rep(c(4,2),n-1),4,1)/3
  return(matrix(c(x,w),2,byrow=TRUE))
}

# Reflection coefficient b(k) using asymptotic formula
bfun<-function(k,M1){
  f<-function(k) exp(-1i*k*1e6)*( u(-1e6,k,M1) - u_0(-1e6,k) )
  return(vapply(k,f,FUN.VALUE=complex(real=0,imaginary=0)))
}
bfun2<-function(k,M1){ # This is parallel version of bfun() above
    f<-function(k) exp(-1i*k*1e6)*( u(-1e6,k,M1) - u_0(-1e6,k) )
  
    cat("Making cluster...")
    cl <- makeCluster(detectCores() - 1)
    clusterEvalQ(cl, source("funcs.R")) 
    clusterExport(cl=cl, varlist=c("M1","supp"))
    cat("Done.\n")
    res = unlist(parLapply(cl,k,f))
    stopCluster(cl)
  
    return(res)
}

# The right hand side g(k) of the linear system Af=g
g  <- function(k,M1) 2*pi*k^3*bfun(k/2,M1)
g2 <- function(k,M1) 2*pi*k^3*bfun2(k/2,M1) # This is parallel version

# Hat functions
# x is variable, j the index of hat function, [a,b] interval of interest and N the number of hat functions
hat_j<-function(x,j,a,b,N){
  xj<-function(j) a+(j-1)*(b-a)/(N-1)
  fx <-function(x,j){
    ifelse(x >= xj(j) & x <xj(j+1), (x-xj(j+1))/(xj(j)-xj(j+1)),
           ifelse(x >= xj(j-1) & x < xj(j),  (x-xj(j-1))/(xj(j)-xj(j-1)),
                  0))
  }
  return(outer(x,j,fx)) #returns a (length(x) by j)-matrix, with different colums corresponding to j:s
}

# Integrals of hat functions
# k is scalar variable, [a,b] is interval of interest, n number of hat functions and j the index of integral
hat.int <- function(k,a,b,n,j){
  h<-(b-a)/(n-1)
  x<- a + (j-1)*h
  
  if(j==1){
    return(  exp(1i*k*a)*(1-exp(1i*k*h)+1i*k*h)/h/k^2) 
  }
  else if(j==n){
    return(  exp(1i*k*b)*(1-exp(-1i*k*h)-1i*k*h)/h/k^2) 
  }
  else
  {
    return(exp(1i*k*x)/(h*k^2)*(2-exp(-1i*k*h)-exp(1i*k*h)))
  }
}

