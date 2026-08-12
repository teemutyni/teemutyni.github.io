
##############################################################################
# This file defines all the helper functions which might be called in parallel.
# Otherwise they could be kept in util.R
#
# T. Tyni, M. Harju, Univ. Oulu 2016
##############################################################################

# potential function V, requires supp to be set before this
V<-function(x){
  idx <- 1*(supp[1] < x & x < supp[2] ) #characteristic function of [0,0.5]
  return(.5*idx)
}

# Iterations u_j(x,k) using Simpson weights and nodes in M1 matrix
u_0<-function(x,k){
  return(exp(1i*k*x))
}
u_j<-function(x,k,j,M1){
  if(j==0){ return(u_0(x,k)) }else{
    G <- function(x,k){
      return( -1/(4*k^3)*(1i*exp(1i*k*abs(x)) - exp(-k*abs(x))) )
      }
    f<-function(x,y){
      return( t(G( outer(x,y,"-") ,k ))*V(y)*u_j(y,k,j-1,M1) )
      }
    return( colSums(M1[2,]*f(x,M1[1,])) )
  }
}

# Total wave u(x,k) as a finite sum. Use more terms for better accuracy but longer computations.
u<-function(x,k,M1) u_0(x,k) + u_j(x,k,1,M1) + u_j(x,k,2,M1) + u_j(x,k,3,M1) + u_j(x,k,4,M1) + u_j(x,k,5,M1)