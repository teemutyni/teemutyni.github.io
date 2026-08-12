
##############################################################################
# This is the main script file for 1D scattering for fourth order operator.
# All the helper (utility) functions are defined in util.R and funcs.R.
#
# HOX! Set working directory to this folder before running.
#
# T. Tyni, M. Harju, Univ. Oulu 2016
##############################################################################

t1 <- Sys.time()

supp = c(0,0.5) # support of V(x)

source("util.R")

# Set up Simpson weights and abscissas for the support of V
M1 <- simpson(supp[1],supp[2],128)  

# Interval of interest is [-a,a]
a = 5

# Plot scattered wave with some wavenumber k
k = 2
plot(function(x) Re(u(x,k,M1)-u_0(x,k)),-a,a,ylab="")
plot(function(x) Im(u(x,k,M1)-u_0(x,k)),-a,a,col='red',add=TRUE,ylab="")
polygon(c(supp[1],supp[1],supp[2],supp[2]),c(-10,10,-10,10),col = adjustcolor("gray", alpha=0.5), border = NA) #corresponds to supp(V)
title(main = "Scattered wave")
legend( x="topleft",legend=c("real", "imag"),  col=c("black","red"),  lty=c(1,1) )

# Prepare for inverse problem by computing data using reflection coefficient
M <- 300		#number of values of k
kmax <- 140  #maximum value of k
kk <- seq(1,kmax,by=(kmax-1)/(M-1)) # vector of wavenumbers
gval <- g2(kk,M1)    #comp data and store result in vector gval

# Small values of k correspond to b(k)=0
k0 <- integrate(V,supp[1],supp[2])$value # pay attention to this, integrate() does not always perform well in R
gval[k <= 2*k0] = 0

# Plot data of inversion
plot(kk,Re(gval),type='l',ylab="",xlab="k",ylim=c(-3,2.5))
lines(Im(gval),col='red',type='l')
title(main="Scaled reflection coeff g(k)")
legend( x="topright",legend=c("real", "imag"),  col=c("black","red"),  lty=c(1,1) )

# Add noise to data
sigma = 0.01 # std of noise is sigma*max(abs(g))
gval = gval +  rnorm(length(gval), mean=0, sd=sigma*max(abs(gval))) +  1i*rnorm(length(gval), mean=0, sd=sigma*max(abs(gval)))

# Prepare matrix A for the linear system by integrating hat functions agains the exponentials
N <- 400		#number of grid points for hat functions
A <- matrix(0,M,N)
for(i in 1:M){
	for(j in 1:N){
		A[i,j] <- hat.int( kk[i], -a, a, N, j )
  }
}
rcond(A)


# Do inversion by TSVD method
s <- svd(A)			#svd of A
plot(s$d,log="y",ylab="Singular values")    #plot the singular values: smaller singular values might be omitted as 0's
stol <- 0.3
s$d[s$d < stol] <- 0	# truncation: replace small singular values by 0
L <- 1/diag(s$d)      # diagonal matrix consisting of reciprocals of singular values
L[L == Inf] <- 0	    # fill with zeros

fj <- s$v %*% L %*%  Conj( t(s$u) ) %*% gval #regularized solution f of Af=g
V_b <- function(x) -1/(2*pi)*Im(c(hat_j(x,1:N,-a,a,N)%*%fj))


# Plot reconstruction and potential
plot(V_b,-a,a, type="l",lty=2, xlab="x", col="red", ylab="", ylim=c(-4,4), n=600)
curve(V,add=TRUE,n=600,col="black")
title("Reconstruction")
legend( x="topleft", legend=c("V", "V_B"), col=c("black","red"), lty=c(1,2) )

cat("main.R run in", Sys.time()-t1 , "seconds\n")




#Plots that are not drawn unless explicitly wanted:
#Select and run the required lines
if(0){ #OPTIONAL PLOTS START
  #plots of u_j for k=5
  #u_0
  plot(function(x) Re(u_0(x,5)),-5,5,col="blue",ylab="")
  plot(function(x) Im(u_0(x,5)),-5,5,col="red",add=TRUE)
  abline(v=c(0,1),col="gray")
  title("u_0")
  legend( x="topleft", legend=c("Re", "Im"),col=c("blue","red"),lty=c(1,1) )
  
  #u_1
  plot(function(x) Re(u_j(x,5,1,M1)),-5,5,col="blue",ylab="")
  plot(function(x) Im(u_j(x,5,1,M1)),-5,5,col="red",add=TRUE)
  abline(v=c(0,1),col="gray")
  title("u_1")
  legend( x="topleft", legend=c("Re", "Im"),col=c("blue","red"),lty=c(1,1) )
  
  #u_2
  plot(function(x) Re(u_j(x,5,2,M1)),-5,5,col="blue",ylab="")
  plot(function(x) Im(u_j(x,5,2,M1)),-5,5,col="red",add=TRUE)
  abline(v=c(0,1),col="gray")
  title("u_2")
  legend( x="topleft", legend=c("Re", "Im"),col=c("blue","red"),lty=c(1,1) )
  
  #u_3
  plot(function(x) Re(u_j(x,5,3,M1)),-5,5,col="blue",ylab="")
  plot(function(x) Im(u_j(x,5,3,M1)),-5,5,col="red",add=TRUE)
  abline(v=c(0,1),col="gray")
  title("u_3")
  legend( x="topleft", legend=c("Re", "Im"),col=c("blue","red"),lty=c(1,1) )
  
  #u_4
  plot(function(x) Re(u_j(x,5,4,M1)),-5,5,col="blue",ylab="")
  plot(function(x) Im(u_j(x,5,4,M1)),-5,5,col="red",add=TRUE)
  abline(v=c(0,1),col="gray")
  title("u_4")
  legend( x="topleft", legend=c("Re", "Im"),col=c("blue","red"),lty=c(1,1) )
  
  #k^3*b(k)
  plot(Re(g.values),ylim=c(-2,2),col="blue",type="l",ylab="")
  lines(Im(g.values),col="red",type="l")
  title("Scaled reflection coefficient")
  legend( x="topright", legend=c("Re", "Im"),col=c("blue","red"),lty=c(1,1) )
  
  #Potential V
  plot(V,-5,5,ylim=c(-1,3),n=400,ylab="")
  title("Potential V(x)")
  
}#OPTIONAL PLOTS END



