# Baze s nulovym integralem
ZsplineBasis = function(knots,k)
{
  library(fda)
  r = length(knots)
  lambda_index = c(0:(r-1)) 
  g = lambda_index[length(lambda_index) - 1]
  
  N = g+(k-1)+1
  
  lambda = c(rep(min(knots),k-1),knots,rep(max(knots),k-1))
  deleni = seq(min(lambda), max(lambda), length = 1000) 
  
  # standard B-spline basis; collocation matrix := C
  splajn.basis = create.bspline.basis(range(knots),nbasis = N , norder = k, breaks = knots)
  C = eval.basis(deleni, splajn.basis)
  
  # Matrix D
  rozdil = lambda[(1+k):(r+2*(k-1))] - lambda[(1:(r+k-2))]
  D = (k)*diag(1/rozdil)
  
  # Matrix L
  L = array(0, c(N,N-1))
  L[1,1]=1
  L[N,N-1]=-1
  
  for (j in (2:(N-1))){
    L[j,j-1] = (-1)
    L[j,j] = 1
  }
  
  # Spline0 basis: collocation matrix C0
  C0 = C%*%D%*%L
  #matplot(deleni,C0, type="l",lty=1, las=1, col=rainbow(N-1), xlab="t", ylab="Z-spline basis")
  #abline(v=knots, col="gray", lty=2)
  
  # Matrix M - function for computing integral
  SLP=function(krok, c){
    integral = krok*(0.5*c[1]+sum(c[2:(length(c)-1)]) +0.5*c[length(c)])
    return (integral)
  }
  
  Deleni = seq(min(lambda), max(lambda),  length = 1000);krok=diff(Deleni[1:2])
  CC = eval.basis(Deleni, splajn.basis)
  
  CC0 = CC%*%D%*%L
  
  M=array(0, c(N-1,N-1))
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      nenulove = c()
      soucin = CC0[,i]*CC0[,j]
      for (m in 1:length(Deleni)){
        if (soucin[m] != 0) {nenulove[m] = soucin[m]}
      }
      M[i,j]=SLP(krok, soucin)
    }
  }
  return(list(C0 = C0, M0 = M, K = L, D = D))
}

#ZsplineBasis(knots,k)
