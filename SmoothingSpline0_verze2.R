# VYHLAZUJICI SPLAJN

# Rozsirena sit je tvorena z puvodni site pridanim #(k-1) NASOBNYCH uzlu.

SmoothingSpline0 = function(knots,t,f,w,k,der,alfa,ch){
  # VSTUP:  knots= sit uzlu splajnu,
  #         t = body aproximace,
  #         f = funkcni hodnoty v bodech aproximace,
  #         w = vahove koeficienty,
  #         k = rad splajnu (order), degree = k-1
  #         der = derivace.
  #         alfa = vyhlazovac? parametr
  #         ch ={1,2};  ch=1: functional with (1-alfa) and alfa
  #                     ch=2: funkcional with alfa   
  
  r = length(knots) 
  library(splines)
  
  Celkova_Delka = 2*(k-1) + r 
  y = c() 
  for (i in 1:(Celkova_Delka)){
    if (i <= k-1){y[i] = knots[1]}
    if ((i > k-1) && (i <= r + k-1)){y[i] = knots[i-(k-1)]}
    if (i > r+ k-1){y[i] = knots[r]}
  }
  
  # Collocation matrix K
  K = splineDesign(y, t, k, outer.ok = TRUE)
  
  # Diag matrix with weights
  W = diag(w)
  
  # Collocation matrix C 
  deleni = seq(min(y), max(y),   length = 1000)    
  lambda = c(0:(r-1)) 
  g = lambda[length(lambda) - 1]
  # Dimension(space of splines)
  N = g+(k-1)+1
  C = array(0, c(length(deleni),N))
  l = c()
  for(i in (1:N)){
    for (j in 1:(k+1)){
      l[j] = y[i+j-1]
    }
    C[ ,i] = splineDesign(l, deleni, k, outer.ok = TRUE) 
  }
  # Verification of full column rank of collocation matrix K
  if (length(t) <= N) stop ('length(t) must be higher then Dimension(space of splines)')
  if (qr(K)$rank != N) stop ('Collocaton matrix does not have full column rank.')
  
  # Matice S
  S = array(0)
  if (der == 0){
    S = diag(1, c(N,N))
  }
  if (der > 0){
    i=der
    while (i>0){
      D = array(0)
      rozdil = y[(1+k):(N+k-i)] - y[(1+i):(N)]
      D = (k-i)*diag(1/rozdil)
      L = array(0, c(N-i,N-i+1))
      for (j in (1:(N-i))){
        L[j,j] = (-1)
        L[j,j+1] = 1
      }
      if (i==der){
        S = D%*%L
      } else{
        S = S%*%D%*%L
      }
      i=i-1
    }
  }
  
  # Matrix M -  order of spline = k-der
  kk = k-der 
  
  # Matrix M - augment knot sequence
  celkova_delka = 2*(kk-1) + r                                         
  Y = c()
  for (i in 1:celkova_delka){
    if (i <= (kk-1)){Y[i] = knots[1]}
    if ((i > kk-1) && (i <= r + kk-1)){Y[i] = knots[i-(kk-1)]}
    if (i > (r+(kk-1))){Y[i] = knots[r]}
  }  
  
  Deleni = seq(min(Y), max(Y),  length = 1000)    
  Lambda = c(0:(r-1)) 
  G = Lambda[length(Lambda) - 1]
  
  # Matrix M - spline space dimension
  NN = G+(kk-1)+1
  
  # Matrix M - collocation matrix KK
  CC = splineDesign(Y, Deleni, kk, outer.ok=TRUE)
  # Matrix M - function for computing integral
  SLP=function(krok, c){
    integral = krok*(0.5*c[1]+sum(c[2:(length(c)-1)]) +0.5*c[length(c)])
    return (integral)
  }
  
  krok=diff(Deleni[1:2])
  
  # Matrix M
  M=array(0, c(NN,NN))
  for (i in 1:NN){
    for (j in 1:NN){
      nenulove = c()
      soucin = CC[,i]*CC[,j]
      for (m in 1:length(Deleni)){
        if (soucin[m] != 0) {nenulove[m] = soucin[m]}
      }
      M[i,j]=SLP(krok, soucin)
    }
  }
  
  # Matrix D
  rozdil = y[(1+k):(r+2*(k-1))] - y[(1:(r+k-2))]
  D = (k)*diag(1/rozdil)
  
  # Matrix K
  KK = array(0, c(N,N-1))
  KK[1,1]=1
  KK[N,N-1]=-1
  
  for (j in (2:(N-1))){
    KK[j,j-1] = (-1)
    KK[j,j] = 1
  }
  KK
  
  # Matrix U
  U = S%*%D%*%KK
  
  # Matrix G, vector g
  if (ch==1){GG = t(U)%*%((1-alfa)*M)%*%U + alfa * t(KK)%*%t(D)%*%t(K)%*%W%*%K%*%D%*%KK}
  if (ch==2){GG = t(U)%*%M%*%U +            alfa * t(KK)%*%t(D)%*%t(K)%*%W%*%K%*%D%*%KK}
  
  gg = alfa*t(KK)%*%t(D)%*%t(K)%*%W%*%f
  
  # vector of B-spline coefficients := z
  z = solve(GG)%*%gg
  
  # B-spline basis in L20
  Bbasis = C%*%D%*%KK
  #matplot(deleni,Bbasis, type="l",lty=1,las=1, xlab="t", col = rainbow(dim(Bbasis)[2]))
  #abline(v=knots, col="gray",lty=2)
  
  # Resulting spline
  spline0 = (C%*%D%*%KK)%*%z
  
  #puvodni - Talska
   matplot(deleni,spline0, type="l",las=1,xlab="t",ylab="",col="darkblue",lwd=2,
          ylim = c(min(c(min(f),min(spline0))),max(c(max(f),max(spline0)))),
          main = paste("Spline of order k =",k))
  matpoints(t,f, pch = 8)
  abline(h=0,col="red",lty=2,lwd=1)
  
  #nove - Pavlu (prevedeni osy)
  #matplot(exp(deleni),spline0, type="l",log="x",las=1,ylab="",col="darkblue",lwd=2,
  #        ylim = c(min(c(min(f),min(spline0))),max(c(max(f),max(spline0)))),
  #        main = paste("Spline of order k =",k),xlab = expression (paste("Particle size (",mu,"m)")))
  #matpoints(exp(t),f, pch = 8)
  #abline(h=0,col="red",lty=2,lwd=1)
  
  #nove - Pavlu2 (prevedeni osy)
  #matplot(exp(deleni),spline0, type="l",log="x",las=1,ylab="clr density",col="darkblue",lwd=3,
  #        ylim = c(min(c(min(f),min(spline0))),max(c(max(f),max(spline0)))),
  #        main = paste("DO 805"),xlab = expression (paste("Particle size (",mu,"m)")))
  #matpoints(exp(t),f, pch = 8)
  #abline(h=0,col="red",lty=2,lwd=1)
  
  integral = SLP(krok,spline0)
  
  if (ch==1) {J = (1-alfa)*t(z)%*%t(U)%*%M%*%U%*%z + alfa*t(f-K%*%D%*%KK%*%z)%*%W%*%(f-K%*%D%*%KK%*%z)}
  if (ch==2) {J =          t(z)%*%t(U)%*%M%*%U%*%z + alfa*t(f-K%*%D%*%KK%*%z)%*%W%*%(f-K%*%D%*%KK%*%z)}
  return(list(J=J,z=z,spline0=spline0))
}

# Smoothing Italy income data
# t = c(6574, 19591, 32608, 45625, 58641, 71658, 84675, 97692, 110709)
# f = c(0.58724119824843, 2.33068209020568, 2.15428355171484, 1.27090963485383, 0.330521351400616, -0.549837371247476, -1.43714056624838, -1.9967563541838, -2.68990353474375)
# knots = c(0, 30000, 70000, 117218)
# w = rep(1,9)
# k = 4
# der = 2
# alfa =0.00001
# ch = 1

# Spline = SmoothingSpline0(knots=knots, t=t.mid, f=f, w=w, k=k, der=der, alfa=alfa, ch=ch)
# Spline$z
# Spline$J

# J = Spline[[1]];J
# z = Spline[[2]];z