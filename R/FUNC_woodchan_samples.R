### Required libraries:
library(fields)
library(FastGP)

#Matern kernel with smoothness nu and length - scale l:
MK  =  function(x,  y , l,  nu)
{
  nu = ifelse(abs(x - y) > 0,  (sqrt(2 * nu) * abs(x - y) / l)^nu / (2^(nu - 1) * gamma(nu)) * besselK(x = abs(x - y) * sqrt(2 * nu) / l,  nu = nu),  1.0)

  return(nu)
}

# Order of the circulant matrix:
# minimum value of g and m so that G can be embedded into C
min_g = function(knot)
{
  N = length(knot)
  g = ceiling(log(2 * N, 2))   #m = 2^g and m >  = 2(n - 1) : Wood & Chan notation; 
  #since we are going upto n and not stopping at (n - 1),  the condition is modified!
  return("g" =g)
}

# forming the circulant matrix:
circulant = function(x)
{
  n = length(x)
  mat = matrix(0,  n,  n)
  for (j in 1:n) 
  {
    mat[j,  ] =  c(x[ - (1:(n + 1 - j))],  x[1:(n + 1 - j)])
  }
  return(mat)
}

# Function for forming the vector of circulant matrix:
circ_vec = function(knot, g, nu, l, tausq)
{
  delta_N = 1 / (length(knot) - 1)
  m = 2**g
  cj = integer()
  for (j in 1:m)
  {
    if (j <= (m / 2))
      cj[j] = (j - 1) * delta_N
    else
      cj[j] = (m - (j - 1)) * delta_N
  }
  x = (tausq * MK(cj, 0, l, nu))
  return(x)
}

# Function for finding a g such that C is nnd:
eig.eval=function(knot,g,nu,l,tausq){
  vec=circ_vec(knot,g,nu,l,tausq)
  C=circulant(vec)
  ev=min(eigen(C)$values)
  return(list("vec" = vec, "min.eig.val" = ev))
}

# Function for finding a g such that C is nnd:
# without forming the circulant matrix and without computing eigen values:
C.eval=function(knot,g,nu,l,tausq){
  vec=circ_vec(knot,g,nu,l,tausq)
  val=fft(vec) # eigenvalues will be real as the circulant matrix formed by the 
  # vector is by construction is symmetric!
  ev=min(Re(val))
  return(list("vec" = vec, "min.eig.val" = ev))
}


nnd_C=function(knot,g,nu,l,tausq){
  C.vec=C.eval(knot,g,nu,l,tausq)$vec
  eval=C.eval(knot,g,nu,l,tausq)$min.eig.val
  if(eval>0)
    return(list("cj" = C.vec,"g" = g))
  else{
    g=g+1
    nnd_C(knot,g,nu,l,tausq)
  }
}


# computing the eigen values of C using FFT:
eigval=function(knot,nu,l,tausq){
  g=min_g(knot)
  c.j=nnd_C(knot,g,nu,l,tausq)$cj
  lambda=Re(fft(c.j))
  if(min(lambda)>0)
    return(lambda)
  else
    stop("nnd condition is NOT satisfied!!!")
}


#################################################################
########## Samples drawn using Wood and Chan Algorithm ##########
#################################################################
# N: (N + 1) is the number of knots
# nu:smoothness parameter of Matern
# l:length - scale parameter of Matern
samp.WC=function(knot, l, nu = 5/2, tausq = 1){
  N=length(knot)
  lambda=eigval(knot,nu,l,tausq)
  m=length(lambda)
  samp.vec=rep(0,N)
  a=rep(0,m)
  a[1]=sqrt(lambda[1])*rnorm(1)/sqrt(m)
  a[(m/2)+1]=sqrt(lambda[(m/2)+1])*rnorm(1)/sqrt(m)
  i=sqrt(as.complex(-1))
  for(j in 2:(m/2)){
    uj=rnorm(1); vj=rnorm(1)
    a[j]=(sqrt(lambda[j])*(uj + i*vj))/(sqrt(2*m))
    a[m+2-j]=(sqrt(lambda[j])*(uj - i*vj))/(sqrt(2*m))
  }
  samp=fft(a)
  samp.vec=Re(samp[1:N])
  return(samp.vec)
}


