# THESE FUNCTIONS IS FOR CONSTRUCTING UNBLOCKED DESIGNS BY 
# USING BAYESIAN CRITERION D - COORDINATE EXCHANGES
# FUNCTION FOR GENERATING INITIAL DESIGN (COORDINATE EXCHANGE) 
DesigniCoord <- function(N,Nlev,Terms,K)
{
  fac <- matrix(1:K,nc=1)
  d0 <- matrix(apply(fac,1,SampleFac,Nlev,N),nr=N) 
  d1 <- matrix(c(rep(1,N)),N,1) 
  dt<-cbind(d1,d0)
  list(d0=d0,dt=dt)
}

# GENERATES A RANDOM COLUMN FOR A FACTOR (COORDINATE EXCHANGE)
SampleFac <- function(f,Nlev,N)
{
  sample(seq(0.1,1,length=Nlev[f]),N,replace=T)
}


####### Utility ###################
# m=1, k=1 

Det_no_gama<-function(alfa, X){
  if(alfa==0){
     D<-length(X[,1])*sum(log(X[,2])^2)- sum(log(X[,2]))*sum(log(X[,2]))
  }
  else{
     D<-(
       sum(((X[,2]^alfa)*log(X[,2]))^2)*(sum(X[,2]^(2*alfa))*length(X[,1]) - sum(X[,2]^(alfa))*sum(X[,2]^(alfa)))+
        2*(sum(X[,2]^alfa)*sum((X[,2]^(2*alfa))*log(X[,2]))*sum((X[,2]^alfa)*log(X[,2])))-
         (sum(X[,2]^(2*alfa))*(sum((X[,2]^alfa)*log(X[,2])))^2) - length(X[,1])*(sum((X[,2]^(2*alfa))*log(X[,2])))^2
    )
  }
  return(log(D))
}
#### DISTRIBUTIONS FOR ALPHA  ###########
# - UNCOMMENT THE DESIRED DISTRIBUTION #

## SYMMETRIC ##
f_alfa<-function(x){
  for (i in length(x)){
    if(x==-1){
      p_a<-0.1
    }
    else if(x==-0.5){
      p_a<-0.2
    }
    else if(x==0){
      p_a<-0.4
    }
    else if(x==0.5){
      p_a<-0.2
    }
    else if(x==1){
      p_a<-0.1
    }
    else{
      p_a<-0
    }
  }
  return(p_a)
}

## UNIFORM ##
# f_alfa<-function(x){
#   for (i in length(x)){
#     if(x==-1){
#       p_a<-0.2
#     }
#     else if(x==-0.5){
#       p_a<-0.2
#     }
#     else if(x==0){
#       p_a<-0.2
#     }
#     else if(x==0.5){
#       p_a<-0.2
#     }
#     else if(x==1){
#       p_a<-0.2
#     }
#     else{
#       p_a<-0
#     }
#   }
#   return(p_a)
# }

## LEFT SKEWED ##
# f_alfa<-function(x){
#   for (i in length(x)){
#     if(x==-1){
#       p_a<-0.03
#     }
#     else if(x==-0.5){
#       p_a<-0.07
#     }
#     else if(x==0){
#       p_a<-0.15
#     }
#     else if(x==0.5){
#       p_a<-0.3
#     }
#     else if(x==1){
#       p_a<-0.45
#     }
#     else{
#       p_a<-0
#     }
#   }
#   return(p_a)
# }

## RIGHT SKEWED ##
# f_alfa<-function(x){
#   for (i in length(x)){
#     if(x==-1){
#       p_a<-0.45
#     }
#     else if(x==-0.5){
#       p_a<-0.3
#     }
#     else if(x==0){
#       p_a<-0.15
#     }
#     else if(x==0.5){
#       p_a<-0.07
#     }
#     else if(x==1){
#       p_a<-0.03
#     }
#     else{
#       p_a<-0
#     }
#   }
#   return(p_a)
# }



f_gama_1<-function(x){
  dnorm(x,1,0.25)
}


f_gama_comp<-function(x){
  log(x^2)* f_gama_1(x)
}


Utility<-function(x, xmin,xmax){
  sum_elements<-0
  list_alfas<-c(-2,-1,-0.5,0,0.5,1,2)
  for(alfa in list_alfas){
    sum_elements<-sum_elements+f_alfa(alfa)*(Det_no_gama(alfa,x))
  }
  integrate_gama_p1<- integrate(f_gama_comp, lower = -0.75, upper = -0.001)
  integrate_gama_p2<- integrate(f_gama_comp, lower = 0.001, upper = 0.75)
  integrate_gama <- integrate_gama_p1$value + integrate_gama_p2$value
  U<-integrate_gama+sum_elements
  U<-sum_elements
  return(U)
}


###############################################

# CONTROLS FOR NON-SINGULAR INITIAL DESIGN (COORDINATE EXCHANGE)
SampleDcoord <- function(K,N,Nlev,Terms,pinv, X,m)
{
  DC<-DesigniCoord(N,Nlev,Terms,K)
  Xagg<-DC$dt
  xmin<-0.1
  xmax<-1
  d<-Utility(Xagg, xmin,xmax)
  list(d=d, Xagg=Xagg) 
}


# LABEL FOR TREATMENTS (POINT EXCHANGE)
TreatLabels <- function(x)
{
  n <- nrow(x)
  Treat <- matrix(c(1,rep(0,n-1)),n,1)
  Label <- 1
  for (i in 2:n)
  {
    Label <- Label+1
    Treat[i] <- Label
    for (j in 1:(i-1))
    {
      if(min(as.numeric(x[i,]==x[j,]))==1) Treat[i] <- Treat[j]
    }
    if(Treat[i] < Label) Label <- Label-1
  }
  Treat <- as.matrix(cbind(Treat,x))
  list(Treat=Treat)
}


#FUNCTION TO ADD ELEMENTS TO LIST
add_element<-function(lista,elemento){
  ind<-which(sapply(lista, is.null))[1]
  
  if(!is.na(ind)){
    lista[[ind]]<-elemento
  }else{
    lista<-append(lista,list(elemento))
  }
  return(lista)
}

# FUNCTION TO DRIVE THE SEARCH ---
SearchTreat <- function(K,Levels,Cubic,N,Ntries, m, 
                        Terms,Npar,AlgType)

{
  inicio <- Sys.time()
  Power <- 1/Npar
  critall <- matrix(0,nr=Ntries)
  Xall<- vector("list",Ntries)
  
  # COORD EXCHANGE BY ROW 
  if(AlgType=="CR")
  {
    Nlev <- matrix(0,nr=K)
    for(i in 1:K)
    {
      Nlev[i] <- max(Levels[[i]])
    }
    for(nt in 1:Ntries) 
    {    
      Des <- SampleDcoord(K,N,Nlev,Terms,Power, X,m) 
      
      Xagg<-Des$Xagg 
      
      Util <- Des$d
      crite <- 0
      criteD <- 0
      
      
      
      criteD <- Util
      crite <- criteD 
      
      improve <- 1 
      
      iteration<-0
      while(improve==1)
      {         
        Xs <- CoordRow(Xagg,Util,iteration,Xall[nt],crite,K,N,Nlev,Terms,Npar,Power, m)  
        Xj<-Xs$D 
        crite <- Xs$crita
        improve <- Xs$improve
        Util <- Xs$Util
        iteration<-Xs$iter
        Xall[[nt]]<-Xs$Xall
        
      }
      critall[nt] <- crite
      if(nt==1)  #first try
      {
        critopt <- crite
        Xopt<-Xj
      }else
      {
        if(crite > critopt)
        {
        
          Xopt<-Xj
          critopt <- crite
          
        }
      }
    }
    Xopt <- TreatLabels(Xopt)$Treat 
  }
  
  # COORD EXCHANGE BY COLUMN
  if(AlgType=="CC")
  {
    Nlev <- matrix(0,nr=K)
    for(i in 1:K)
    {
      Nlev[i] <- max(Levels[[i]])
    }
    for(nt in 1:Ntries) 
    {    
      Des <- SampleDcoord(K,N,Nlev,Terms,Power, X,m) 
      Xagg<-Des$Xagg
      Util <- Des$d
      crite <- 0
      criteD <- 0
       
      criteD <- Util
      crite <- criteD
      improve <- 1 
      iteration<-0
      while(improve==1)
      {         
        Xs <- CoordCol(Xagg,Util,iteration,Xall[nt],crite,K,N,Nlev,Terms,W,Npar,Power,Ncali, m) 
        Xj <- Xs$D
        crite <- Xs$crita
        improve <- Xs$improve
        
        Util <- Xs$Util
        iteration<-Xs$iter
        Xall[[nt]]<-Xs$Xall
      }
      critall[nt] <- crite
      if(nt==1)
      {
        critopt <- crite
        Xopt<-Xj
      }else
      {
        if(crite > critopt)
        {
          Xopt<-Xj
          critopt <- crite
        }
      }
    }
    Xopt <- TreatLabels(Xopt)$Treat
  }
  
  tfim <- Sys.time()-inicio
  dfbest <- N-nlevels(as.factor(Xopt[,1]))
  xmin<-0.1
  xmax<-1
  Util<-Utility(Xopt[,-1],xmin,xmax)
  list(Xopt=Xopt, critopt=critopt, 
       dfbest=dfbest,Util=Util,
       critall=critall,inicio=inicio,tfim=tfim,Xall=Xall,iter=iteration) 
}


# FUNCTION TO EXCHANGE COORDINATES BY ROW 
CoordRow <- function(D,Util,iteration,Xall,crita,K,N,Nlev,Terms,Npar,pinv, m) #linha 7 a 17?  Xr,D
{
  improve <- 0
  for(i in 1:N)
    
  {
    calibrate <- 0
    for(j in 1:K)
    {
      for(x in seq(0.1,1,length=Nlev[j]))
      {
        if(D[i,(j+1)]!=x)
        { 
          xout <- matrix(D[i,],nc=1)
          xin <- t(xout[2:(K+1)])
          xin[j] <- x
          xin_tot<-matrix(c(1,xin),nc=1)
          Xc <- D
          Xc[i,] <- t(xin_tot)
          critc <- 0
          ns <- 0
          xmin<-0.1
          xmax<-1
          
          Utilup <- Utility(Xc,xmin,xmax)

          critD <- Utilup
            
          critc <- critD
          
          if(critc>crita)
          {
            D <- Xc 
            crita <- critc
            busca <- 1
            Util <- Utilup
            Improve <- 1
            iteration<-iteration+1
            Xall<-add_element(Xall,list(iteration = iteration, criterio = crita)) 
          } 
        }
      }
    }
  }
  list(D=D,crita=crita,improve=improve, Util=Util, Xall=Xall) 
}

# FUNCTION TO EXCHANGE COORDINATES BY COLUMN         
CoordCol<- function(D,Util,iteration,Xall,crita,K,N,Nlev,Terms,Npar,pinv, m)
{
  improve <- 0
  for(j in 1:K)
  {
    calibrate <- 0
    for(i in 1:N)
    {
      for(x in seq(0.1,1,length=Nlev[j]))
      {
        if(D[i,(j+1)]!=x) 
        { 
          xout <- matrix(D[i,],nc=1)
          xin <- t(xout[2:(K+1)])
          xin[j] <- x
          xin_tot<-matrix(c(1,xin),nc=1)
          Xc <- D
          Xc[i,] <- t(xin_tot)
          critc <- 0
          ns <- 0
          xmin<-0.1
          xmax<-1
          Utilup <- Utility(Xc,xmin,xmax)
          critD <- Utilup
          critc <- critD
          
          if(critc>crita)
          {
            D <- Xc
            crita <- critc
            busca <- 1
            Util <- Utilup
            Improve <- 1
            iteration<-iteration+1
            Xall<-add_element(Xall,list(iteration = iteration, criterio = crita))
          } 
        }
      }
    }
  }
  list(D=D,crita=crita,improve=improve,Util=Util, Xall=Xall)
}



