# THESE FUNCTIONS IS FOR CONSTRUCTING UNBLOCKED DESIGNS BY 
# USING CRITERION D - POINT AND COORDINATE EXCHANGES
# FUNCTION FOR GENERATING INITIAL DESIGN (COORDINATE EXCHANGE) 
DesigniCoord <- function(N,Nlev,Terms,K, alfa)
{
  fac <- matrix(1:K,nc=1)
  d0 <- matrix(apply(fac,1,SampleFac,Nlev,N),nr=N) 
  di <- ModelMatrix(d0,Terms,K, alfa)$X
  dt<-cbind(di[,1],d0,di[,-1]) 
  list(d0=d0,di=di,dt=dt)
}

# GENERATES A RANDOM COLUMN FOR A FACTOR (COORDINATE EXCHANGE)
SampleFac <- function(f,Nlev,N)
{
  sample(seq(0.1,1,length=Nlev[f]),N,replace=T)
}

####################################################
# INFORMATION MATRIX FUNCTION #
# m - polynomial order
# K - number of factors
# beta - coefficients vector
# alfa - power vector - for fat = 1, alfa is numeric

# X[,1] - intercept
# X[,2] - refers to x1 for 2 factors and x for 1 factor
# X[,3] - refers to x2 for 2 factors 

fisher_matrix <- function(X,m,K, beta, alfa){
  if (m==1 & K ==1){
    if(alfa[1]!=0){
      M<- matrix(c(length(X[,1]),sum(X[,2]^alfa[1]), beta[2]*sum((X[,2]^alfa[1])*log(X[,2])),
                   sum(X[,2]^alfa[1]), sum(X[,2]^(2*alfa[1])), beta[2]*sum(X[,2]^(2*alfa[1])*log(X[,2])),
                   beta[2]*sum((X[,2]^alfa[1])*log(X[,2])), beta[2]*sum(X[,2]^(2*alfa[1])*log(X[,2])), beta[2]^2*sum((X[,2]^(2*alfa[1]))*(log(X[,2])^2))),
                 byrow = TRUE, nrow = 3, ncol = 3) 
    }
    else if(alfa[1]==0){
      M<-matrix(c(length(X[,1]),sum(log(X[,2])),
                  sum(log(X[,2])),sum(log(X[,2])^2)), byrow = TRUE, nrow = 2, ncol = 2) 
      }
    
  }
  
  else if (m==1 & K ==2){
    if(alfa[1]!=0 & alfa[2]!=0){
      M<-matrix(c(length(X[,1]),sum(X[,2]^alfa[1]), sum(X[,3]^alfa[2]), sum((X[,2]^alfa[1])*(X[,3]^alfa[2])),beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*(log(X[,2]))), beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+beta[4]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*(log(X[,3]))),
                sum(X[,2]^alfa[1]), sum(X[,2]^(2*alfa[1])), sum((X[,2]^alfa[1])*(X[,3]^alfa[2])), sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])), beta[3]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,3]))+ beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])* log(X[,3])),
                sum(X[,3]^alfa[2]), sum((X[,2]^alfa[1])*(X[,3]^alfa[2])), sum(X[,3]^(2*alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))), beta[2]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))* log(X[,2])), beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[4]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))* log(X[,3])),
                sum((X[,2]^alfa[1])*(X[,3]^alfa[2])), sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))), sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))), beta[2]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2]))+beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])), beta[3]*sum((X[,3]^(2*alfa[2]))*(X[,2]^alfa[1])*log(X[,3]))+beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])),
                
                beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])* log(X[,2])), beta[2]*sum((X[,2]^(alfa[1]))*(X[,3]^(alfa[2]))*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))* log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2]))+beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])),
                beta[2]^2*sum((X[,2]^(2*alfa[1]))*log(X[,2])^2)+ 2*beta[2]*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])^2)+ (beta[4]^2)*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])^2),beta[2]*beta[3]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2])*log(X[,3]))+ beta[2]*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3]))+ beta[4]*beta[3]*sum((X[,2]^alfa[1])*log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[4]^2*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])),  
                
                beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+ beta[4]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,3])),  beta[3]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,3]))+ beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])* log(X[,3])), beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[4]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))* log(X[,3])), beta[3]*sum((X[,3]^(2*alfa[2]))*(X[,2]^alfa[1])*log(X[,3]))+beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])),
                beta[2]*beta[3]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2])*log(X[,3]))+ beta[2]*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3]))+ beta[4]*beta[3]*sum((X[,2]^alfa[1])*log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[4]^2*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])),  
                (beta[3]^2)*sum((X[,3]^(2*alfa[2]))*log(X[,3])^2)+2*beta[3]*beta[4]*sum((X[,3]^(2*alfa[2]))*(X[,2]^alfa[1])*log(X[,3])^2)+ (beta[4]^2)*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])^2)
    ),byrow = TRUE, nrow = 6, ncol = 6)
    }
    else if(alfa[1]==0 & alfa[2]!=0){
      M<-matrix(c(length(X[,1]),sum(log(X[,2])), sum(X[,3]^alfa[2]), sum(log(X[,2])*(X[,3]^alfa[2])), beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+beta[4]*sum(log(X[,2])*(X[,3]^alfa[2])*(log(X[,3]))),
                  sum(log(X[,2])), sum(log(X[,2])^2), sum(log(X[,2])*(X[,3]^alfa[2])), sum(log(X[,2])^2*(X[,3]^alfa[2])),beta[3]*sum(log(X[,2])*(X[,3]^alfa[2])*log(X[,3]))+ beta[4]*sum((log(X[,2])^2)*(X[,3]^alfa[2])* log(X[,3])),
                  sum(X[,3]^alfa[2]), sum(log(X[,2])*(X[,3]^alfa[2])), sum(X[,3]^(2*alfa[2])), sum(log(X[,2])*(X[,3]^(2*alfa[2]))), beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[4]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))* log(X[,3])),
                  sum(log(X[,2])*(X[,3]^alfa[2])), sum((log(X[,2]))^2*(X[,3]^alfa[2])), sum(log(X[,2])*(X[,3]^(2*alfa[2]))), sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))),  beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3]))+beta[4]*sum(log(X[,2])^2*(X[,3]^(2*alfa[2]))*log(X[,3])),
                  beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+ beta[4]*sum(log(X[,2])*(X[,3]^alfa[2])*log(X[,3])),  beta[3]*sum(log(X[,2])*(X[,3]^alfa[2])*log(X[,3]))+ beta[4]*sum(log(X[,2])^2*(X[,3]^alfa[2])* log(X[,3])), beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[4]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))* log(X[,3])), beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3]))+beta[4]*sum(log(X[,2])^2*(X[,3]^(2*alfa[2]))*log(X[,3])),
                   
                  (beta[3]^2)*sum((X[,3]^(2*alfa[2]))*log(X[,3])^2)+2*beta[3]*beta[4]*sum((X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])^2)+ (beta[4]^2)*sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))*log(X[,3])^2)
                  ),byrow = TRUE, nrow = 5, ncol = 5)            
    }
    else if(alfa[1]!=0 & alfa[2]==0){
      M<-matrix(c(length(X[,1]),sum(X[,2]^alfa[1]), sum(log(X[,3])), sum((X[,2]^alfa[1])*log(X[,3])),beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*log(X[,3])*(log(X[,2]))), 
                  sum(X[,2]^alfa[1]), sum(X[,2]^(2*alfa[1])), sum((X[,2]^alfa[1])*log(X[,3])), sum((X[,2]^(2*alfa[1]))*log(X[,3])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2])),
                  sum(log(X[,3])), sum((X[,2]^alfa[1])*log(X[,3])), sum(log(X[,3])^2), sum((X[,2]^alfa[1])*(log(X[,3])^2)), beta[2]*sum((X[,2]^alfa[1])*log(X[,3])*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*(log(X[,3])^2)* log(X[,2])),
                  sum((X[,2]^alfa[1])*log(X[,3])), sum((X[,2]^(2*alfa[1]))*log(X[,3])), sum((X[,2]^alfa[1])*log(X[,3])^2), sum((X[,2]^(2*alfa[1]))*log(X[,3])^2), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2]))+beta[4]*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])),
                  
                  beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*log(X[,3])*log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,3])* log(X[,2])), beta[2]*sum((X[,2]^(alfa[1]))*log(X[,3])*log(X[,2]))+ beta[4]*sum((X[,2]^alfa[1])*log(X[,3])^2* log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2]))+beta[4]*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])),
                  (beta[2]^2)*sum((X[,2]^(2*alfa[1]))*log(X[,2])^2)+ 2*beta[2]*beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2])^2)+ (beta[4]^2)*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])^2)
                  ),byrow = TRUE, nrow = 5, ncol = 5)  
    }
    else if(alfa[1]==0 & alfa[2]==0){
      M<-matrix(c(length(X[,1]),sum(log(X[,2])), sum(log(X[,3])), sum(log(X[,2])*log(X[,3])),
                  sum(log(X[,2])), sum(log(X[,2])^2), sum(log(X[,2])*log(X[,3])), sum(log(X[,2])^2*log(X[,3])),
                  sum(log(X[,3])), sum(log(X[,2])*log(X[,3])), sum(log(X[,3])^2), sum(log(X[,2])*(log(X[,3])^2)),
                  sum(log(X[,2])*log(X[,3])), sum((log(X[,2])^2)*log(X[,3])), sum(log(X[,2])*log(X[,3])^2), sum((log(X[,2])^2)*log(X[,3])^2), 
                  ),byrow = TRUE, nrow = 4, ncol = 4) 
                  }
  }
  else if (m==2 & K ==1){
    if(alfa[1]!=0){
      M<-matrix(c(length(X[,1]),sum(X[,2]^alfa[1]), sum(X[,2]^(2*alfa[1])), beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+2*beta[3]*sum((X[,2]^(2*alfa[1]))*log(X[,2])),
                  sum(X[,2]^alfa[1]), sum(X[,2]^(2*alfa[1])),  sum(X[,2]^(3*alfa[1])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+2*beta[3]*sum((X[,2]^(3*alfa[1]))*log(X[,2])),
                  sum(X[,2]^(2*alfa[1])),  sum(X[,2]^(3*alfa[1])), sum(X[,2]^(4*alfa[1])), beta[2]*sum((X[,2]^(3*alfa[1]))*log(X[,2]))+2*beta[3]*sum((X[,2]^(4*alfa[1]))*log(X[,2])),
                  beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+2*beta[3]*sum((X[,2]^(2*alfa[1]))*log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+2*beta[3]*sum((X[,2]^(3*alfa[1]))*log(X[,2])), beta[2]*sum((X[,2]^(3*alfa[1]))*log(X[,2]))+2*beta[3]*sum((X[,2]^(4*alfa[1]))*log(X[,2])),  
                  (beta[2]^2)*sum((X[,2]^(2*alfa[1])*(log(X[,2]))^2))+4*beta[2]*beta[3]*sum((X[,2]^(3*alfa[1]))*log(X[,2])^2)+4*(beta[3]^2)*sum((X[,2]^(4*alfa[1]))*(log(X[,2]))^2)
      ),byrow = TRUE, nrow = 4, ncol = 4)
      
    }
    else if(alfa[1]==0){
      M<-matrix(c(length(X[,1]),sum(log(X[,2])), sum(log(X[,2])^2), 
                  sum(log(X[,2])), sum(log(X[,2])^2),  sum(log(X[,2])^3),
                  sum(log(X[,2])^2), sum(log(X[,2])^3), sum(log(X[,2])^4), 
                  ),byrow = TRUE, nrow = 4, ncol = 4)
    }
    
  }
  else if (m==2 & K ==2){
    if(alfa[1]!=0 & alfa[2]!=0){
      M<-matrix(c(length(X[,1]),sum(X[,2]^alfa[1]),sum(X[,3]^alfa[2]),sum(X[,2]^(2*alfa[1])),sum(X[,3]^(2*alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^alfa[2])), beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2])), beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+ 2*beta[5]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[6]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,3])), 
                sum(X[,2]^alfa[1]), sum(X[,2]^(2*alfa[1])), sum((X[,2]^alfa[1])*(X[,3]^alfa[2])),sum(X[,2]^(3*alfa[1])), sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))), sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])), beta[3]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,3])) + 2*beta[5]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,3])),
                sum(X[,3]^alfa[2]), sum((X[,2]^alfa[1])*(X[,3]^alfa[2])), sum(X[,3]^(2*alfa[2])), sum(X[,2]^(2*alfa[1])*(X[,3]^alfa[2])), sum(X[,3]^(3*alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))), beta[2]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2]))+ beta[6]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])), beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])),
                sum(X[,2]^(2*alfa[1])), sum(X[,2]^(3*alfa[1])), sum(X[,2]^(2*alfa[1])*(X[,3]^alfa[2])), sum(X[,2]^(4*alfa[1])), sum(X[,2]^(2*alfa[1])*(X[,3]^(2*alfa[2]))), sum(X[,2]^(3*alfa[1])*(X[,3]^alfa[2])), beta[2]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(4*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])), beta[3]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(alfa[2]))*log(X[,3]))+ 2*beta[5]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])*log(X[,3])),
                sum(X[,3]^(2*alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))), sum(X[,3]^(3*alfa[2])), sum(X[,2]^(2*alfa[1])*(X[,3]^(2*alfa[2]))), sum(X[,3]^(4*alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^(3*alfa[2]))), beta[2]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])) + beta[6]*sum((X[,2]^(alfa[1]))*(X[,3]^(3*alfa[2]))*log(X[,2])), beta[3]*sum((X[,3]^(3*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,3]^(4*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(alfa[1]))*(X[,3]^(3*alfa[2]))*log(X[,3])),
                sum((X[,2]^alfa[1])*(X[,3]^alfa[2])), sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))),sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])), sum((X[,2]^alfa[1])*(X[,3]^(3*alfa[2]))), sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))), beta[2]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(alfa[2]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*(X[,3]^(alfa[2]))*log(X[,2]))+ beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])), beta[3]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,2]^(alfa[1]))*(X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])),
                
                beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])),  beta[2]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2]))+ beta[6]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])), beta[2]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(4*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])),
                beta[2]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])) + beta[6]*sum((X[,2]^(alfa[1]))*(X[,3]^(3*alfa[2]))*log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(alfa[2]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*(X[,3]^(alfa[2]))*log(X[,2]))+ beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])), (beta[2]^2)*sum((X[,2]^(2*alfa[1]))*log(X[,2])^2)+4*beta[2]*beta[4]*sum((X[,2]^(3*alfa[1]))*(log(X[,2])^2))+ 2*beta[2]*beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*(log(X[,2])^2))+ 4*(beta[4]^2)*sum((X[,2]^(4*alfa[1]))*(log(X[,2])^2))+4*beta[4]*beta[6]*sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])*(log(X[,2])^2))+ (beta[6]^2)*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*(log(X[,2])^2)),                  beta[2]*beta[3]*sum((X[,2]^(alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3])) + 2*beta[2]*beta[5]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])) +
                  beta[2]*beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3]))+ 2*beta[3]*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3])) + 4*beta[4]*beta[5]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])) + 2*beta[4]*beta[6]*sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3])) + beta[6]*beta[3]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])) + 2*beta[6]*beta[5]*sum((X[,2]^(alfa[1]))*(X[,3]^(3*alfa[2]))*log(X[,2])*log(X[,3])) + beta[6]^2*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])) ,
                
                beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+ 2*beta[5]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[6]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,3])), beta[3]*sum((X[,2]^alfa[1])*(X[,3]^alfa[2])*log(X[,3])) + 2*beta[5]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,3])),  beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])), beta[3]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(alfa[2]))*log(X[,3]))+ 2*beta[5]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])*log(X[,3])),
                beta[3]*sum((X[,3]^(3*alfa[2]))*log(X[,3]))+ 2*beta[5]*sum((X[,3]^(4*alfa[2]))*log(X[,3]))+ beta[6]*sum((X[,2]^alfa[1])*(X[,3]^(3*alfa[2]))*log(X[,3])),
                beta[3]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,2]^(alfa[1]))*(X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,3])),  
                
                beta[2]*beta[3]*sum((X[,2]^(alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3])) + 2*beta[2]*beta[5]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])) +
                  beta[2]*beta[6]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3]))+ 2*beta[3]*beta[4]*sum((X[,2]^(2*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3])) + 4*beta[4]*beta[5]*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])) + 2*beta[4]*beta[6]*sum((X[,2]^(3*alfa[1]))*(X[,3]^alfa[2])*log(X[,2])*log(X[,3])) + beta[6]*beta[3]*sum((X[,2]^(alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])) + 2*beta[6]*beta[5]*sum((X[,2]^(alfa[1]))*(X[,3]^(3*alfa[2]))*log(X[,2])*log(X[,3])) + (beta[6]^2)*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*log(X[,2])*log(X[,3])),
                (beta[3]^2)*sum((X[,3]^(2*alfa[2]))*(log(X[,3])^2))+4*beta[3]*beta[5]*sum((X[,3]^(3*alfa[2]))*(log(X[,3])^2))+2*beta[3]*beta[6]*sum((X[,2]^alfa[1])*(X[,3]^(2*alfa[2]))*(log(X[,3])^2))+ 4*(beta[5]^2)*sum((X[,3]^(4*alfa[2]))*(log(X[,3])^2))+4*beta[5]*beta[6]*sum((X[,2]^alfa[1])*(X[,3]^(3*alfa[2]))*(log(X[,3])^2))+ (beta[6]^2)*sum((X[,2]^(2*alfa[1]))*(X[,3]^(2*alfa[2]))*(log(X[,3])^2))
    ),
    byrow = TRUE, nrow = 8, ncol = 8)
    }
    else if(alfa[1]==0 & alfa[2]!=0){
      M<-matrix(c(length(X[,1]),sum(log(X[,2])),sum(X[,3]^alfa[2]),sum(log(X[,2])^2),sum(X[,3]^(2*alfa[2])), sum(log(X[,2])*(X[,3]^alfa[2])), beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+ 2*beta[5]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[6]*sum(log(X[,2])*(X[,3]^alfa[2])*log(X[,3])), 
                  sum(log(X[,2])), sum(log(X[,2])^2), sum(log(X[,2])*(X[,3]^alfa[2])),sum(log(X[,2])^3), sum(log(X[,2])*(X[,3]^(2*alfa[2]))), sum((log(X[,2])^2)*(X[,3]^alfa[2])), beta[3]*sum(log(X[,2])*(X[,3]^alfa[2])*log(X[,3])) + 2*beta[5]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((log(X[,2]))^2*(X[,3]^alfa[2])*log(X[,3])),
                  sum(X[,3]^alfa[2]), sum(log(X[,2])*(X[,3]^alfa[2])), sum(X[,3]^(2*alfa[2])), sum((log(X[,2])^2)*(X[,3]^alfa[2])), sum(X[,3]^(3*alfa[2])), sum(log(X[,2])*(X[,3]^(2*alfa[2]))),beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3])),
                  sum(log(X[,2])^2),sum(log(X[,2])^3), sum(log(X[,2])^2*(X[,3]^alfa[2])), sum(log(X[,2])^4), sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))), sum((log(X[,2])^3)*(X[,3]^alfa[2])), beta[3]*sum((log(X[,2])^2)*(X[,3]^(alfa[2]))*log(X[,3]))+ 2*beta[5]*sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((log(X[,2])^3)*(X[,3]^alfa[2])*log(X[,3])),
                  sum(X[,3]^(2*alfa[2])), sum(log(X[,2])*(X[,3]^(2*alfa[2]))), sum(X[,3]^(3*alfa[2])), sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))), sum(X[,3]^(4*alfa[2])), sum(log(X[,2])*(X[,3]^(3*alfa[2]))),  beta[3]*sum((X[,3]^(3*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,3]^(4*alfa[2]))*log(X[,3])) + beta[6]*sum(log(X[,2])*(X[,3]^(3*alfa[2]))*log(X[,3])),
                  sum(log(X[,2])*(X[,3]^alfa[2])), sum(log(X[,2])^2*(X[,3]^alfa[2])), sum(log(X[,2])*(X[,3]^(2*alfa[2]))),sum(log(X[,2])^3*(X[,3]^alfa[2])), sum(log(X[,2])*(X[,3]^(3*alfa[2]))), sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))), beta[3]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum(log(X[,2])*(X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))*log(X[,3])),
                  
                  beta[3]*sum((X[,3]^alfa[2])*log(X[,3]))+ 2*beta[5]*sum((X[,3]^(2*alfa[2]))*log(X[,3]))+ beta[6]*sum(log(X[,2])*(X[,3]^alfa[2])*log(X[,3])), beta[3]*sum(log(X[,2])*(X[,3]^alfa[2])*log(X[,3])) + 2*beta[5]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((log(X[,2])^2)*(X[,3]^alfa[2])*log(X[,3])),  beta[3]*sum((X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum((X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3])), beta[3]*sum((log(X[,2])^2)*(X[,3]^(alfa[2]))*log(X[,3]))+ 2*beta[5]*sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))*log(X[,3])) + beta[6]*sum((log(X[,2])^3)*(X[,3]^alfa[2])*log(X[,3])),
                  beta[3]*sum((X[,3]^(3*alfa[2]))*log(X[,3]))+ 2*beta[5]*sum((X[,3]^(4*alfa[2]))*log(X[,3]))+ beta[6]*sum(log(X[,2])*(X[,3]^(3*alfa[2]))*log(X[,3])),
                  beta[3]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3])) + 2*beta[5]*sum(log(X[,2])*(X[,3]^(3*alfa[2]))*log(X[,3])) + beta[6]*sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))*log(X[,3])),  
                  
                  (beta[3]^2)*sum((X[,3]^(2*alfa[2]))*log(X[,3])^2)+4*beta[3]*beta[5]*sum((X[,3]^(3*alfa[2]))*log(X[,3])^2)+2*beta[3]*beta[6]*sum(log(X[,2])*(X[,3]^(2*alfa[2]))*log(X[,3])^2)+ 4*(beta[5]^2)*sum((X[,3]^(4*alfa[2]))*log(X[,3])^2)+4*beta[5]*beta[6]*sum(log(X[,2])*(X[,3]^(3*alfa[2]))*log(X[,3])^2)+ (beta[6])^2*sum((log(X[,2])^2)*(X[,3]^(2*alfa[2]))*log(X[,3])^2)
                                                                                                                                                  
                  ),
      byrow = TRUE, nrow = 7, ncol = 7)            
    }
    else if(alfa[1]!=0 & alfa[2]==0){
      M<-matrix(c(length(X[,1]),sum(X[,2]^alfa[1]),sum(log(X[,3])),sum(X[,2]^(2*alfa[1])),sum(log(X[,3])^2), sum((X[,2]^alfa[1])*log(X[,3])), beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^alfa[1])*log(X[,3])*log(X[,2])),  
                  sum(X[,2]^alfa[1]), sum(X[,2]^(2*alfa[1])), sum((X[,2]^alfa[1])*log(X[,3])),sum(X[,2]^(3*alfa[1])), sum((X[,2]^alfa[1])*log(X[,3])^2), sum((X[,2]^(2*alfa[1]))*log(X[,3])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + beta[6]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2])), 
                  sum(log(X[,3])), sum((X[,2]^alfa[1])*log(X[,3])), sum(log(X[,3])^2), sum(X[,2]^(2*alfa[1])*log(X[,3])), sum(log(X[,3])^3), sum((X[,2]^alfa[1])*log(X[,3])^2), beta[2]*sum((X[,2]^alfa[1])*log(X[,3])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2]))+ beta[6]*sum((X[,2]^(alfa[1]))*(log(X[,3])^2)*log(X[,2])), 
                  sum(X[,2]^(2*alfa[1])), sum(X[,2]^(3*alfa[1])), sum(X[,2]^(2*alfa[1])*log(X[,3])), sum(X[,2]^(4*alfa[1])), sum(X[,2]^(2*alfa[1])*log(X[,3])^2), sum(X[,2]^(3*alfa[1])*log(X[,3])), beta[2]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(4*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^(3*alfa[1]))*log(X[,3])*log(X[,2])), 
                  sum(log(X[,3])^2), sum((X[,2]^alfa[1])*(log(X[,3])^2)), sum(log(X[,3])^3), sum(X[,2]^(2*alfa[1])*log(X[,3])^2), sum(log(X[,3])^4), sum((X[,2]^alfa[1])*log(X[,3])^3), beta[2]*sum((X[,2]^alfa[1])*(log(X[,3])^2)*log(X[,2])) + 2*beta[4]*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])) + beta[6]*sum((X[,2]^(alfa[1]))*(log(X[,3])^3)*log(X[,2])), 
                  sum((X[,2]^alfa[1])*log(X[,3])), sum((X[,2]^(2*alfa[1]))*log(X[,3])), sum((X[,2]^alfa[1])*log(X[,3])^2),sum((X[,2]^(3*alfa[1]))*log(X[,3])), sum((X[,2]^alfa[1])*log(X[,3])^3), sum((X[,2]^(2*alfa[1]))*log(X[,3])^2), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*log(X[,3])*log(X[,2]))+ beta[6]*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])), 
                  
                  beta[2]*sum((X[,2]^alfa[1])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^alfa[1])*log(X[,3])*log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + beta[6]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2])),  beta[2]*sum((X[,2]^alfa[1])*log(X[,3])*log(X[,2]))+ 2*beta[4]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2]))+ beta[6]*sum((X[,2]^(alfa[1]))*log(X[,3])^2*log(X[,2])), beta[2]*sum((X[,2]^(3*alfa[1]))*log(X[,2])) + 2*beta[4]*sum((X[,2]^(4*alfa[1]))*log(X[,2]))+ beta[6]*sum((X[,2]^(3*alfa[1]))*log(X[,3])*log(X[,2])),
                  beta[2]*sum((X[,2]^alfa[1])*(log(X[,3])^2)*log(X[,2])) + 2*beta[4]*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])) + beta[6]*sum((X[,2]^(alfa[1]))*(log(X[,3])^3)*log(X[,2])), beta[2]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2])) + 2*beta[4]*sum((X[,2]^(3*alfa[1]))*log(X[,3])*log(X[,2]))+ beta[6]*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])), (beta[2]^2)*sum((X[,2]^(2*alfa[1]))*log(X[,2])^2)+4*beta[2]*beta[4]*sum((X[,2]^(3*alfa[1]))*log(X[,2])^2)+ 2*beta[2]*beta[6]*sum((X[,2]^(2*alfa[1]))*log(X[,3])*log(X[,2])^2)+ 4*(beta[4]^2)*sum((X[,2]^(4*alfa[1]))*log(X[,2])^2)+4*beta[4]*beta[6]*sum((X[,2]^(3*alfa[1]))*log(X[,3])*log(X[,2])^2)+ (beta[6]^2)*sum((X[,2]^(2*alfa[1]))*(log(X[,3])^2)*log(X[,2])^2)
                  
                  ),
      byrow = TRUE, nrow = 7, ncol = 7)              
    }
    else if(alfa[1]==0 & alfa[2]==0){
      M<-matrix(c(length(X[,1]),sum(log(X[,2])),sum(log(X[,3])),sum(log(X[,2])^2),sum(log(X[,3])^2), sum(log(X[,2])*(log(X[,3]))), beta[2]*sum((log(X[,2]))^2)+ 2*beta[4]*sum((log(X[,2]))^3)+ beta[6]*sum(((log(X[,2]))^2)*(log(X[,3]))), beta[3]*sum((log(X[,3]))^2)+ 2*beta[5]*sum((log(X[,3]))^3)+ beta[6]*sum((log(X[,2]))*(log(X[,3]))^2), 
                  sum(log(X[,2])), sum(log(X[,2])^2), sum((log(X[,2]))*(log(X[,3]))),sum(log(X[,2])^3), sum((log(X[,2]))*(log(X[,3])^2)), sum((log(X[,2])^2)*(log(X[,3]))), beta[2]*sum((log(X[,2])^3)) + 2*beta[4]*sum((log(X[,2]))^4) + beta[6]*sum((log(X[,2])^3)*(log(X[,3]))), beta[3]*sum((log(X[,2]))*(log(X[,3]))^2) + 2*beta[5]*sum((log(X[,2]))*(log(X[,3])^3)) + beta[6]*sum((log(X[,2])^2)*(log(X[,3])^2)),
                  sum(log(X[,3])), sum((log(X[,2]))*(log(X[,3]))), sum(log(X[,3])^2), sum((log(X[,2])^2)*(log(X[,3]))), sum(log(X[,3])^3), sum((log(X[,2]))*(log(X[,3])^2)), beta[2]*sum((log(X[,2])^2)*(log(X[,3])))+ 2*beta[4]*sum((log(X[,2])^3)*(log(X[,3])))+ beta[6]*sum((log(X[,2])^2)*(log(X[,3])^2)), beta[3]*sum(log(X[,3])^3) + 2*beta[5]*sum(log(X[,3])^4) + beta[6]*sum((log(X[,2]))*(log(X[,3])^3)),
                  sum(log(X[,2])^2), sum(log(X[,2])^3), sum((log(X[,2])^2)*(log(X[,3]))), sum(log(X[,2])^4), sum((log(X[,2])^2)*(log(X[,3])^2)), sum((log(X[,2])^3)*(log(X[,3]))), beta[2]*sum((log(X[,2])^4)) + 2*beta[4]*sum((log(X[,2])^5))+ beta[6]*sum((log(X[,2])^4)*log(X[,3])), beta[3]*sum((log(X[,2])^2)*(log(X[,3])^2))+ 2*beta[5]*sum((log(X[,2])^2)*(log(X[,3])^3)) + beta[6]*sum((log(X[,2])^3)*(log(X[,3])^2)),
                  sum(log(X[,3])^2), sum((log(X[,2]))*(log(X[,3])^2)), sum(log(X[,3])^3), sum((log(X[,2])^2)*(log(X[,3])^2)), sum((log(X[,3])^4)), sum((log(X[,2]))*(log(X[,3])^3)), beta[2]*sum((log(X[,2])^2)*(log(X[,3])^2)) + 2*beta[4]*sum((log(X[,2])^3)*(log(X[,3])^2)) + beta[6]*sum((log(X[,2])^2)*(log(X[,3])^3)), beta[3]*sum(log(X[,3])^4) + 2*beta[5]*sum(log(X[,3])^5) + beta[6]*sum((log(X[,2]))*(log(X[,3])^4)),
                  sum((log(X[,2]))*(log(X[,3]))), sum((log(X[,2])^2)*(log(X[,3]))), sum((log(X[,2]))*(log(X[,3])^2)),sum((log(X[,2])^3)*(log(X[,3]))), sum((log(X[,2]))*(log(X[,3])^3)), sum((log(X[,2])^2)*(log(X[,3])^2)), beta[2]*sum((log(X[,2])^3)*(log(X[,3]))) + 2*beta[4]*sum((log(X[,2])^4)*(log(X[,3])))+ beta[6]*sum((log(X[,2])^3)*(log(X[,3])^2)), beta[3]*sum((log(X[,2]))*(log(X[,3])^3)) + 2*beta[5]*sum((log(X[,2]))*(log(X[,3])^4)) + beta[6]*sum((log(X[,2])^2)*(log(X[,3])^3)),
                  
                  beta[2]*sum((log(X[,2]))^2)+ 2*beta[4]*sum((log(X[,2]))^3)+ beta[6]*sum(((log(X[,2]))^2)*(log(X[,3]))), beta[2]*sum((log(X[,2])^3)) + 2*beta[4]*sum((log(X[,2]))^4) + beta[6]*sum((log(X[,2])^3)*(log(X[,3]))),   beta[2]*sum((log(X[,2])^2)*(log(X[,3])))+ 2*beta[4]*sum((log(X[,2])^3)*(log(X[,3])))+ beta[6]*sum((log(X[,2])^2)*(log(X[,3])^2)),
                  beta[2]*sum((log(X[,2])^4)) + 2*beta[4]*sum((log(X[,2])^5))+ beta[6]*sum((log(X[,2])^4)*log(X[,3])), beta[2]*sum((log(X[,2])^2)*(log(X[,3])^2)) + 2*beta[4]*sum((log(X[,2])^3)*(log(X[,3])^2)) + beta[6]*sum((log(X[,2])^2)*(log(X[,3])^3)), beta[2]*sum((log(X[,2])^3)*(log(X[,3]))) + 2*beta[4]*sum((log(X[,2])^4)*(log(X[,3])))+ beta[6]*sum((log(X[,2])^3)*(log(X[,3])^2)),
                  (beta[2]^2)*sum((log(X[,2]))^4)+4*beta[2]*beta[4]*sum((log(X[,2]))^5)+2*beta[2]*beta[6]*sum((log(X[,2])^4)*log(X[,3]))+4*(beta[4]^2)*sum((log(X[,2]))^6)+4*beta[4]*beta[6]*sum((log(X[,2])^5)*log(X[,3]))+(beta[6]^2)*sum((log(X[,2])^4)*(log(X[,3])^2)),
                  beta[2]*beta[3]*sum((log(X[,2])^2)*(log(X[,3])^2)) + 2*beta[2]*beta[5]*sum((log(X[,2])^2)*(log(X[,3])^3)) +
                  beta[2]*beta[6]*sum((log(X[,2])^3)*(log(X[,3])^2))+ 2*beta[3]*beta[4]*sum((log(X[,2])^3)*(log(X[,3])^2)) + 4*beta[4]*beta[5]*sum((log(X[,2])^3)*(log(X[,3])^3)) + 2*beta[4]*beta[6]*sum((log(X[,2])^4)*(log(X[,3])^2)) + beta[6]*beta[3]*sum((log(X[,2])^2)*(log(X[,3])^3)) + 2*beta[6]*beta[5]*sum((log(X[,2])^2)*(log(X[,3])^4)) + (beta[6]^2)*sum((log(X[,2])^3)*(log(X[,3])^3)),  
                  
                  
                  beta[3]*sum((log(X[,3]))^2)+ 2*beta[5]*sum((log(X[,3]))^3)+ beta[6]*sum((log(X[,2]))*(log(X[,3]))^2), beta[3]*sum((log(X[,2]))*(log(X[,3]))^2) + 2*beta[5]*sum((log(X[,2]))*(log(X[,3])^3)) + beta[6]*sum((log(X[,2])^2)*(log(X[,3])^2)), beta[3]*sum(log(X[,3])^3) + 2*beta[5]*sum(log(X[,3])^4) + beta[6]*sum((log(X[,2]))*(log(X[,3])^3)), beta[3]*sum((log(X[,2])^2)*(log(X[,3])^2))+ 2*beta[5]*sum((log(X[,2])^2)*(log(X[,3])^3)) + beta[6]*sum((log(X[,2])^3)*(log(X[,3])^2)),
                  beta[3]*sum(log(X[,3])^4) + 2*beta[5]*sum(log(X[,3])^5) + beta[6]*sum((log(X[,2]))*(log(X[,3])^4)), beta[3]*sum((log(X[,2]))*(log(X[,3])^3)) + 2*beta[5]*sum((log(X[,2]))*(log(X[,3])^4)) + beta[6]*sum((log(X[,2])^2)*(log(X[,3])^3)),
                  beta[2]*beta[3]*sum((log(X[,2])^2)*(log(X[,3])^2)) + 2*beta[2]*beta[5]*sum((log(X[,2])^2)*(log(X[,3])^3)) +
                  beta[2]*beta[6]*sum((log(X[,2])^3)*(log(X[,3])^2))+ 2*beta[3]*beta[4]*sum((log(X[,2])^3)*(log(X[,3]))^2) + 4*beta[4]*beta[5]*sum((log(X[,2])^3)*(log(X[,3])^3)) + 2*beta[4]*beta[6]*sum((log(X[,2])^4)*(log(X[,3])^2)) + beta[6]*beta[3]*sum((log(X[,2])^2)*(log(X[,3])^3)) + 2*beta[6]*beta[5]*sum((log(X[,2])^2)*(log(X[,3])^4)) + (beta[6]^2)*sum((log(X[,2])^3)*(log(X[,3])^3)),  
                  (beta[3]^2)*sum(log(X[,3])^4)+4*beta[3]*beta[5]*sum(log(X[,3])^5)+2*beta[3]*beta[6]*sum((log(X[,2]))*(log(X[,3])^4))+ 4*(beta[5]^2)*sum((log(X[,3])^6))+4*beta[5]*beta[6]*sum((log(X[,2]))*(log(X[,3])^5))+ (beta[6]^2)*sum((log(X[,2])^2)*(log(X[,3])^4))
      ),
      byrow = TRUE, nrow = 8, ncol = 8)  
      
      
      
                                                                        }
  }
  
    
  list(M=M)
}

###############################################
# FUNCTION FOR GENERATING INITIAL DESIGN (POINT EXCHANGE)
dinicial <- function(comb,Ncand,N,pinv, m,K, beta, alfa)
{
  ind <- matrix(sample(seq(1:Ncand),N,replace=TRUE),N,1)
  ind <- as.matrix(ind[order(ind)])
  di <- comb[ind,]
  M_matrix <- fisher_matrix(di[,-1],m,K, beta, alfa) 
  M <- M_matrix$M
  d <- det(M)
  list(di=di,d=d,M=M)
}

# CONTROLS FOR NO-SINGULAR INITIAL DESIGN (POINT EXCHANGE)
SampleD <- function(comb,Ncand,N,pinv,m,K, beta, alfa)
{
  d <- 0
  while(d<10^(-6))
  {
    dini <- dinicial(comb,Ncand,N,pinv,m,K, beta, alfa)
    d <- dini$d
  }
  X <- as.matrix(dini$di)
  M <- dini$M
  list(X=X,d=d,M=M)
}
# CONTROLS FOR NON-SINGULAR INITIAL DESIGN (COORDINATE EXCHANGE)
SampleDcoord <- function(K,N,Nlev,Terms,pinv, X,m, beta, alfa)  
{
  d <- 0
  while(d<10^(-6))
  {
    DC<-DesigniCoord(N,Nlev,Terms,K, alfa)
    X <- DC$di
    Xagg<-DC$dt
    M <- fisher_matrix(Xagg,m,K, beta, alfa)$M 
    d <- det(M)
  }
  list(X=X, M=M, d=d, Xagg=Xagg)
}


# FUNCTION FOR EXPANDING THE X MATRIX 

ModelMatrix <- function(x,Terms,K, alfa)
{
  Nr <- nrow(x)
  X <- matrix(c(rep(1,Nr)),Nr,1)  
  contador<-0
  if(K ==1 & m==1){
    if(Terms[1]>0){X <- cbind(X,x[,1])^(alfa[1])}
  } 
  else if(K ==1 & m==2){
    if(Terms[1]>0){X <- cbind(X,x[,1])^(alfa[1])}
    if(Terms[1]>0){X <- cbind(X,(x[,1])^(2*alfa[1]))}
  }
  else if(K ==2 & m==1){
    if(Terms[1]>0){X <- cbind(X,x[,1]^(alfa[1]))}
    if(Terms[1]>0){X <- cbind(X,(x[,2])^(alfa[2]))}
    if(Terms[1]>0){X <- cbind(X,(x[,1]^(alfa[1]))*(x[,2]^(alfa[2])))}
  }
  else if(K ==2 & m==2){
    if(Terms[1]>0){X <- cbind(X,(x[,1])^(alfa[1]))}
    if(Terms[1]>0){X <- cbind(X,(x[,2])^(alfa[2]))}
    if(Terms[1]>0){X <- cbind(X,(x[,1])^(2*alfa[1]))}
    if(Terms[1]>0){X <- cbind(X,(x[,2])^(2*alfa[2]))}
    if(Terms[1]>0){X <- cbind(X,(x[,1]^(alfa[1]))*(x[,2]^(alfa[2])))}
  }
  
  
  X <- matrix(X,nrow=Nr)
  list(X=X)
}

# FUNCTION TO EXCHANGE ROWS (POINT EXCHANGE)
swap <- function(cand,D,M,V,Det,crita,W,Ncand,N,pinv,Ncali, m,K, beta, alfa) #linha 7 a 17
{
  improve <- 0
  calibrate <- -1
  for(i in 1:Ncand)
  {
    if(D[1,1]!=cand[i,1])
    {  
      Xc <- D
      xin <- matrix(cand[i,-1],nc=1)
      xout <- matrix(Xc[1,-1],nc=1)
      Xc[1,] <- cand[i,]
      critc <- 0
      ns <- 0
      calibrate <- calibrate+1
      if(calibrate>=Ncali)
      {
        
        Mup <- fisher_matrix(Xc,m,K, beta, alfa)$M
        Detup<-det(Mup)
        calibrate <- -1
        if(Detup > 10^(-2)) 
        {
          ns <- 1
          Vup <- 0
        }
      }
      
      if(ns>0)
      {
        critD <- Detup
        
        critc <- critD
      }
      if(critc>crita)
      {
        D <- Xc
        crita <- critc
        busca <- 1
        M <- Mup
        V <- Vup
        Det <- Detup
      } 
    }
  }
  calibrate <- -1
  for(l in 2:N) 
  {
    if(D[l,1]!=D[(l-1),1])
    {  
      Xc <- D 
      for(i in 1:Ncand)
      {
        if(D[l,1]!=cand[i,1])
        {
          Xc <- D
          xin <- matrix(cand[i,-1],nc=1)
          xout <- matrix(Xc[l,-1],nc=1)
          Xc[l,] <- cand[i,]
          critc <- 0
          calibrate <- calibrate+1 
          if(calibrate>=Ncali)
          {
            Mup <- fisher_matrix(Xc,m,K, beta, alfa)$M
            
            Detup<-log(det(Mup))
            calibrate <- -1
            if(Detup > 10^(-2)) 
            {
              ns <- 1
              Vup <- 0
            }
          }
          
          
          if(ns>0)
          {
            critD <- Detup
            
            critc <- critD
          }
          if(critc>crita)
          {
            D <- Xc
            crita <- critc
            improve <- 1
            Det <- Detup
            M <- Mup
            V <- Vup
          } 
        }
      }
    }
  }
  list(D=D,crita=crita,improve=improve,Det=Det,M=M,V=V)
}

# FUNCTION FOR CODING THE LEVELS (POINT EXCHANGE)
codifica <- function(x)
{
  xc <- (x-(max(x)+min(x))/2)/((max(x)-min(x))/2)
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

# FUNCTION TO DRIVE THE SEARCH 
SearchTreat <- function(K,Levels,Cubic,N,Ntries, m, beta, alfa,
                        Terms,Npar,AlgType)

{
  start <- Sys.time()
  Power <- 1/Npar
  critall <- matrix(0,nr=Ntries)
  Xall<- vector("list",Ntries)
  # POINT EXCHANGE
  if(AlgType=="P")
  {
    cand <-as.matrix(expand.grid(Levels)) # matriz candidata
    candc <- apply(cand,2,codifica)
    if(Cubic=='N') candc <- as.matrix(Sphcand(candc)$candc) 
    candc <- TreatLabels(candc)$Treat
    cand <- cbind(candc[,1],ModelMatrix(candc[,-1],Terms,K, alfa)$X)
    Ncand <- nrow(cand)
    for(nt in 1:Ntries) 
    {
      Des <- SampleD(cand,Ncand,N,Power,m,K, beta, alfa)
      X <- Des$X
      M <- Des$M
      Det <- Des$d
      crite <- 0
      criteD <- 0
      
      if(Det>10^(-6))
      {
        V <- 0
        criteD <- Det
      }
      
      crite <- criteD
      improve<-1
      while(improve==1) 
      {          
        Xs <- swap(cand,X,M,V,Det,crite,W,Ncand,N,Power,Ncali,m,K, beta, alfa)
        X <- Xs$D
        crite <- Xs$crita
        improve <- Xs$improve
        M <- Xs$M
        V <- Xs$V
        Det <- Xs$Det
      }
      critall[nt] <- crite
      if(nt==1)
      {
        critopt <- crite
        Xopt <- X
        Vopt <- V
      }else
      {
        if(crite > critopt)
        {
          Xopt <- X
          critopt <- crite
          Vopt <- V
        }
      }
    }
  }
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
      Des <- SampleDcoord(K,N,Nlev,Terms,Power, X,m, beta, alfa) 
      X <- Des$X
      Xagg<-Des$Xagg 
      M <- Des$M
      Det <- Des$d
      crite <- 0
      criteD <- 0
      if(Det>10^(-6))
      {
        
        criteD <- Det
        
      }
      
      crite <- criteD
      improve <- 1 
      iteration<-0
      while(improve==1)
      {         
        Xs <- CoordRow(Xagg,M,Det,iteration,Xall[nt],crite,K,N,Nlev,Terms,Npar,Power, m, beta, alfa)
        X <- Xs$Xmod_orig 
        Xj<-Xs$D 
        crite <- Xs$crita
        improve <- Xs$improve
        M <- Xs$M
        Det <- Xs$Det
        iteration<-Xs$iter
        Xall[[nt]]<-Xs$Xall
        
      }
      critall[nt] <- crite
       
      if(nt==1)  
      {
        critopt <- crite
        Xopt <- X
        Xopt_ref<-Xj
        
      }else
      {
        if(crite > critopt)
        {
          Xopt <- X
          Xopt_ref<-Xj
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
      Des <- SampleDcoord(K,N,Nlev,Terms,Power, X,m, beta, alfa) 
      X <- Des$X
      Xagg<-Des$Xagg
      M <- Des$M
      Det <- Des$d
      crite <- 0
      criteD <- 0
        
      if(Det>10^(-6))
      {
     
        criteD <-Det
      }
      
      crite <- criteD
      improve <- 1 
      iteration<-0
      while(improve==1)
      {         
        Xs <- CoordCol(Xagg,M,Det,iteration,Xall[nt],crite,K,N,Nlev,Terms,Npar,Power, m, beta, alfa) 
        X <- Xs$Xmod_orig 
        Xj <- Xs$D
        crite <- Xs$crita
        improve <- Xs$improve
        M <- Xs$M
        Det <- Xs$Det
        iteration<-Xs$iter
        Xall[[nt]]<-Xs$Xall
      }
      critall[nt] <- crite
      if(nt==1)
      {
        critopt <- crite
        Xopt <- X
        Xopt_ref<-Xj
        
      }else
      {
        if(crite > critopt)
        {
          Xopt <- X
          Xopt_ref<-Xj
          critopt <- crite
          
        }
      }
    }
    Xopt <- TreatLabels(Xopt)$Treat
  }
  
  tend <- Sys.time()-start
  
  M_matrix <- fisher_matrix(Xopt_ref,m,K, beta, alfa)
  M <- M_matrix$M  
 
  D<-det(M)
  

  list(Xopt=Xopt, Xopt_ref=Xopt_ref,critopt=critopt,D=D,
      critall=critall,start=start,tend=tend,M=M,Xall=Xall,iter=iteration) 
}


# FUNCTION TO EXCHANGE COORDINATES BY ROW 
CoordRow <- function(D,M,Det,iteration,Xall,crita,K,N,Nlev,Terms,Npar,pinv, m, beta, alfa) 
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
          xin_mod <- matrix((ModelMatrix(xin,Terms,K, alfa)$X)[,-1],nc=1)
          xin_tot<-matrix(c(1,xin[1:K],xin_mod),nc=1)
          Xmod_orig<-cbind(D[,1],D[,((K+2):ncol(D))])
          Xc <- D
          Xc[i,] <- t(xin_tot)
          critc <- 0
          ns <- 0
          Xmod<-cbind(Xc[,1],Xc[,((K+2):ncol(Xc))])
          
          Mup<- fisher_matrix(Xc,m,K, beta, alfa)$M 
          
          Detup<-det(Mup)
          
          
          if(Detup > 10^(-6)) 
          {
            ns <- 1
            Detup<-det(Mup)
            
          }
        
          if(ns>0)
          {
            
            critD <- Detup
            critc <- critD
            
          }
          if(critc>crita)
          {
            D <- Xc 
            Xmod_orig<-Xmod
            crita <- critc
            busca <- 1
            M <- Mup
            Det <- Detup
            Improve <- 1
            iteration<-iteration+1
      
            Xall<-add_element(Xall,list(iteration = iteration, criterio = crita)) 
          } 
        }
      }
    }
  }
  list(D=D,crita=crita,improve=improve,M=M,Det=Det, Xmod_orig=Xmod_orig, Xall=Xall)
}

# FUNCTION TO EXCHANGE COORDINATES BY COLUMN         
CoordCol<- function(D,M,Det,iteration,Xall,crita,K,N,Nlev,Terms,Npar,pinv,m, beta, alfa)
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
          xin_mod <- matrix(ModelMatrix(xin,Terms,K, alfa)$X[,-1],nc=1)
          xin_tot<-matrix(c(1,xin[1:K],xin_mod),nc=1)
          Xmod_orig<-cbind(D[,1],D[,((K+2):ncol(D))])
          Xc <- D
          Xc[i,] <- t(xin_tot)
          critc <- 0
          ns <- 0
          Xmod<-cbind(Xc[,1],Xc[,((K+2):ncol(Xc))])
          
          Mup<- fisher_matrix(Xc,m,K, beta, alfa)$M
           
          Detup<-det(Mup)
            
          if(Detup > 10^(-6)) 
          {
            ns <- 1
          }
          
          if(ns>0)
          {
           
            critD <- Detup
            
            critc <- critD
          }
          if(critc>crita)
          {
            D <- Xc
            Xmod_orig<-Xmod
            crita <- critc
            busca <- 1
            M <- Mup
            Det <- Detup
            Improve <- 1
            iteration<-iteration+1
            Xall<-add_element(Xall,list(iteration = iteration, criterio = crita))
          } 
        }
      }
    }
  }
  list(D=D,crita=crita,improve=improve,M=M,Det=Det, Xmod_orig=Xmod_orig, Xall=Xall)
}




