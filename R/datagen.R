
####################################################################### 
# =====================================================================
# DATA generation functions
#
#
# =====================================================================
####################################################################### 

dataPA <- function(N,case = 1,alpha){
	d1   <- rep(0,N)
	d2   <- rep(0,N)
	d3   <- rep(0,N)
	Ri   <- rep(Inf,N)
	Li   <- rep(0,N)


	if(case == 1){#d3 = %
		base  <- function(t){t^2/5+t/5}
		ppois <- 1
		pexp  <- 0.5
	}
	if(case == 2){#d3 = %
		base  <- function(t){t} 
		ppois <- 3
		pexp  <- 0.5
	}
	if(case == 3){#d3 = %
		base  <- function(t){log(1+3*t)+t/3}
		ppois <- 3
		pexp  <- 0.5
	}
	
	b  <- matrix(c(-1,-1,-1,1,1,-1),3,2,byrow=TRUE)[case,]
	Z  <- cbind(rbinom(N,1,.5),rnorm(N,0,1))
	xb <- Z%*%b   # Generate under partially linear model

	# Assume the survival time data modeled by transformation model with
	# the cumulative baseline hazard function ``base''
	# Generate survival time data

	fail  <- function(t){
		res <- base(t)+temp
		return(res)
	}

	ft <- rep(-99,N)
	for(i in 1:N){
		u     <- runif(1)
		if(alpha==0){
			temp  <- log(1-u)/exp(xb[i]) # PH
		}else{
			temp  <- (-1)/alpha*((1-u)^{-alpha}-1)*exp(-xb[i]) # PAny
		}
		
		#temp  <- (-1)*u/(1-u)*exp(-xb[i]) # PO
		#temp  <- (-2)*((1-u)^{-1/2}-1)*exp(-xb[i]) # PM
		
		ft[i] <-uniroot(fail,c(.000000000000001,400000000))$root
	}

	for(i in 1:N){
		obs   <- c(0,cumsum(rexp(1+rpois(1,ppois),1/pexp)),Inf)
		Li[i] <- max(obs[obs<ft[i]])
		Ri[i] <- min(obs[obs>=ft[i]])      
		if(Li[i]==0)           {d1[i]<-1}
		if(Ri[i]==Inf)         {d3[i]<-1}
		if(Li[i]>0 & Ri[i]<Inf){d2[i]<-1}
	}

	return(list("d1"=d1,"d2"=d2,"d3"=d3,"Li"=Li,"Ri"=Ri,"Z"=Z))

}

