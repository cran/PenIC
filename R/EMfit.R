
# =========================================================================
# *************************************************************************
# Fit the model using Isotonic Regression with EM algorithm
EM_fit <- function(g0,b0,d1,d2,d3,Li,Ri,Z,nsub,alpha,qn,order,t.seq,tol=1e-5,itmax=500,lamu=1e5){

  gf <- function(x,alpha){ #alpha>0
		log(((1-x)^(-alpha)-1)/alpha)
	}

	gfinv <- function(x,alpha){ #alpha>0
		1-(1+alpha*exp(x))^(-1/alpha)
	}

	#######################################################################
	# =====================================================================
	# Functions provided:
	# Log Likelihood function with penalty
	#
	# =====================================================================
	#######################################################################
	llh_em <- function(theta,alpha,lambda){
		gamma <- tail(theta,K_0)
		if(alpha==0){ # PH model
			a     <- as.numeric(exp(ldmat%*%theta))
			b     <- as.numeric(exp(rdmat%*%theta))
			res1  <- sum(log(1-exp(-b[d1==1])))
			cof2  <- exp(-a[d2==1])-exp(-b[d2==1])
			if(sum(is.na(cof2))>=1){print(theta)}
			res2  <- sum(log(cof2))
			res3  <- sum(-a[d3==1])
		}else{        # any alpha value
			HL    = as.numeric(-log(1-gfinv(ldmat%*%theta,alpha)))
			HR    = as.numeric(-log(1-gfinv(rdmat%*%theta,alpha)))
			res1  <- sum(log(1-exp(-HR[d1==1])))
			cof2  <- exp(-HL[d2==1])-exp(-HR[d2==1])
			if(sum(is.na(cof2))>=1){print(theta)}
			res2  <- sum(log(cof2))
			res3  <- sum(-HL[d3==1])
		}

		res   <- res1+res2+res3-lambda*t(gamma)%*%DMD%*%gamma/2
		return(-as.numeric(res))
	}

	#######################################################################
	# =====================================================================
	# EM functions
	# Using B-splines to estimate
	#
	# =====================================================================
	#######################################################################

	BQfun <- function(theta,theta_p,alpha,lambda){
		gamma <- tail(theta,K_0)

		if(alpha==0){ # PH model
			HR    <- as.numeric(exp(rdmat%*%theta_p))
			HL    <- as.numeric(exp(ldmat%*%theta_p))
			llam1 <- as.numeric(ldmat%*%theta)
			llam2 <- as.numeric(rdmat%*%theta)
			elam1 <- exp(llam1)
			elam2 <- exp(llam2)
			lam1  <- elam1
			lam2  <- elam2-elam1
			EY    <- rep(0,length(d1))
			EW    <- rep(0,length(d1))
			EY[d1==1] <- HR[d1==1]/(1-exp(-HR[d1==1]))
			EWv   <- (HR[d2==1]-HL[d2==1])/(1-exp(HL[d2==1]-HR[d2==1]))
			#idd   <- abs(HR[d2==1]-HL[d2==1]) <1e-6
			#EWv[idd] <- 1
			EW[d2==1] <- EWv
			EY    <- EY*d1
			EW    <- EW*d2
			res1  <- sum(EY*llam1-lam1)
			res2  <- sum((EW*log(lam2)-lam2)[d2==1])
		}else{
			HR    <- as.numeric(-log(1-gfinv(rdmat%*%theta_p,alpha)))
			HL    <- as.numeric(-log(1-gfinv(ldmat%*%theta_p,alpha)))
			EY    <- rep(0,length(d1))
			EW    <- rep(0,length(d1))
			EY[d1==1] <- HR[d1==1]/(1-exp(-HR[d1==1]))
			EWv   <- (HR[d2==1]-HL[d2==1])/(1-exp(HL[d2==1]-HR[d2==1]))
			#idd   <- abs(HR[d2==1]-HL[d2==1]) <1e-6
			#EWv[idd] <- 1
			EW[d2==1] <- EWv
			EY    <- EY*d1
			EW    <- EW*d2
			lam1  <- as.numeric(-log(1-gfinv(ldmat%*%theta,alpha)))
			lam2  <- as.numeric(-log(1-gfinv(rdmat%*%theta,alpha))) - lam1
			par1  <- EY*log(lam1)
			res1  <- sum(par1[d1==1])- sum(lam1)
			res2  <- sum((EW*log(lam2)-lam2)[d2==1])
		}

		res   <- res1 + res2 - lambda*t(gamma)%*%DMD%*%gamma/2
		return(-res)
	}

	EMoptim <- function(par,alpha,lambda,tol=1e-3,ctmax=50){
		environment(BQfun)   <- environment()
		environment(llh_em)     <- environment()
		#environment(llh_PO)     <- environment()

		cons_mat <- cbind(matrix(0,K_0-1,p),diff(diag(K_0)))
		cons_vec <- rep(0,K_0-1)

		theta0   <- par
		abs.diff <- 1
		count    <- 1
		while(count<ctmax&abs.diff>tol){

			fit      <- constrOptim(theta0,f=BQfun,method="Nelder-Mead",ui=cons_mat,ci=cons_vec,
										hessian=FALSE,theta_p=theta0,alpha=alpha,lambda=lambda)
			#print(fit$convergence)

			theta    <- fit$par
			abs.diff <- max(abs(theta-theta0))
			theta0   <- theta
			count     <- count + 1
		}

		hesmat <- hessian(llh_em,theta,alpha=alpha,lambda=lambda)

		##################################################################################
		return(list(par=theta0,hessian=hesmat,convergence=as.numeric(count>ctmax)))
	}

	#------------Matrices for left and right------------;
	Li[d1 == 1] <- Ri[d1 == 1]
	Ri[d3 == 1] <- Li[d3 == 1]

	# Prepare for model fitting
  Y        <- c(Li[d1 == 0], Ri[d3 == 0])
	prob     <- seq(0,1,1/(qn+1))
	psi      <- as.numeric(unique(quantile(Y,prob)))
	knots    <- c(rep(min(Y),3),psi,rep(max(Y),3))
	lmatrix  <- splineDesign(knots = knots,Li,ord = order+1,outer.ok = TRUE)
	rmatrix  <- splineDesign(knots = knots,Ri,ord = order+1,outer.ok = TRUE)
	bt       <- splineDesign(knots = knots,t.seq,ord = order+1,outer.ok = TRUE)
	ldmat    <- cbind(Z,lmatrix)
	rdmat    <- cbind(Z,rmatrix)
	p        <- ncol(Z)
	K_0      <- ncol(lmatrix)

	D2       <- diff(diag(K_0),difference = 2)
	DMD      <- crossprod(D2)

	environment(EMoptim)    <- environment()

	#------------Set initial values-------------;
	#theta0   <- c(rep(0,p), sort(runif(K_0,-1,1)))
	theta0   <- c(b0, g0)

	#------------Conduct the optimization to find theta with given lambda------------
	#------------Unpenalized spline method------------;
	lmbd0    <- 0

	out0     <- EMoptim(theta0,alpha,tol=tol*10,lambda=lmbd0)

	#------------Penalized spline method------------;
	diff     <- 1
	par.n    <- theta0
	lambda.n <- 0.1
	S        <- as.matrix(bdiag(matrix(0,p,p),DMD))
	S.inv    <- ginv(S)
	fit0     <- EMoptim(par.n,alpha=alpha,tol=tol*10,lambda=lambda.n)
	hess     <- fit0$hessian

	iter     <- 0
	while(diff>tol&iter<=itmax){
		#print(par.n)

		par.c       <- par.n
		lambda.c    <- lambda.n
    quad        <- t(par.c)%*%S%*%par.c
		hess.inv    <- solve(hess)
		tr.2        <- sum(diag(hess.inv%*%S))
		tr.1        <- sum(diag(ginv(S*lambda.c)%*%S))
		lambda.n    <- as.numeric((tr.1-tr.2)/quad)*lambda.c
		if(lambda.n > lamu){lambda.n = lamu}
		fit1        <- EMoptim(par.c,alpha=alpha,tol=tol*10,lambda=lambda.n)
    par.n       <- fit1$par
		hess        <- fit1$hessian
		diff        <- max(abs(par.n - par.c))
		#print(diff)
		iter        <- iter + 1
	}

	#------Save results ------;

	theta1  <- par.n
	beta1   <- theta1[1:p]
	gama1   <- tail(theta1,K_0)
  se1     <- sqrt(diag(solve(hess)))[1:p]
	base1   <- as.numeric(bt%*%gama1)

	return(list(b=beta1,g=gama1,se=se1,theta=theta1,base=base1,
				lambda=lambda.n,flag=as.numeric(iter>itmax)))

}

