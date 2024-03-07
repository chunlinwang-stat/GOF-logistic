####################################################
## Calculate the p_i & q_i under H0: 
## logistic model or exponential tilting DRM 
####################################################
est.H0 <- function(x0,x1){
        n0 <- length(x0)
        n1 <- length(x1)
        y <- c(rep(0,n0),rep(1,n1))
        y <- factor(y)
        n <- n0+n1
        x <- c(x0,x1)      
        dat <- data.frame(y=y, x=x)
        dat <- dat[order(x),]
        x <- sort(x) 
        res <- glm(dat[,1]~dat[,2], family=binomial(link="logit"))
        coeff <- coef(res)
        beta <- coeff[2]
        alpha <- coeff[1]+log(n0/n1) 
        pi <- 1/(n0+n1*exp(alpha+beta*x))+ 1e-9
        qi <- exp(alpha+beta*x)*pi+ 1e-9
        return(list(dat=x, dG=pi, dH=qi))
}

####################################################
## Calculate the p_i & q_i under H1: monotonic DRM
####################################################
est.H1 <- function(x0,x1){      
        n0 <- length(x0)
        n1 <- length(x1)
        y <- c(rep(0,n0),rep(1,n1))
        n <- n0+n1
        x <- c(x0,x1)       
        lambda <- n1/n
        ri <- y[order(x)]
        di <- 1-ri
        x <- sort(x)
        PAVAy <- ri/(di+ri)
        PAVAw <- di+ri
        phi.hat <- (di+ri)/n
        
        ## theta(x) is increasing;
        theta.hat.inc <- Iso::pava(PAVAy,PAVAw)   
        qi.inc <- phi.hat*theta.hat.inc/lambda + 1e-9
        pi.inc <- phi.hat*(1-theta.hat.inc)/(1-lambda) + 1e-9
        lik.inc <- sum(di*log(pi.inc)+ri*log(qi.inc))        
        pi.inc <- pi.inc[order(x)]
        qi.inc <- qi.inc[order(x)]
        return(list(dG.inc=pi.inc, dH.inc=qi.inc, lik=lik.inc))
}

####################################################
# This function is to calculate three KS Statistics:
# (1) Qin & Zhang (1997) with ECDF under H1; 
# (2) monotonic increasing (with bete>0) DRM constraint under H1; 
# and (3) general monotonic DRM constraint under H1;
####################################################

KS.all <- function(x0,x1){
        xx <- c(x0,x1)
        x <- sort(xx)
        n <- length(x)
        n0 <- length(x0)
        n1 <- length(x1)
        #### empirical distribution Gn(x);
        Gn <- function(z){ (1/n0)*sum(x0<=z) }        
        #### MLE of G(x) under exponential tilting DRM;
        g.drm <- est.H0(x0,x1)$dG
        G.drm <- function(z){ sum(g.drm*(x<=z)) }      
        #### MLE of G(x) under monotonic DRM;
        res.H1.inc <- est.H1(x0,x1)
        dG.inc <- res.H1.inc$dG.inc
        Ghat.inc <- function(z){ sum(dG.inc*(x<=z))}   
        #### with no info on beta sign; 
        res.H1.dec <- est.H1(x1,x0)
        dG.dec <- res.H1.dec$dH.inc   
        if(res.H1.inc$lik > res.H1.dec$lik){dG.mono <- dG.inc} else {dG.mono <- dG.dec}
        Ghat.mono <- function(z){ sum(dG.mono*(x<=z))}
        
        #### KS Statistic of Qin & Zhang (1997);
        Delta.i <- function(z){ sqrt(n)*abs(Gn(z)-G.drm(z)) }
        Delta <- sapply(x, Delta.i)
        KS.QZ <- max(Delta)
        #### Proposed KS Statistic;
        Delta.i.inc <- function(z){ sqrt(n)*abs(Ghat.inc(z)-G.drm(z)) }
        Delta.i.mono <- function(z){ sqrt(n)*abs(Ghat.mono(z)-G.drm(z)) }
        Delta.inc <- sapply(x, Delta.i.inc)
        Delta.mono <- sapply(x, Delta.i.mono)
      
        KS.inc <- max(Delta.inc)
        KS.mono <- max(Delta.mono)
        KS.all <- c(KS.QZ, KS.inc, KS.mono)
        names(KS.all) <- c("KS.QZ", "KS.inc", "KS.mono")
        return(list(KS.all=KS.all))
}

####################################################
## This is the main function; using bootstrap approximation of the p-values;
####################################################
KS.bootstrap <- function(x0,x1, B=999){   
        obs.KS.all <- KS.all(x0,x1)$KS.all     
        res.H0 <- est.H0(x0,x1)
        g0 <- res.H0$dG
        g1 <- res.H0$dH
        x <- res.H0$dat
        n0 <- length(x0)
        n1 <- length(x1)
        n  <- n0+n1 
        x0star <- matrix(sample(x, replace = T, prob = g0, size = n0*B), B, n0)
        x1star <- matrix(sample(x, replace = T, prob = g1, size = n1*B), B, n1)
        x.star <- cbind(x0star, x1star)
        KS.test.stat <- function(x){ xx0<-x[1:n0]; xx1 <- x[(n0+1):n]; KS.all(xx0,xx1)$KS.all }
        Delta.star <- apply(x.star, 1, KS.test.stat)
        pval.all <- apply(cbind(obs.KS.all, Delta.star), 1, function(y){ mean(y[-1] > y[1])} ) 
        return(list(obs.KS.all=obs.KS.all, pvalue.all=pval.all))  
}

####################################################
## Examples for illustration;
####################################################
## the null case;
# set.seed(202403)
# control <- rnorm(50,0,1)
# case <- rnorm(50,2,1)
# KS.bootstrap(control, case)

## the alternative case;
# set.seed(202403)
# n0.mix <- rbinom(1, 50, 1-0.9)
# n1.mix <- rbinom(1, 50, 1-0.7)
# control <- c(rnorm(50-n0.mix,0,1), rnorm(n0.mix,2,1))
# case <- c(rnorm(n1.mix,0,1), rnorm(50-n1.mix,2,1))
# KS.bootstrap(control, case)
