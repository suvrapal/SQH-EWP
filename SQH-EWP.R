library(numDeriv)

# data generation function

sim_data_Wei=function(n,g1,g2,delta,eta1,p_low,p_high,phi){
  
  t=rep(NA,n)
  d=rep(NA,n)
  z=rep(NA,n)#eta
  x=rep(NA,n)#tumor thickness
  u=runif(n,0,1)
  
  for(i in 1 :n){
    if(u[i] <= 0.44){
      z[i]=1
      x[i]=rweibull(1,shape=1.36,scale=4.74)
    }else{
      z[i]=0
      x[i]=rexp(1,rate=0.55)
    }
  }
  
  x_min=min(x)
  x_max=max(x)
  
  b21=log(eta1)
  
  b11=(log(p_high/(1-p_high))-log(p_low/(1-p_low)))/(x_max-x_min)
  b10=log(p_low/(1-p_low))-(b11*x_min)
  
  #print("b21")
  #print(b21)
  #print("b10")
  #print(b10)
  #print("b11")
  #print(b11)
  
  eta=rep(NA,n)
  p=rep(NA,n)
  D=rep(NA,n)
  C=rexp(n,rate=delta)
  count=0
  
  
  for(i in 1:n){
    eta[i]=exp((b21*z[i]))
    p[i]=exp(b10+(b11*x[i]))/(1+exp(b10+(b11*x[i])))
    
    m=rpois(1,lambda=(eta[i]*exp(phi)))
    if(m==0){
      D[i]=0
    }else{
      D[i]=rbinom(1,size=m,prob=p[i])
    }
    
    #m=rpois(1,lambda=(eta[i]*p[i]))
    #D[i]=rpois(1,lambda=(eta[i]*p[i]*exp(phi)))
    
    if(D[i]==0){
      count=count+1
      t[i]=C[i]
      d[i]=0
    }else{
      y=min(rweibull(D[i],shape=1/g1,scale=1/g2))
      t[i]=min(y,C[i])
      if(min(y,C[i])==C[i]){
        d[i]=0
      }else{
        d[i]=1
      }
    }
  }#end of for
  
  #print(count)
  data=data.frame(cbind(t,d,z,x))
  return(data)
}


# observed data log-likelihood function for SQH

log.lik = function(b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc){ 
  
  etat = exp((b21*zt))
  etac = exp((b21*zc))
  pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
  pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
  Ft = 1-exp(-((g2*data_obs)^(1/g1)))
  Fc = 1-exp(-((g2*data_cens)^(1/g1)))
  ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  Spt=exp(-etat*pt*Ft*exp(phi))
  Spc=exp(-etac*pc*Fc*exp(phi))
  fpt=etat*pt*exp(phi)*Spt*ft
  
  log.lik.fn = sum(log(fpt))+sum(log(Spc))
  
  return(log.lik.fn)
  
}#end of lik function


# Augmented log-likelihood w.r.t. b21

l.aug.b21 = function(par=c(bb21),b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc,eps){ 
  
  etat = exp((par[1]*zt))
  etac = exp((par[1]*zc))
  pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
  pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
  Ft = 1-exp(-((g2*data_obs)^(1/g1)))
  Fc = 1-exp(-((g2*data_cens)^(1/g1)))
  ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  Spt=exp(-etat*pt*Ft*exp(phi))
  Spc=exp(-etac*pc*Fc*exp(phi))
  fpt=etat*pt*exp(phi)*Spt*ft
  
  log.lik.aug.b21 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b21)^2)
  
  return(-log.lik.aug.b21)
  
}

# Augmented log-likelihood w.r.t. b10

l.aug.b10 = function(par=c(bb10),b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc,eps){ 
  
  etat = exp((b21*zt))
  etac = exp((b21*zc))
  pt = (exp(par[1]+(b11*xt)))/(1+exp(par[1]+(b11*xt)))
  pc = (exp(par[1]+(b11*xc)))/(1+exp(par[1]+(b11*xc)))
  Ft = 1-exp(-((g2*data_obs)^(1/g1)))
  Fc = 1-exp(-((g2*data_cens)^(1/g1)))
  ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  Spt=exp(-etat*pt*Ft*exp(phi))
  Spc=exp(-etac*pc*Fc*exp(phi))
  fpt=etat*pt*exp(phi)*Spt*ft
  
  log.lik.aug.b10 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b10)^2)
  
  return(-log.lik.aug.b10)
  
}

# Augmented log-likelihood w.r.t. b11

l.aug.b11 = function(par=c(bb11),b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc,eps){ 
  
  etat = exp((b21*zt))
  etac = exp((b21*zc))
  pt = (exp(b10+(par[1]*xt)))/(1+exp(b10+(par[1]*xt)))
  pc = (exp(b10+(par[1]*xc)))/(1+exp(b10+(par[1]*xc)))
  Ft = 1-exp(-((g2*data_obs)^(1/g1)))
  Fc = 1-exp(-((g2*data_cens)^(1/g1)))
  ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  Spt=exp(-etat*pt*Ft*exp(phi))
  Spc=exp(-etac*pc*Fc*exp(phi))
  fpt=etat*pt*exp(phi)*Spt*ft
  
  log.lik.aug.b11 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b11)^2)
  
  return(-log.lik.aug.b11)
  
}


# Augmented log-likelihood w.r.t. g1

l.aug.g1 = function(par=c(gg1),b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc,eps){ 
  
  etat = exp((b21*zt))
  etac = exp((b21*zc))
  pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
  pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
  Ft = 1-exp(-((g2*data_obs)^(1/par[1])))
  Fc = 1-exp(-((g2*data_cens)^(1/par[1])))
  ft = (1/(data_obs*par[1]))*((g2*data_obs)^(1/par[1]))*(1-Ft)
  Spt=exp(-etat*pt*Ft*exp(phi))
  Spc=exp(-etac*pc*Fc*exp(phi))
  fpt=etat*pt*exp(phi)*Spt*ft
  
  log.lik.aug.g1 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g1)^2)
  
  return(-log.lik.aug.g1)
  
}


# Augmented log-likelihood w.r.t. g2

l.aug.g2 = function(par=c(gg2),b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc,eps){ 
  
  etat = exp((b21*zt))
  etac = exp((b21*zc))
  pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
  pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
  Ft = 1-exp(-((par[1]*data_obs)^(1/g1)))
  Fc = 1-exp(-((par[1]*data_cens)^(1/g1)))
  ft = (1/(data_obs*g1))*((par[1]*data_obs)^(1/g1))*(1-Ft)
  Spt=exp(-etat*pt*Ft*exp(phi))
  Spc=exp(-etac*pc*Fc*exp(phi))
  fpt=etat*pt*exp(phi)*Spt*ft
  
  log.lik.aug.g2 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g2)^2)
  
  return(-log.lik.aug.g2)
  
}


# Augmented log-likelihood w.r.t. phi

l.aug.phi = function(par=c(phi1),b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc,eps){ 
  
  etat = exp((b21*zt))
  etac = exp((b21*zc))
  pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
  pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
  Ft = 1-exp(-((g2*data_obs)^(1/g1)))
  Fc = 1-exp(-((g2*data_cens)^(1/g1)))
  ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  Spt=exp(-etat*pt*Ft*exp(par[1]))
  Spc=exp(-etac*pc*Fc*exp(par[1]))
  fpt=etat*pt*exp(par[1])*Spt*ft
  
  log.lik.aug.phi = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-phi)^2)
  
  return(-log.lik.aug.phi)
  
}


max.fn = function(p.init,b21,b10,b11,g1,g2,phi,data_obs,data_cens,zt,zc,xt,xc,eps){
  
  # Augmented log-likelihood w.r.t. b21
  
  l.aug.b21 = function(par=c(bb21)){ 
    
    etat = exp((par[1]*zt))
    etac = exp((par[1]*zc))
    pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
    pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
    Ft = 1-exp(-((g2*data_obs)^(1/g1)))
    Fc = 1-exp(-((g2*data_cens)^(1/g1)))
    ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
    Spt=exp(-etat*pt*Ft*exp(phi))
    Spc=exp(-etac*pc*Fc*exp(phi))
    fpt=etat*pt*exp(phi)*Spt*ft
    
    log.lik.aug.b21 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b21)^2)
    
    return(-log.lik.aug.b21)
    
  }
  
  # Augmented log-likelihood w.r.t. b10
  
  l.aug.b10 = function(par=c(bb10)){ 
    
    etat = exp((b21*zt))
    etac = exp((b21*zc))
    pt = (exp(par[1]+(b11*xt)))/(1+exp(par[1]+(b11*xt)))
    pc = (exp(par[1]+(b11*xc)))/(1+exp(par[1]+(b11*xc)))
    Ft = 1-exp(-((g2*data_obs)^(1/g1)))
    Fc = 1-exp(-((g2*data_cens)^(1/g1)))
    ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
    Spt=exp(-etat*pt*Ft*exp(phi))
    Spc=exp(-etac*pc*Fc*exp(phi))
    fpt=etat*pt*exp(phi)*Spt*ft
    
    log.lik.aug.b10 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b10)^2)
    
    return(-log.lik.aug.b10)
    
  }
  
  # Augmented log-likelihood w.r.t. b11
  
  l.aug.b11 = function(par=c(bb11)){ 
    
    etat = exp((b21*zt))
    etac = exp((b21*zc))
    pt = (exp(b10+(par[1]*xt)))/(1+exp(b10+(par[1]*xt)))
    pc = (exp(b10+(par[1]*xc)))/(1+exp(b10+(par[1]*xc)))
    Ft = 1-exp(-((g2*data_obs)^(1/g1)))
    Fc = 1-exp(-((g2*data_cens)^(1/g1)))
    ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
    Spt=exp(-etat*pt*Ft*exp(phi))
    Spc=exp(-etac*pc*Fc*exp(phi))
    fpt=etat*pt*exp(phi)*Spt*ft
    
    log.lik.aug.b11 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b11)^2)
    
    return(-log.lik.aug.b11)
    
  }
  
  
  # Augmented log-likelihood w.r.t. g1
  
  l.aug.g1 = function(par=c(gg1)){ 
    
    etat = exp((b21*zt))
    etac = exp((b21*zc))
    pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
    pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
    Ft = 1-exp(-((g2*data_obs)^(1/par[1])))
    Fc = 1-exp(-((g2*data_cens)^(1/par[1])))
    ft = (1/(data_obs*par[1]))*((g2*data_obs)^(1/par[1]))*(1-Ft)
    Spt=exp(-etat*pt*Ft*exp(phi))
    Spc=exp(-etac*pc*Fc*exp(phi))
    fpt=etat*pt*exp(phi)*Spt*ft
    
    log.lik.aug.g1 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g1)^2)
    
    return(-log.lik.aug.g1)
    
  }
  
  
  # Augmented log-likelihood w.r.t. g2
  
  l.aug.g2 = function(par=c(gg2)){ 
    
    etat = exp((b21*zt))
    etac = exp((b21*zc))
    pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
    pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
    Ft = 1-exp(-((par[1]*data_obs)^(1/g1)))
    Fc = 1-exp(-((par[1]*data_cens)^(1/g1)))
    ft = (1/(data_obs*g1))*((par[1]*data_obs)^(1/g1))*(1-Ft)
    Spt=exp(-etat*pt*Ft*exp(phi))
    Spc=exp(-etac*pc*Fc*exp(phi))
    fpt=etat*pt*exp(phi)*Spt*ft
    
    log.lik.aug.g2 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g2)^2)
    
    return(-log.lik.aug.g2)
    
  }
  
  
  # Augmented log-likelihood w.r.t. phi
  
  l.aug.phi = function(par=c(phi1)){ 
    
    etat = exp((b21*zt))
    etac = exp((b21*zc))
    pt = (exp(b10+(b11*xt)))/(1+exp(b10+(b11*xt)))
    pc = (exp(b10+(b11*xc)))/(1+exp(b10+(b11*xc)))
    Ft = 1-exp(-((g2*data_obs)^(1/g1)))
    Fc = 1-exp(-((g2*data_cens)^(1/g1)))
    ft = (1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
    Spt=exp(-etat*pt*Ft*exp(par[1]))
    Spc=exp(-etac*pc*Fc*exp(par[1]))
    fpt=etat*pt*exp(par[1])*Spt*ft
    
    log.lik.aug.phi = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-phi)^2)
    
    return(-log.lik.aug.phi)
    
  }
  

  b21.new = tryCatch({nlm(p=p.init[1],f=l.aug.b21)$estimate
  },error=function(e){
    b21.new = c(0)
    return(b21.new)
  }
  )
  
  
  b10.new = tryCatch({nlm(p=p.init[2],f=l.aug.b10)$estimate
  },error=function(e){
    b10.new = c(0)
    return(b10.new)
  }
  )
  
  b11.new = tryCatch({nlm(p=p.init[3],f=l.aug.b11)$estimate
  },error=function(e){
    b11.new = c(0)
    return(b11.new)
  }
  )
  
  g1.new = tryCatch({nlm(p=p.init[4],f=l.aug.g1)$estimate
  },error=function(e){
    g1.new = c(0)
    return(g1.new)
  }
  )
  
  g2.new = tryCatch({nlm(p=p.init[5],f=l.aug.g2)$estimate
  },error=function(e){
    g2.new = c(0)
    return(g2.new)
  }
  )
  
  phi.new = tryCatch({nlm(p=p.init[6],f=l.aug.phi)$estimate
  },error=function(e){
    phi.new = c(0)
    return(phi.new)
  }
  )
  
  out = c(b21.new,b10.new,b11.new,g1.new,g2.new,phi.new)
  return(out)
}


# Development of the SQH function

SQH_EWP_Wei=function(data_obs,data_cens,zt,zc,xt,xc,tol,maxit,b21,b10,b11,g1,g2,phi,epsilon1,lambda1,zeta1,rho1){
  
  p.new=rep(0,6)
  p.old=rep(0,6)
  
  p.old = c(b21,b10,b11,g1,g2,phi)
  eps = epsilon1
  
  continue = TRUE
  iter = 1
  
  while(continue){
    #print("iter")
    #print(iter)
    
    p.new = max.fn(p.init=p.old,b21=p.old[1],b10=p.old[2],b11=p.old[3],g1=p.old[4],g2=p.old[5],phi=p.old[6],
                   data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc,eps=eps)
    
    if(0%in%p.new | p.new[4]<0 | p.new[5]<0 ){
      p.new = c(0,0,0,0,0,0)
      continue = FALSE
    }else{
      tau = sum((p.new-p.old)^2)
      if((log.lik(b21=p.new[1],b10=p.new[2],b11=p.new[3],g1=p.new[4],g2=p.new[5],phi=p.new[6],data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc) -
          log.lik(b21=p.old[1],b10=p.old[2],b11=p.old[3],g1=p.old[4],g2=p.old[5],phi=p.old[6],data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc)) < (rho1*tau)){
        eps = eps*lambda1
      }else{
        eps = eps*zeta1
        p.old = p.new
        continue = (tau>tol) & (iter<maxit)
        iter = iter+1
        if(iter==maxit){
          p.new=c(0,0,0,0,0,0)
        }
      }#end of else
    }#end of else
  }#end of while
  
  result = p.new
  
  return(result)
  
}# end of the main function



######  NCG functions  #############

log.lik.fun = function(pnew=c(b21,b10,b11,g1,g2,phi),data_obs,data_cens,zt,zc,xt,xc){ 
  
  etat = exp((pnew[1]*zt))
  etac = exp((pnew[1]*zc))
  pt = (exp(pnew[2]+(pnew[3]*xt)))/(1+exp(pnew[2]+(pnew[3]*xt)))
  pc = (exp(pnew[2]+(pnew[3]*xc)))/(1+exp(pnew[2]+(pnew[3]*xc)))
  Ft = 1-exp(-((pnew[5]*data_obs)^(1/pnew[4])))
  Fc = 1-exp(-((pnew[5]*data_cens)^(1/pnew[4])))
  Spt=exp(-etat*pt*Ft*exp(pnew[6]))
  Spc=exp(-etac*pc*Fc*exp(pnew[6]))
  ft=(1/(data_obs*pnew[4]))*((pnew[5]*data_obs)^(1/pnew[4]))*(1-Ft)
  fpt=etat*pt*exp(pnew[6])*Spt*ft
  
  log.lik=sum(log(fpt))+sum(log(Spc))
  
  return(-log.lik)
  
}#end of lik function


lambda = function(pars,d.k,g.k,data_obs,data_cens,zt,zc,xt,xc,del){
  
  k = 1
  cont = TRUE
  
  while(cont){
    lam.k = 1/(2^(k-1))
    pp = pars+(d.k*lam.k)
    pp.new = c(pp[1],pp[2],pp[3],max(pp[4],0.01),max(pp[5],0.01),pp[6])
    cc1 = log.lik.fun(pp.new,data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc)
    cc2 = (log.lik.fun(pars,data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc)+(del*lam.k*sum(d.k*g.k))) 
    
    #cat("s = ", k ,",", "cc1 = ", cc1,",", "cc2 = ", cc2, "\n")
    
    if(  (cc1 != "NaN" & cc1 != Inf & cc2 != "NaN" & cc2 != "Inf") ){
      cont = (cc1 > cc2) & (k < 21) 
      k = k + 1
    }else{
      k = k + 1
      if(k>20){
        cont = FALSE
      }
    }
    
  }#end of while
  return(c(lam.k,k))
}


NCG_EWP_Wei=function(data_obs,data_cens,zt,zc,xt,xc,tol,maxit,b21,b10,b11,gam1,gam2,phi,del){
  
  p.new=rep(0,6)
  p.old=rep(0,6)
  d.old=rep(0,6)
  d.new=rep(0,6)
  
  p.old = c(b21,b10,b11,gam1,gam2,phi)
  d.old = -1*grad(log.lik.fun,p.old,method="Richardson",data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc)
  
  
  continue = TRUE
  iter=1
  
  while(continue){
    
    #print(iter)
    
    g.old = tryCatch({ grad(log.lik.fun,p.old,method="Richardson",data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc)
    },error=function(e){
      g.old = c(0,0,0,0,0,0)
      return(g.old)
    }
    )
    
    #g.old = grad(log.lik.fun,p.old,method="Richardson",data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc)
    
    lam.vec =  lambda(p.old,d.old,g.old,data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc,del=del) # lam calculated from line search algorithm
    
    
    p.new = p.old + (lam.vec[1]*d.old)
    
    dum1.new = tryCatch({ (-1*grad(log.lik.fun,p.new,method="Richardson",data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc))
    },error=function(e){
      dum1.new = c(0,0,0,0,0,0)
      return(dum1.new)
    }
    )
    
    if(lam.vec[2]>20){
      p.new = c(0,0,0,0,0,0)
      
      #cat("lam = ", lam.vec[1],",", "k = ", lam.vec[2], ",", "pnew = ", p.new, ",","norm = ",sqrt(sum(g.old*g.old)), "," , "value of fn = ", 
      #log.lik.fun(p.old,data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc), "\n")
      
      continue = FALSE
    }else{
      
      #dum1.new = (-1*grad(log.lik.fun,p.new,method="Richardson",data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc))
      g.new = -dum1.new
      
      #cat("lam = ", lam.vec[1], ",", "k = ", lam.vec[2], ",", "pnew = ", p.new, ",", "norm = ", sqrt(sum(g.new*g.new)), ",", "value of fn = ", 
      #log.lik.fun(p.new,data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc), "\n")
      
      y =  g.new - g.old
      dum2 = y - (2*d.old*sum(y*y)/sum(d.old*y))
      HG = sum(dum2*g.new)/sum(d.old*y)
      d.new = dum1.new + (HG*d.old)
      
      iter = iter + 1
      
      continue=(abs((p.new[1]-p.old[1])/p.old[1])>tol | abs((p.new[2]-p.old[2])/p.old[2])>tol | abs((p.new[3]-p.old[3])/p.old[3])>tol | 
                  abs((p.new[4]-p.old[4])/p.old[4])>tol | abs((p.new[5]-p.old[5])/p.old[5])>tol | abs((p.new[6]-p.old[6])/p.old[6])>tol) & (iter<maxit)
      
      p.old = p.new
      d.old = d.new
      #g.old = g.new
      
    }# end of else
    
  }#end of while
  
  if(iter==maxit){
    p.new=matrix(c(0,0,0,0,0,0))
  }
  
  
  return(c(p.new,iter))
  
}# end of the main function



# calling the functions

m=200 # sample size
n=500 # MC run
e1=3
pl=0.3
ph=0.9
g1_true=0.316
g2_true=0.179
phi_true= 0.7

tol=0.001
maxit=500
incr=0.20
c=0.10
del = 0.1

beta21_hat_sqh=rep(NA,n)
temp_b21_sqh=rep(NA,n)
beta10_hat_sqh=rep(NA,n)
temp_b10_sqh=rep(NA,n)
beta11_hat_sqh=rep(NA,n)
temp_b11_sqh=rep(NA,n)
g1_hat_sqh=rep(NA,n)
temp_g1_sqh=rep(NA,n)
g2_hat_sqh=rep(NA,n)
temp_g2_sqh=rep(NA,n)
phi_hat_sqh=rep(NA,n)
temp_phi_sqh=rep(NA,n)


beta21_hat_ncg=rep(NA,n)
temp_b21_ncg=rep(NA,n)
beta10_hat_ncg=rep(NA,n)
temp_b10_ncg=rep(NA,n)
beta11_hat_ncg=rep(NA,n)
temp_b11_ncg=rep(NA,n)
g1_hat_ncg=rep(NA,n)
temp_g1_ncg=rep(NA,n)
g2_hat_ncg=rep(NA,n)
temp_g2_ncg=rep(NA,n)
phi_hat_ncg=rep(NA,n)
temp_phi_ncg=rep(NA,n)


b21_true = log(e1)
b11_true = rep(NA,n)
b10_true = rep(NA,n)

count_sqh=0
count_ncg=0

for(i in 1:n){
  print(i)
  set.seed(100+i)
  
  data=sim_data_Wei(m,g1_true,g2_true,0.15,e1,pl,ph,phi_true)
  
  obs_data=data[data$d==1,]
  cens_data=data[data$d==0,]
  data_obs=obs_data$t
  data_cens=cens_data$t
  zt=obs_data$z
  zc=cens_data$z
  xt=obs_data$x
  xc=cens_data$x
  
  xmin=min(data$x)
  xmax=max(data$x)
  
  b11_true[i] = (log(ph/(1-ph))-log(pl/(1-pl)))/(xmax-xmin)
  b10_true[i] = log(pl/(1-pl))-(b11_true[i]*xmin)
  
  b21_init = sample(seq((b21_true-(incr*abs(b21_true))),(b21_true+(incr*abs(b21_true))),by=0.01),1)
  b10_init = sample(seq((b10_true[i]-(incr*abs(b10_true[i]))),(b10_true[i]+(incr*abs(b10_true[i]))),by=0.01),1)
  b11_init = sample(seq((b11_true[i]-(incr*abs(b11_true[i]))),(b11_true[i]+(incr*abs(b11_true[i]))),by=0.01),1)
  g1_init = sample(seq((g1_true-(incr*abs(g1_true))),(g1_true+(incr*abs(g1_true))),by=0.01),1)
  g2_init = sample(seq((g2_true-(incr*abs(g2_true))),(g2_true+(incr*abs(g2_true))),by=0.01),1)
  phi_init = sample(seq((phi_true-(incr*abs(phi_true))),(phi_true+(incr*abs(phi_true))),by=0.01),1)
  
  
  nr_sqh = SQH_EWP_Wei(data_obs=data_obs,data_cens=data_cens,zt=zt,zc=zc,xt=xt,xc=xc,tol=tol,maxit=1000,b21=b21_init,b10=b10_init,b11=b11_init,
                       g1=g1_init,g2=g2_init,phi=phi_init,epsilon1=100,lambda1=100,zeta1=0.50,rho1=100)
  
  nr_ncg = NCG_EWP_Wei(data_obs,data_cens,zt,zc,xt,xc,tol,maxit,b21_init,b10_init,b11_init,g1_init,g2_init,phi_init,del)
           
  
  if(0%in%nr_sqh){
    
    count_sqh=count_sqh+1
    
    beta21_hat_sqh[i]=0
    temp_b21_sqh[i]=0
    beta10_hat_sqh[i]=0
    temp_b10_sqh[i]=0
    beta11_hat_sqh[i]=0
    temp_b11_sqh[i]=0
    g1_hat_sqh[i]=0
    temp_g1_sqh[i]=0
    g2_hat_sqh[i]=0
    temp_g2_sqh[i]=0
    phi_hat_sqh[i]=0
    temp_phi_sqh[i]=0
    
  }else{
    
    beta21_hat_sqh[i]=nr_sqh[1]
    temp_b21_sqh[i]=beta21_hat_sqh[i]-b21_true
    beta10_hat_sqh[i]=nr_sqh[2]
    temp_b10_sqh[i]=beta10_hat_sqh[i]-b10_true[i]
    beta11_hat_sqh[i]=nr_sqh[3]
    temp_b11_sqh[i]=beta11_hat_sqh[i]-b11_true[i]
    g1_hat_sqh[i]=nr_sqh[4]
    temp_g1_sqh[i]=g1_hat_sqh[i]-g1_true
    g2_hat_sqh[i]=nr_sqh[5]
    temp_g2_sqh[i]=g2_hat_sqh[i]-g2_true
    phi_hat_sqh[i]=nr_sqh[6]
    temp_phi_sqh[i]=phi_hat_sqh[i]-phi_true
    
  }#end of else
  
  if(0%in%nr_ncg){
    
    count_ncg=count_ncg+1
    
    beta21_hat_ncg[i]=0
    temp_b21_ncg[i]=0
    beta10_hat_ncg[i]=0
    temp_b10_ncg[i]=0
    beta11_hat_ncg[i]=0
    temp_b11_ncg[i]=0
    g1_hat_ncg[i]=0
    temp_g1_ncg[i]=0
    g2_hat_ncg[i]=0
    temp_g2_ncg[i]=0
    phi_hat_ncg[i]=0
    temp_phi_ncg[i]=0

    
  }else{
    
    beta21_hat_ncg[i]=nr_ncg[1]
    temp_b21_ncg[i]=beta21_hat_ncg[i]-b21_true
    beta10_hat_ncg[i]=nr_ncg[2]
    temp_b10_ncg[i]=beta10_hat_ncg[i]-b10_true[i]
    beta11_hat_ncg[i]=nr_ncg[3]
    temp_b11_ncg[i]=beta11_hat_ncg[i]-b11_true[i]
    g1_hat_ncg[i]=nr_ncg[4]
    temp_g1_ncg[i]=g1_hat_ncg[i]-g1_true
    g2_hat_ncg[i]=nr_ncg[5]
    temp_g2_ncg[i]=g2_hat_ncg[i]-g2_true
    phi_hat_ncg[i]=nr_ncg[6]
    temp_phi_ncg[i]=phi_hat_ncg[i]-phi_true
    
  }#end of else
  
  
}#end of for


  
avg_b21_sqh=sum(beta21_hat_sqh)/(n-count_sqh)
avg_b10_sqh=sum(beta10_hat_sqh)/(n-count_sqh)
avg_b11_sqh=sum(beta11_hat_sqh)/(n-count_sqh)
avg_g1_sqh=sum(g1_hat_sqh)/(n-count_sqh)
avg_g2_sqh=sum(g2_hat_sqh)/(n-count_sqh)
avg_phi_sqh = sum(phi_hat_sqh)/(n-count_sqh)


bias_b21_sqh=sum(temp_b21_sqh)/(n-count_sqh)
bias_b10_sqh=sum(temp_b10_sqh)/(n-count_sqh)
bias_b11_sqh=sum(temp_b11_sqh)/(n-count_sqh)
bias_g1_sqh=sum(temp_g1_sqh)/(n-count_sqh)
bias_g2_sqh=sum(temp_g2_sqh)/(n-count_sqh)
bias_phi_sqh=sum(temp_phi_sqh)/(n-count_sqh)


rmse_b21_sqh=sqrt(sum(temp_b21_sqh^2)/(n-count_sqh-1))
rmse_b10_sqh=sqrt(sum(temp_b10_sqh^2)/(n-count_sqh-1))
rmse_b11_sqh=sqrt(sum(temp_b11_sqh^2)/(n-count_sqh-1))
rmse_g1_sqh=sqrt(sum(temp_g1_sqh^2)/(n-count_sqh-1))
rmse_g2_sqh=sqrt(sum(temp_g2_sqh^2)/(n-count_sqh-1))
rmse_phi_sqh=sqrt(sum(temp_phi_sqh^2)/(n-count_sqh-1))



avg_b21_ncg=sum(beta21_hat_ncg)/(n-count_ncg)
avg_b10_ncg=sum(beta10_hat_ncg)/(n-count_ncg)
avg_b11_ncg=sum(beta11_hat_ncg)/(n-count_ncg)
avg_g1_ncg=sum(g1_hat_ncg)/(n-count_ncg)
avg_g2_ncg=sum(g2_hat_ncg)/(n-count_ncg)
avg_phi_ncg = sum(phi_hat_ncg)/(n-count_ncg)


bias_b21_ncg=sum(temp_b21_ncg)/(n-count_ncg)
bias_b10_ncg=sum(temp_b10_ncg)/(n-count_ncg)
bias_b11_ncg=sum(temp_b11_ncg)/(n-count_ncg)
bias_g1_ncg=sum(temp_g1_ncg)/(n-count_ncg)
bias_g2_ncg=sum(temp_g2_ncg)/(n-count_ncg)
bias_phi_ncg=sum(temp_phi_ncg)/(n-count_ncg)


rmse_b21_ncg=sqrt(sum(temp_b21_ncg^2)/(n-count_ncg-1))
rmse_b10_ncg=sqrt(sum(temp_b10_ncg^2)/(n-count_ncg-1))
rmse_b11_ncg=sqrt(sum(temp_b11_ncg^2)/(n-count_ncg-1))
rmse_g1_ncg=sqrt(sum(temp_g1_ncg^2)/(n-count_ncg-1))
rmse_g2_ncg=sqrt(sum(temp_g2_ncg^2)/(n-count_ncg-1))
rmse_phi_ncg=sqrt(sum(temp_phi_ncg^2)/(n-count_ncg-1))


bias_sqh = c(bias_b21_sqh,bias_b10_sqh,bias_b11_sqh,bias_g1_sqh,bias_g2_sqh,bias_phi_sqh)
bias_ncg = c(bias_b21_ncg,bias_b10_ncg,bias_b11_ncg,bias_g1_ncg,bias_g2_ncg,bias_phi_ncg)
rmse_sqh = c(rmse_b21_sqh,rmse_b10_sqh,rmse_b11_sqh,rmse_g1_sqh,rmse_g2_sqh,rmse_phi_sqh)
rmse_ncg = c(rmse_b21_ncg,rmse_b10_ncg,rmse_b11_ncg,rmse_g1_ncg,rmse_g2_ncg,rmse_phi_ncg)

round(data.frame(bias_sqh,bias_ncg,rmse_sqh,rmse_ncg),3)


print("count_sqh is:")
print(count_sqh)

print("count_ncg is:")
print(count_ncg)




