library(bayesPO)
library(bayesplot)
library(ggplot2)
library(coda)





geraPO = function(b,d,ls = 1e4,covcom=0,link="logit"){

  nb = length(b)
  nd = length(d)
  ## Gerando os pontos totais
  Nt = rpois(1,ls)
  # Espalhando os pontos na regi?o
  xT = runif(Nt)
  yT = runif(Nt)
  pT = cbind(xT,yT)

  ## Gerando valores da covari?vel Z
  ZT = matrix(c(rep(1,Nt),rnorm(Nt*(nb-1))),nrow=Nt,ncol=nb,byrow=F)

  ## Reamostra de acordo com a fun??o liga??o
  if (link=="logit") probY = -log(1+exp(-ZT%*%b)) else probY = ZT%*%b-log(ls)
  selY = 1*(log(runif(Nt))<=probY)
  Y = pT[which(selY==1),]
  ZY = as.matrix(ZT[which(selY==1),])
  NY = nrow(Y)
  U = pT[which(selY==0),]
  ZU = as.matrix(ZT[which(selY==0),])
  NU = nrow(U)

  ## Gerando valores da covari?vel W
  if (covcom) WY = matrix(c(rep(1,NY),ZY[,2],rnorm(NY*(nd-2))),nrow=NY,ncol=nd,byrow=F)
  else WY = matrix(c(rep(1,NY),rnorm(NY*(nd-1))),nrow=NY,ncol=nd,byrow=F)

  ## Reamostra para separar quem foi observado e quem n?o (amostra PO)
  probX = -log(1+exp(-WY%*%d))
  selX = 1*(log(runif(NY))<=probX)
  X = Y[which(selX==1),]
  ZX = as.matrix(ZY[which(selX==1),])
  WX = as.matrix(WY[which(selX==1),])
  NX = nrow(X)
  Xl = Y[which(selX==0),]
  ZXl = as.matrix(ZY[which(selX==0),])
  WXl = as.matrix(WY[which(selX==0),])
  NXl = nrow(Xl)

    #	return(caminho)
  return(list(zx=ZX[,-1],wx=WX[,-1],zXp=ZXl[,-1],wXp=WXl[,-1],zU=ZU[,-1]))
}

set.seed(123)
background = matrix(rnorm(4e6),ncol=2)
betas = 1:2; deltas = 1:2; truels = 1e3
obs = geraPO(betas,deltas,truels)
obs_mat = cbind(obs$zx,obs$wx)
areaD = 1

model = bayesPO_model(po = obs_mat,intensitySelection = 1,
                      observabilitySelection = 2,
                      initial_values = 4)

output = fit_bayesPO(model,background,mcmc_setup = list(burnin = 5e3,n_iter = 5e3))

output

post = as.array(output)
mcmc_trace(post)
mcmc_dens(post)
mcmc_pairs(post)

post_df = as.data.frame(output)
ggplot(post_df,aes(beta_0,beta_1))+theme_bw()+geom_point(aes(color=chain),size=0.3)

