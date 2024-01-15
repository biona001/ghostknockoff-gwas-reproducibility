#Enable command line arguments
args = commandArgs(TRUE)
out_type = c('C','D')[as.numeric(args[1])]

print(out_type)

#ml R/4.0.2

#Load SKAT haplotypes
library(SKAT)
library(Matrix)
library(numDeriv)
library(glmnet)
library(corpcor)
library(susieR)
library(CVXR)
library(ghostbasil)
library(MASS)

source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')
#load('/oak/stanford/groups/zihuai/ESGWAS_lasso/sysdata.rda')

library(seqminer)
library(data.table)

#load saved data
data.all<-read.table(paste0('/oak/stanford/groups/zihuai/ESGWAS_lasso/Simulation/data_GhostKnockoff_EUR.txt'),header=T)

result.all<-c();result.detail<-c()
a<-1;i<-1;N.causal<-10
for (i in 1:1000){
  if(i %in%seq(1,100,by=10)){
    #sample a subset of individuals and variants
    n<-6000
    
    sample.index<-sample(1:nrow(data.all),n)
    data.total<-data.all[sample.index,,drop=F]
    data<-data.total[1:3000,]
    data.ref<-data.total[-(1:3000),]
    cor.G<-matrix(as.numeric(corpcor::cor.shrink(data.ref,verbose=F)), nrow=ncol(data.ref))
    
    #basic operation
    MAF<-apply(data,2,mean)/2
    data[,MAF>0.5]<-2-data[,MAF>0.5]
    pos<-as.numeric(sub('.*[.]','',colnames(data)))
    G<-as.matrix(data)
    MAF<-apply(G,2,mean)/2
    MAC<-apply(G,2,sum)
    sd.G<-apply(G,2,sd)
    N.SNP<-ncol(G)
    G<-Matrix(G,sparse=T)

    #multiple knockoffs
    fit.prelim.M.sdp<-GhostKnockoff.prelim(cor.G,M=1,method='sdp')

    temp.G<-scale(G)
    
    #multiple second order knockoffs - sdp
    mu.G<-colMeans(temp.G)
    sd.G<-apply(temp.G,2,sd)
    scale.G<-t((t(temp.G)-mu.G)/sd.G)
    M<-1
    E<-matrix(rnorm(nrow(scale.G)*M*ncol(scale.G)),nrow(scale.G),M*ncol(scale.G))
    scale.G_k1.M<-t(apply(scale.G%*%t(fit.prelim.M.sdp$P.each),1,rep,times=M))+E%*%t(fit.prelim.M.sdp$V.left)
    temp.G_k1.M<-t(rep(sd.G,M)*t(scale.G_k1.M)+rep(mu.G,M))
    
    #generate causal
    causal.index<-(1:N.SNP)%in%sample(1:N.SNP,N.causal)
    beta<-as.matrix(sqrt(a/sum(causal.index)/apply(G,2,var))*causal.index)
  }
  
  #Generate phenotype
  if(out_type=='C'){
    X1<-NULL
    mu<-rnorm(nrow(G),0,3)
    corariates<-X1
    y<-as.matrix(mu+G%*%beta)
  }else{
    Prevalence<-0.1;#;OR.f1<-9;OR.f2<-1#2
    b0_pheno<-log(Prevalence/(1-Prevalence))#for binary trait
    X1<-NULL
    mu<-G%*%beta+rnorm(nrow(G),0,1)
    mu<-mu-mean(mu)
    y<-as.matrix(rbinom(nrow(G),1, prob=as.numeric(exp(mu+b0_pheno)/(1+exp(mu+b0_pheno)))))
  }
  temp.X1<-X1
  temp.y<-scale(y)

  result.prelim<-KS.prelim(y, X=X1, out_type=out_type)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  
  #second order multiple knockoffs - sdp
  Zscore_0<-t(temp.G)%*%temp.y/sqrt(n)
  Zscore_k<-matrix(t(temp.G_k1.M)%*%temp.y/sqrt(n),length(Zscore_0),M)
  T_0<-Zscore_0^2
  T_k<-Zscore_k^2
  MK.stat<-MK.statistic(T_0,T_k)
  q1<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)
  
  #second order multiple knockoffs - lasso
  cv.fit<-cv.glmnet(scale(cbind(as.matrix(temp.G),as.matrix(temp.G_k1.M),temp.X1)),scale(temp.y),standardize = T,alpha=1)
  lasso.fit<-glmnet(scale(cbind(as.matrix(temp.G),as.matrix(temp.G_k1.M),temp.X1)),scale(temp.y),standardize = T,alpha=1,lambda=cv.fit$lambda.min)
  T_0<-abs(lasso.fit$beta[1:ncol(temp.G)])
  T_k<-abs(matrix(lasso.fit$beta[-c(1:ncol(temp.G))],ncol(temp.G),M))
  MK.stat<-MK.statistic(T_0,T_k)
  q2<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)#rep(1,length(q1))#
  
  #new method - SummaryStat multiple knockoffs
  n<-length(Y.res)
  Zscore_0<-t(temp.G)%*%temp.y/sqrt(n)
  #new method - SummaryStat sdp knockoff + marginal
  ES.stat<-GhostKnockoff.fit(Zscore_0=Zscore_0,N.effect=n,fit.prelim=fit.prelim.M.sdp,method='marginal')
  T_0<-ES.stat$T_0[[1]]
  T_k<-ES.stat$T_k[[1]]
  MK.stat<-MK.statistic(T_0,T_k)
  q3<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)
  
  #new method - SummaryStat sdp knockoff + marginal - zhaomeng's implementation
  marginal_ghost <- function(Z, y_norm, n, P, V){
    return(c(abs(Z), abs(as.vector(t(P)%*%Z) + y_norm*mvrnorm(1,rep(0,nrow(V)),V))))
  }
  Z_zhaomeng<-t(temp.G)%*%temp.y
  y_norm<-sqrt(sum(temp.y^2))
  n<-length(Y.res)
  P<-t(fit.prelim.M.sdp$P.each)
  V<-fit.prelim.M.sdp$V.left%*%t(fit.prelim.M.sdp$V.left)
  T_all<-matrix(marginal_ghost(Z=Z_zhaomeng,y_norm=y_norm,n=n,P=P,V=V),ncol=M+1)^2
  T_0<-T_all[,1,drop=F]
  T_k<-T_all[,-1,drop=F]
  MK.stat<-MK.statistic(T_0,T_k)
  q4<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)
  
  #new method - SummaryStat sdp knockoff + lasso
  ES.stat<-GhostKnockoff.fit(Zscore_0=Zscore_0,N.effect=n,fit.prelim=fit.prelim.M.sdp,method='lasso')
  T_0<-ES.stat$T_0[[1]]
  T_k<-ES.stat$T_k[[1]]
  MK.stat<-MK.statistic(T_0,T_k)
  q5<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)
  
  #new method - SummaryStat sdp knockoff + lasso - zhaomeng's implementation
  get_sigma <- function(Z, y_norm, Sigma_0_inv, n, p){
    Z_full = Z
    return (sqrt(max((p+n+1)*(y_norm^2)/(n^2+n) - t(Z_full)%*%Sigma_0_inv%*%Z_full/(n^2+n),0)))
  }
  get_lambda_1 <- function(sigma, chol_Sigma_0, M = 1, n, p, kappa = 0.3){
    X = matrix(rnorm(n*(M+1)*p),n,(M+1)*p)%*%chol_Sigma_0
    n = nrow(X)
    eps = rnorm(n)
    return (kappa*max(abs(t(X)%*%eps))*sigma/n)
  }
  get_lambda <- function(sigma, chol_Sigma_0, M = 1, n, p, kappa = 0.3, num=10){
    mean(replicate(num, get_lambda_1(sigma, chol_Sigma_0, M, n, p, kappa)))
  }
  lasso_ghost_new <- function(Z, y_norm, cor.G, Sigma_0, Sigma_0_inv,chol_Sigma_0, P, V, M = 1, n, p, kappa = 0.3){
    Z_full = rbind(Z, t(P)%*%Z + y_norm*mvrnorm(1,rep(0,p),V))
    beta = Variable((M+1)*p)
    sigma = get_sigma(Z, y_norm, Sigma_0_inv=solve(cor.G), n, p)
    lambda0 = get_lambda(sigma, chol_Sigma_0, M, n, p, kappa)
    obj = 1/2*quad_form(beta, Sigma_0) - t(Z_full)%*%beta/n + lambda0*cvxr_norm(beta, 1)
    prob = Problem(Minimize(obj))
    result = solve(prob)
    beta_ghost = result$getValue(beta)
    return(abs(beta_ghost))
  }
  Z_zhaomeng<-t(temp.G)%*%temp.y
  y_norm<-sqrt(sum(temp.y^2))
  n<-length(Y.res)
  p<-length(Z_zhaomeng)
  P<-t(fit.prelim.M.sdp$P)
  V<-fit.prelim.M.sdp$V.left%*%t(fit.prelim.M.sdp$V.left)
  Sigma_0<-fit.prelim.M.sdp$A+diag(0.01,nrow(fit.prelim.M.sdp$A))
  Sigma_0_inv = solve(Sigma_0)
  chol_Sigma_0 = chol(Sigma_0)
  ans = lasso_ghost_new(Z_zhaomeng, y_norm, cor.G, Sigma_0, Sigma_0_inv,chol_Sigma_0, P, V, M = 1, n, p, kappa = 0.6)
  
  T_all<-matrix(ans,ncol=M+1)
  T_0<-T_all[,1,drop=F]
  T_k<-T_all[,-1,drop=F]
  MK.stat<-MK.statistic(T_0,T_k)
  q6<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)
  
  #new method - SummaryStat sdp knockoff + lasso - zhaomeng's tunning parameter + ghostbasil
  Sigma_0<-as.matrix(fit.prelim.M.sdp$A)
  shrink<-0.01
  Sigma_0<-Sigma_0+diag(shrink,nrow(Sigma_0)) #shrinkage
  Sigma_0_inv = solve(Sigma_0)
  chol_Sigma_0 = chol(Sigma_0)
  sigma = get_sigma(Z_zhaomeng, y_norm, Sigma_0_inv=solve(cor.G), n, p)
  lambda = get_lambda(sigma, chol_Sigma_0, M = 1, n, p, kappa = 0.6, num=10)
  #generate knockoff z-scores
  Normal_k<-matrix(fit.prelim.M.sdp$Normal_50Studies[,1],nrow=p)
  GK.Zscore_0<-Zscore_0
  GK.Zscore_k<-as.vector(fit.prelim.M.sdp$P%*%GK.Zscore_0)+Normal_k
  r<-GK.Zscore_0/sqrt(n)#sqrt(N.effect-1+GK.Zscore_0^2)
  r_k<-as.vector(GK.Zscore_k/sqrt(n))#sqrt(N.effect-1+GK.Zscore_k^2))
  r_all<-as.matrix(c(r,r_k))
  #Get lambda sequence
  fit.basil<-ghostbasil(Sigma_0, r_all, alpha=1)
  lambda.all<-fit.basil$lmdas
  lambda.seq <- lambda.all[lambda.all > lambda]
  lambda.seq <- c(lambda.seq, lambda)
  #fit basil
  fit.basil<-ghostbasil(Sigma_0, r_all,user.lambdas=lambda.seq, alpha=1)
  beta.fit<-fit.basil$betas[,ncol(fit.basil$betas)]
  #knockoff filter
  T_all<-matrix(abs(beta.fit),ncol=M+1)
  T_0<-T_all[,1,drop=F]
  T_k<-T_all[,-1,drop=F]
  MK.stat<-MK.statistic(T_0,T_k)
  q7<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)
  
  #new method - SummaryStat sdp knockoff + susie
  ES.stat<-GhostKnockoff.fit(Zscore_0=Zscore_0,N.effect=n,fit.prelim=fit.prelim.M.sdp,method='susie')
  T_0<-ES.stat$T_0[[1]]
  T_k<-ES.stat$T_k[[1]]
  MK.stat<-MK.statistic(T_0,T_k)
  q8<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],ncol(T_k),Rej.Bound=10000)
  
  Q<-cbind(q1,q2,q3,q4,q5,q6,q7,q8);Q[is.na(Q)]<-0

  result<-c()
  for(fdr in seq(0,0.2,by=0.01)){
    Get.select<-function(x){Q[,x] <= fdr}
    selected<-sapply(1:ncol(Q),Get.select);selected[is.na(selected)]<-F
    signal<-(beta!=0)
    
    power<-function(x){sum(signal&x)/max(1,sum(signal))}
    fdp<-function(x){sum((!signal)&x)/max(1,sum(x))}
    
    result<-rbind(result,c(apply(selected,2,power),apply(selected,2,fdp)))
  }
  if(i==1){result.fdr<-result}else{result.fdr<-result.fdr+result}
  print(i)
  temp.print<-cbind(N.causal,seq(0,0.2,by=0.01),result.fdr/i)
  print(temp.print)
  
  method.name<-c('SecondOrderM_sdp_marginal','SecondOrderM_sdp_lasso',
                 'SummaryStatM_sdp_marginal','SummaryStatM_sdp_marginal-zhaomeng',
                 'SummaryStatM_sdp_lasso','SummaryStatM_sdp_lasso-zhaomeng','SummaryStatM_sdp_lasso-zhaomeng tuning + basil',
                 'SummaryStatM_sdp_susie')
  temp.result<-cbind(i,a,N.causal,seq(0,0.2,by=0.01),result.fdr/i)
  colnames(temp.result)<-c('replicate','a','N.causal','fdr',
                           paste0(method.name,'_power'),
                           paste0(method.name,'_fdr'))
  temp.result<-round(temp.result,digits=3)
  write.table(temp.result,paste0('/oak/stanford/groups/zihuai/ESGWAS_lasso/Simulation/result_GhostKnockoff_',out_type,'.txt'),row.names=F,col.names=T,quote=F,sep='\t')
  
  result.detail<-rbind(result.detail,result[6,])
  colnames(result.detail)<-c(paste0(method.name,'_power'),
                             paste0(method.name,'_fdr'))
  write.table(result.detail,paste0('/oak/stanford/groups/zihuai/ESGWAS_lasso/Simulation/result_detail_GhostKnockoff_FDR05_',out_type,'.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}




