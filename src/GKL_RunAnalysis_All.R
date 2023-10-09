args = commandArgs(TRUE)
replicate = as.numeric(args[1])
zscore.file = args[2]
SummaryStat.info.file = args[3]
cS2G.file = args[4]
out.dir = args[5]

MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  
  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1
  
  if(method=='max'){tau<-apply(T.temp,1,max)-apply(T.temp,1,max.nth,n=2)}
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}

MK.q.byStat<-function (kappa,tau,M,clusters=1:length(kappa),Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    G.factor<-max(table(clusters[b][1:i]))
    temp_ratio<-(1/M*G.factor+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));
  if(length(which(tau[b]>0))!=0){
    index_bound<-max(which(tau[b]>0))
    for(i in 1:length(b)){
      temp.index<-i:min(length(b),Rej.Bound,index_bound)
      if(length(temp.index)==0){next}
      q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
      if(i>Rej.Bound){break}
    }
    q[q>1]<-1
  }
  
  return(q)
}

#ancestry = c('EUR','AFR','EAS')[as.numeric(args[2])]

print('start; July 2023')

#ml R/4.0.2 cmake/3.24.2 harfbuzz/1.4.8 fribidi/1.0.12 libgit2/1.1.0 openssl/3.0.7
#install.packages("devtools", repos='http://cran.us.r-project.org', Ncpus=4)
#load_all('/oak/stanford/groups/zihuai/ESGWAS_lasso/ghostbasil/R/')
#install.packages("libfribidi-dev", repos='http://cran.us.r-project.org', Ncpus=4)

#libcurl4-openssl-dev

library(data.table)
library(ghostbasil)
library(Matrix)
library(susieR)
library(rhdf5)
# source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
# source('/oak/stanford/groups/zihuai/ESGWAS/ESGWAS.R')
# source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')

#load z-scores
M<-5 # number of knockoffs
# Zscores<-as.data.frame(fread('/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/AD_Zscores_Meta.txt'))
Zscores<-as.data.frame(fread(zscore.file))
Zscores<-Zscores[,-ncol(Zscores)]
colnames(Zscores)[8+replicate]<-'Z'
Zscores<-Zscores[!is.na(Zscores[,'Z']),]
#exclude APOE
#LD_block<-read.table('/oak/stanford/groups/zihuai/GeneticsResources/LD_block/EUR/fourier_ls-all.bed',header=T)
#LD_block[which(LD_block[,1]=='chr19' & LD_block[,2]<=45409053 & LD_block[,3]>=45412650),]
#exclude.start<-as.numeric(LD_block[which(LD_block[,1]=='chr19' & LD_block[,2]<=45409053 & LD_block[,3]>=45412650),2])
#exclude.end<-LD_block[which(LD_block[,1]=='chr19' & LD_block[,2]<=45409053 & LD_block[,3]>=45412650),3]
#exclude.chr<-LD_block[which(LD_block[,1]=='chr19' & LD_block[,2]<=45409053 & LD_block[,3]>=45412650),1]
#Zscores<-Zscores[-which(Zscores[,3]=='19' & Zscores[,4]>=exclude.start & Zscores[,4]<=exclude.end),]

# SummaryStat.info<-read.table('/oak/stanford/groups/zihuai/ESGWAS_lasso/AD_Analysis/SummaryStatInfo.txt',header=T,sep='\t')
SummaryStat.info<-read.table(SummaryStat.info.file,header=T,sep='\t')
if(replicate<=10){N.effect<-round(SummaryStat.info[replicate,'SampleSize'])}else{
  N.effect<-506200 # effective sample size
}
print(N.effect)

#lasso path
lambda_max <- max(abs(Zscores[,'Z']/sqrt(N.effect)))
epsilon <- .0001
K <- 100
lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                            length.out = K)), digits = 10)

# match reference panel variants
chr<-Zscores[,'hg19.chr']
pos<-Zscores[,'hg19.pos']
p.chr<-c()
for(j in 1:22){
  print(j)
  bim.filename<-paste0('/oak/stanford/groups/zihuai/UKB_data/genotyped_call/', "ukb_cal_chr", j, "_v2.bim")
  bim.data<-read.table(bim.filename)
  tested.pos<-as.matrix(bim.data[,4])
  p.chr<-c(p.chr,sum(pos[chr==j]%in%tested.pos))
}
p.genome<-sum(p.chr)
#p.genome<-784256 # currently this is the total number of UKB typed variants

#expectation of extreme values
lambda.sample<-c()
for(j in 1:10){
  lambda.sample<-c(lambda.sample,0.6*max(abs(rnorm(p.genome*(M+1))))/sqrt(N.effect))
}
lambda<-mean(lambda.sample)
print(lambda)

lambda.seq <- lambdapath[lambdapath > lambda]
lambda.seq <- c(lambda.seq, lambda)

#cS2G file
cS2G.QC<-fread(cS2G.file)
# cS2G.QC<-fread('/oak/stanford/groups/zihuai/XinranQi/cS2G_UKBB/topcS2G_allVariants/topcS2GGene_allVariants.csv')

#options(warn=2)

hg<-'hg19'
options(scipen=999)
group.results.all<-c()
results.all<-c();total.time<-0
mat.list<-list();id<-c()
z<-c();z_A<-c();z_B<-c();
susie.error<-0
#for zhaomeng validation
sigma_stats<-c()
extreme_stats<-c()
heritability<-c()

for(initial.chr in 21:22){
  print(paste0('chr',initial.chr))
  
  #KF.dir<-paste0('/oak/stanford/groups/zihuai/pan_ukb_group_knockoffs/maxent_hc/chr',initial.chr,'/')
  KF.dir<-paste0('/oak/stanford/groups/zihuai/pan_ukb_group_knockoffs/EUR/chr',initial.chr,'/')
  #KF.dir<-paste0('/oak/stanford/groups/zihuai/pan_ukb_group_knockoffs/EUR_hc_completeLinkage/chr',initial.chr,'/')
  KF.SLD.names<-list.files(KF.dir,pattern = 'LD')#paste0('SLD',sub('.txt','',sub('Info','',KF.info.names)),'.gds')
  KF.info.names<-paste0('Info',sub('.h5','',sub('LD','',KF.SLD.names)),'.csv')#list.files(KF.dir,pattern = 'Info')
  N_set<-length(KF.info.names)
  #set<-64
  
  for(set in 1:N_set){
    t1<-proc.time()
    print(set)
    
    results<-c()
    #load h5 data
    #KF.data = H5Fopen(paste0(KF.dir,KF.SLD.names[set]))
    D<-h5read(paste0(KF.dir,KF.SLD.names[set]),'D')#KF.data$D
    S<-h5read(paste0(KF.dir,KF.SLD.names[set]),'S')#KF.data$S
    Sigma<-h5read(paste0(KF.dir,KF.SLD.names[set]),'Sigma')#KF.data$Sigma
    group_reps<-h5read(paste0(KF.dir,KF.SLD.names[set]),'group_reps')#KF.data$group_reps
    SigmaInv<-h5read(paste0(KF.dir,KF.SLD.names[set]),'SigmaInv')#KF.data$SigmaInv
    groups<-h5read(paste0(KF.dir,KF.SLD.names[set]),'groups')#KF.data$groups
    #H5Fclose(KF.data)
    
    KF.info<-as.data.frame(fread(paste0(KF.dir,KF.info.names[set])))
    if(sum(!is.na(KF.info[,paste0('pos_',hg)]))==0){next}
    chr<-KF.info[!is.na(KF.info[,'chr']),'chr'][1]
    pos<-KF.info[,paste0('pos_',hg)]
    ##preparations - info
    clusters<-groups
    ifRep<-rep(0,nrow(KF.info))
    ifRep[group_reps]<-1
    KF.info<-cbind(KF.info,clusters,ifRep)
    start<-min(KF.info[,'pos_hg19']);end<-max(KF.info[,'pos_hg19'])
    #add functional scores
    # temp_dir<-'/oak/stanford/groups/zihuai/ESGWAS/AD_results/TempScores/'
    # regBase.filename<-'/oak/stanford/groups/zihuai/fredlu/regBase/V1.1/regBase_Common.V1.1.gz'
    # job_id<-paste0(chr,':',start,'-',end)
    # dir.string<-paste0('tabix ',regBase.filename,' ',sub(".*chr", "", chr),':',start,'-',end,' > ',temp_dir,'Temp_FA_AD_',replicate,'.txt')
    # system(dir.string)
    # score<-try(data.frame(fread(paste0(temp_dir,'Temp_FA_AD_',replicate,'.txt'))),silent = T)
    # info<-score[,c(1,3,4,5)]
    # score<-as.matrix(score[,c(17,19,21,23,25,29,31,35)])
    # colnames(score)<-c('CADD','DANN','FATHMM-MKL','FunSeq2','Eigen','GenoCanyon','FIRE','LINSIGHT')
    # index<-match(paste0(KF.info[,'chr'],':',KF.info[,'pos_hg19'],'-',KF.info[,'ref'],'-',KF.info[,'alt']),
    #              paste0(info[,1],':',info[,2],'-',info[,3],'-',info[,4]))
    # score<-score[index,]
    # KF.info<-cbind(KF.info,score)
    #add cS2G
    temp.cS2G.QC<-data.frame(cS2G.QC[which(cS2G.QC[,'hg19.chr']==chr & cS2G.QC[,'hg19.pos']<= end & cS2G.QC[,'hg19.pos']>= start),])
    KF.info<-cbind(KF.info,temp.cS2G.QC[match(KF.info[,'pos_hg19'],as.numeric(temp.cS2G.QC[,'hg19.pos'])),5:7])

    #load Zscores
    Zscores.chr<-Zscores[Zscores[,paste0(hg,'.chr')]==chr,]
    temp.Zscores<-Zscores.chr[which(Zscores.chr[,paste0(hg,'.pos')] >= min(pos,na.rm=T) & Zscores.chr[,paste0(hg,'.pos')] <= max(pos,na.rm=T)),]
    
    #match Zscores
    match.index<-match(pos,temp.Zscores[,paste0(hg,'.pos')])
    KF.Z<-temp.Zscores[match.index,'Z'] #meta-analysis Z-scores
    #match ref and alt allele
    reverse.index<-which(temp.Zscores[match.index,'ref']==KF.info[,'alt'] & temp.Zscores[match.index,'alt']==KF.info[,'ref'])
    KF.Z[reverse.index]<--KF.Z[reverse.index]
    NA.index<-which(!((temp.Zscores[match.index,'ref']==KF.info[,'ref'] & temp.Zscores[match.index,'alt']==KF.info[,'alt']) | (temp.Zscores[match.index,'ref']==KF.info[,'alt'] & temp.Zscores[match.index,'alt']==KF.info[,'ref'])))
    KF.Z[NA.index]<-NA
    KF.Z[is.na(KF.info[,'chr'])]<-NA
    ##remove large inconsistencies (optional; being tested; heuristic procedure)
    #option 1 (susie paper)
    #omega<-solve(Sigma[!is.na(KF.Z),!is.na(KF.Z)])
    #diag.omega<-diag(omega)
    #diag(omega)<-0
    #t.stat<-sqrt(diag.omega)*(KF.Z[!is.na(KF.Z)]+omega%*%KF.Z[!is.na(KF.Z)]/diag.omega)
    #KF.Z[which(!is.na(KF.Z))[abs(t.stat)>3 & abs(KF.Z[!is.na(KF.Z)])>2]]<-NA
    #option 2 (DENTIST paper)
    #t.stat<-(KF.Z-Sigma[,which.max(KF.Z)]*KF.Z[which.max(KF.Z)])^2/(1-Sigma[,which.max(KF.Z)]^2)
    #p.stat<-pchisq(t.stat,lower.tail=F)
    #KF.Z[which(p.stat<1e-4 & abs(Sigma[,which.max(KF.Z)])>0.6)]<-NA
    
    #update info file
    KF.info<-cbind(KF.info,KF.Z)
    
    #remove variants that do not show up
    index<-which(!is.na(KF.info[,'KF.Z']))
    rep.index<-match(which(KF.info[,'ifRep']==1 & !is.na(KF.info[,'KF.Z'])),which(KF.info[,'ifRep']==1))
    if(length(index)<=1 | length(rep.index)<=1){next}
    KF.info<-KF.info[index,]
    
    KF.info[,'clusters']<-match(KF.info[,'clusters'],unique(KF.info[,'clusters']))
    KF.Z<-KF.info[,'KF.Z']
    KF.Rep.Z<-KF.info[which(KF.info[,'ifRep']==1),'KF.Z']

    #all data matrices
    KF.All.data<-cbind(D,Sigma)
    KF.All.data<-KF.All.data[index,c(index,index+nrow(KF.All.data))]
    #representative variants data matrices
    KF.Rep.data<-cbind(S,Sigma[group_reps,group_reps])
    KF.Rep.data<-KF.Rep.data[rep.index,c(rep.index,rep.index+nrow(KF.Rep.data))]
    
    #shrinkage for consistency
    temp.z<-as.matrix(KF.Rep.Z)
    eigen.fit<-eigen(KF.Rep.data[,(nrow(KF.Rep.data)+1):(2*nrow(KF.Rep.data))])
    #temp.z<-as.matrix(KF.Z)
    #eigen.fit<-eigen(KF.All.data[,(nrow(KF.All.data)+1):(2*nrow(KF.All.data))])
    Eloglik<-function(lambda,eigen.fit,z){
      znull<-t(eigen.fit$vectors)%*%z
      return(-sum(log((1-lambda)*eigen.fit$values+lambda))-t(znull)%*%(((1-lambda)*eigen.fit$values+lambda)^(-1)*znull))
    }
    shrinkage.lambda<-optimize(Eloglik,c(0,1),eigen.fit=eigen.fit,z=temp.z,maximum=T)$maximum
    print(shrinkage.lambda)
    #shrinkage.lambda<-0.01
    
    #update matrices
    #KF.All.data[,1:nrow(KF.All.data)]<-(1-shrinkage.lambda)*KF.All.data[,1:nrow(KF.All.data)]+shrinkage.lambda*(M+1)/M*diag(1,nrow(KF.All.data))
    #KF.All.data[,(nrow(KF.All.data)+1):(2*nrow(KF.All.data))]<-(1-shrinkage.lambda)*KF.All.data[,(nrow(KF.All.data)+1):(2*nrow(KF.All.data))]+shrinkage.lambda*diag(1,nrow(KF.All.data))
    
    #KF.Rep.data[,1:nrow(KF.Rep.data)]<-(1-shrinkage.lambda)*KF.Rep.data[,1:nrow(KF.Rep.data)]+shrinkage.lambda*(M+1)/M*diag(1,nrow(KF.Rep.data))
    #KF.Rep.data[,(nrow(KF.Rep.data)+1):(2*nrow(KF.Rep.data))]<-(1-shrinkage.lambda)*KF.Rep.data[,(nrow(KF.Rep.data)+1):(2*nrow(KF.Rep.data))]+shrinkage.lambda*diag(1,nrow(KF.Rep.data))
    
    #prepare initial inverse for efficient calculation
    inv.Rep<-as.matrix(solve(KF.Rep.data[,(nrow(KF.Rep.data)+1):(2*nrow(KF.Rep.data))]))
    
    ## Generate knockoff Z-scores for the representative variants
    
    #correlation matrix for representative variables
    cor.G<-KF.Rep.data[,(nrow(KF.Rep.data)+1):(2*nrow(KF.Rep.data))]
    cor.G<-as.matrix(cor.G)
    cor.G<-as.matrix(forceSymmetric(cor.G))
    if(length(rep.index)==nrow(inv.Rep)){
      inv.cor.G<-inv.Rep
    }else{
      inv.cor.G<-inv.Rep[rep.index,rep.index]-inv.Rep[rep.index,-rep.index,drop=F]%*%solve(inv.Rep[-rep.index,-rep.index])%*%inv.Rep[-rep.index,rep.index,drop=F]
    }
    inv.cor.G<-as.matrix(forceSymmetric(inv.cor.G))
    
    #S matrix
    diag_s<-KF.Rep.data[,1:nrow(KF.Rep.data)]
    
    #update S
    #temp.c<-min(eigen((M+1)/M*cor.G-diag_s)$values)
    #diag(diag_s)<-diag(diag_s)+max(temp.c-1e-5,0)
    
    #Knockoff preparation
    Sigma<-cor.G
    SigmaInv<-inv.cor.G
    P.each<-diag(1,ncol(Sigma))-diag_s%*%SigmaInv
    Sigma_k<-2*diag_s - diag_s%*%SigmaInv%*%diag_s
    #new trick by Jiaqi to improve computing time
    L1.left<-t(chol(Sigma_k-(M-1)/M*diag_s))
    Normal_1<-L1.left%*%matrix(rnorm(ncol(L1.left)),ncol(L1.left),1)
    L2.left<-t(chol(diag_s))
    temp<-as.matrix(L2.left%*%matrix(rnorm(ncol(L2.left)*M),ncol(L2.left),M))
    Normal_2<-temp-apply(temp,1,mean)
    Normal_k<-as.matrix(as.vector(Normal_1)+Normal_2)
    
    #Knockoff Z-score
    GK.Rep.Zscore_0<-as.matrix(KF.Rep.Z)
    GK.Rep.Zscore_k<-as.vector(P.each%*%GK.Rep.Zscore_0)+Normal_k
    colnames(GK.Rep.Zscore_k)<-paste0('GK.Zscore_',1:ncol(GK.Rep.Zscore_k))
    
    ## Generate knockoff Z-scores for the remaining variants
    cor.All.G<-KF.All.data[,(nrow(KF.All.data)+1):(2*nrow(KF.All.data)),drop=F]
    cor.All.G<-as.matrix(cor.All.G)
    cor.All.G<-as.matrix(forceSymmetric(cor.All.G))
    All.diag_s<-forceSymmetric(as.matrix(KF.All.data[,1:nrow(KF.All.data)]))
    
    temp.index<-which(KF.info[,'ifRep']==1)
    if(length(temp.index)!=nrow(cor.All.G)){
      Sigma11<-cor.All.G[temp.index,temp.index]
      Sigma22<-cor.All.G[-temp.index,-temp.index]
      Sigma12<-cor.All.G[temp.index,-temp.index]
      
      #random part of knockoff - remaining variables
      C<-Sigma22-t(Sigma12)%*%inv.cor.G%*%Sigma12
      C.left<-t(chol(C))
      E2<-as.matrix(C.left%*%matrix(rnorm(ncol(C.left)*M),ncol(C.left),M))
      E1<-t(Sigma12)%*%inv.cor.G
      GK.Rem.Zscore_k<-E1%*%GK.Rep.Zscore_k+E2
      
      #All z-scores
      GK.Zscore_0<-as.matrix(KF.Z)
      GK.Zscores_k<-matrix(0,nrow(GK.Zscore_0),M)
      GK.Zscores_k[temp.index,]<-GK.Rep.Zscore_k
      GK.Zscores_k[-temp.index,]<-GK.Rem.Zscore_k
    }else{
      GK.Zscore_0<-as.matrix(KF.Z)
      GK.Zscores_k<-GK.Rep.Zscore_k
    }
    colnames(GK.Zscores_k)<-paste0('GK.Zscore_',1:ncol(GK.Zscores_k))
    
    ## marginal
    temp.results<-cbind(KF.info,GK.Zscore_0,GK.Zscores_k)
    temp.results[,'clusters']<-paste0(chr,':',gsub('end','',gsub('start','',sub('.csv.*','',sub(".*Info_", "", KF.info.names[set])))),'-',temp.results[,'clusters'])
    results<-temp.results
    
    ## variant level
    T_0<-results[,'GK.Zscore_0',drop=F]^2
    T_k<-results[,grep('GK.Zscore_1',colnames(temp.results)):grep('GK.Zscore_5',colnames(temp.results)),drop=F]^2
    MK.stat<-MK.statistic(T_0,T_k)
    variant.kappa.marginal<-MK.stat[,'kappa']
    variant.tau.marginal<-MK.stat[,'tau']
    variant.W.marginal<-(variant.kappa.marginal==0)*variant.tau.marginal
    
    ## group level
    clusters<-results[,'clusters']
    Get.group.stat<-function(x){tapply(x,clusters,sum,na.rm=T)}
    T.group_0<-Get.group.stat(T_0[,1])
    T.group_k<-matrix(apply(T_k,2,Get.group.stat),ncol=ncol(T_k))
    MK.stat<-MK.statistic(T.group_0,T.group_k)
    group.kappa.marginal<-MK.stat[,'kappa']
    group.tau.marginal<-MK.stat[,'tau']
    group.W.marginal<-(group.kappa.marginal==0)*group.tau.marginal
    n.cluster<-as.numeric(table(clusters))
    p.marginal<-pnorm(abs(results[,'KF.Z']),lower.tail=F)*2
    s<-diag(as.matrix(KF.All.data[,1:nrow(KF.All.data)]))
    
    ##ghostbasil (split)
    temp<-results
    temp.id<-paste0(rep(paste0(temp[,'chr'],':',temp[,paste0('pos_',hg)],'-',temp[,'ref'],'-',temp[,'alt']),M+1),'_',
                    rep(0:M,each=nrow(temp)))
    id<-append(id,temp.id)
    
    temp.r<-as.vector(as.matrix(1/sqrt(N.effect)*temp[,grep('GK.Zscore_',colnames(temp))]))
    temp.z<-as.vector(as.matrix(temp[,grep('GK.Zscore_',colnames(temp))]))
    z<-append(z,temp.z)
    
    #add diagonal matrix for numerical performance
    alpha<-0.01
    D<-as.matrix(All.diag_s)
    S<-as.matrix(cor.All.G)
    S<-S-D
    D<-D+alpha*diag(1,nrow(D))
    
    D <- BlockMatrix(list(D))
    gmat <- BlockGroupGhostMatrix(S, D, M+1)
    mat.list <- append(mat.list, list(gmat))
    
    #zhaomeng validation
    #sigma part
    temp.sigma_stats<-t(as.matrix(temp.z[1:nrow(temp)]))%*%solve((1-shrinkage.lambda)*cor.All.G+diag(shrinkage.lambda,nrow(cor.All.G)))%*%as.matrix(temp.z[1:nrow(temp)])
    sigma_stats<-c(sigma_stats,temp.sigma_stats)
    n<-N.effect;p<-length(temp.z)
    sigma<-sqrt(min(1,max(0.5,(p+n+1)/(n+1)-sum(temp.sigma_stats)/(n+1))))
    #record heritability
    temp.heritability<-c(chr,start,end,shrinkage.lambda,1-sigma^2,lambda)
    heritability<-rbind(heritability,temp.heritability)
    
    #split basil
    #find beta corresponding to optimal lambda
    fit.basil<-ghostbasil(gmat, r=temp.z/sqrt(N.effect),user.lambdas=lambda.seq, delta.strong.size = 500, max.strong.size = nrow(gmat),n.threads=1,use.strong.rule=F)
    temp.beta<-fit.basil$betas[,ncol(fit.basil$betas)]
    #rescale beta so that different splits are comparable
    #temp.beta<-temp.beta/max(abs(temp.beta))*max(abs(KF.Z))
    
    T_0<-as.matrix(abs(temp.beta[grep('_0',temp.id)]))
    T_k<-abs(cbind(temp.beta[grep('_1',temp.id)],temp.beta[grep('_2',temp.id)],temp.beta[grep('_3',temp.id)],temp.beta[grep('_4',temp.id)],temp.beta[grep('_5',temp.id)]))
    MK.stat<-MK.statistic(T_0,T_k)
    variant.kappa.lassosplit<-MK.stat[,'kappa']
    variant.tau.lassosplit<-MK.stat[,'tau']
    variant.W.lassosplit<-(variant.kappa.lassosplit==0)*variant.tau.lassosplit
    
    #group level
    clusters<-results[,'clusters']
    Get.group.stat<-function(x){tapply(x,clusters,sum,na.rm=T)}
    T.group_0<-matrix(Get.group.stat(T_0[,1]),ncol=1)
    T.group_k<-matrix(apply(T_k,2,Get.group.stat),ncol=ncol(T_k))
    MK.stat<-MK.statistic(T.group_0,T.group_k)
    group.kappa.lassosplit<-MK.stat[,'kappa']
    group.tau.lassosplit<-MK.stat[,'tau']
    group.W.lassosplit<-(group.kappa.lassosplit==0)*group.tau.lassosplit
    
    ##susie
    A.each<-cor.All.G-All.diag_s
    A<-matrix(1,M+1,M+1)%x%A.each+diag(1,M+1)%x%All.diag_s
    #is_posdef(A)
    A<-Matrix(forceSymmetric(A))
    z.susie<-as.vector(as.matrix(results[,grep('GK.Zscore_0',colnames(results)):grep('GK.Zscore_5',colnames(results))]))
    fitted_rss <- try(suppressMessages(susie_rss(z=z.susie, R=as.matrix(A), L = min(10,nrow(results)))),silent=T)
    if(class(fitted_rss)!="try-error"){
      fitted_vars<-summary(fitted_rss)$vars
      beta<-fitted_vars[order(fitted_vars[,1]),2]#*z.susie^2
    }else{beta<-rep(0,length(z.susie));susie.error<-susie.error+1}
    ##variant level
    T_0<-as.matrix(abs(beta[1:nrow(results)]))
    T_k<-abs(matrix(beta[-(1:nrow(results))],nrow(results),M))
    MK.stat<-MK.statistic(T_0,T_k)
    variant.kappa.susie<-MK.stat[,'kappa']
    variant.tau.susie<-MK.stat[,'tau']
    variant.W.susie<-(variant.kappa.susie==0)*variant.tau.susie
    
    ##group level
    clusters<-results[,'clusters']
    Get.group.stat<-function(x){tapply(x,clusters,sum,na.rm=T)}
    T.group_0<-Get.group.stat(T_0[,1])
    T.group_k<-matrix(apply(T_k,2,Get.group.stat),ncol=ncol(T_k))
    MK.stat<-MK.statistic(T.group_0,T.group_k)
    group.kappa.susie<-MK.stat[,'kappa']
    group.tau.susie<-MK.stat[,'tau']
    group.W.susie<-(group.kappa.susie==0)*group.tau.susie
    
    #conventional susie
    z.susie<-as.vector(as.matrix(results[,grep('GK.Zscore_0',colnames(results))]))
    fitted_rss <- try(suppressMessages(susie_rss(z=z.susie, R=as.matrix(cor.All.G), L = min(10,nrow(results)))),silent=T)
    if(class(fitted_rss)!="try-error"){
      fitted_vars<-summary(fitted_rss)$vars
      susie.pip<-fitted_vars[order(fitted_vars[,1]),2]
      susie.cs<-fitted_vars[order(fitted_vars[,1]),3]
      susie.cs[susie.cs>=1]<-paste0(chr,':',gsub('end','',gsub('start','',sub('.csv.*','',sub(".*Info_", "", KF.info.names[set])))),'-',susie.cs[susie.cs>=1])
      susie.cs[susie.cs<1]<-NA
    }else{
      fitted_rss <- try(suppressMessages(susie_rss(z=z.susie, R=as.matrix((1-shrinkage.lambda)*cor.All.G+shrinkage.lambda*diag(nrow(cor.All.G))), L = min(10,nrow(results)))),silent=T)
      if(class(fitted_rss)!="try-error"){
        fitted_vars<-summary(fitted_rss)$vars
        susie.pip<-fitted_vars[order(fitted_vars[,1]),2]
        susie.cs<-fitted_vars[order(fitted_vars[,1]),3]
        susie.cs[susie.cs>=1]<-paste0(chr,':',gsub('end','',gsub('start','',sub('.csv.*','',sub(".*Info_", "", KF.info.names[set])))),'-',susie.cs[susie.cs>=1])
        susie.cs[susie.cs<1]<-NA
      }else{
        susie.pip<-rep(0,length(z.susie))
        susie.cs<-rep(NA,length(z.susie))
      }
    }
    
    group.results<-cbind(data.frame(names(T.group_0)),n.cluster,
                         group.kappa.marginal,group.tau.marginal,group.W.marginal,
                         group.kappa.susie,group.tau.susie,group.W.susie,
                         group.kappa.lassosplit,group.tau.lassosplit,group.W.lassosplit)
    group.results.all<-rbind(group.results.all,group.results)
    
    results<-cbind(results,p.marginal,susie.pip,susie.cs,s,
                   variant.kappa.marginal,variant.tau.marginal,variant.W.marginal,
                   variant.kappa.susie,variant.tau.susie,variant.W.susie,
                   variant.kappa.lassosplit,variant.tau.lassosplit,variant.W.lassosplit)
    results.all<-rbind(results.all,results)
    
    t2<-proc.time()
    
    print((t2-t1)[3])
    total.time<-total.time+(t2-t1)[3]
  }
}
print(total.time)

#genome-wide lasso fit via basil
A <- BlockBlockGroupGhostMatrix(mat.list)

#zhaomeng tunning
p<-length(z)
n<-N.effect
sigma<-sqrt(min(1,max(0.5,(p+n+1)/(n+1)-sum(sigma_stats)/(n+1))))

#estimate heritability
print(paste0('The estimated total heritability is ',1-sigma^2))
temp.heritability<-c(NA,NA,NA,NA,1-sigma^2,lambda)
heritability<-rbind(heritability,temp.heritability)
colnames(heritability)<-c('chr','start','end','LD.shrinkage.lambda','heritability','lasso.lambda')

temp.filename<-paste0(out.dir,'results_AD_Meta_All_heritability',replicate,'.txt')
write.table(heritability,temp.filename,col.names=T,row.names=F,sep='\t',quote=F)

#find beta corresponding to optimal lambda
fit.basil<-ghostbasil(A, r=z/sqrt(n),user.lambdas=lambda.seq, delta.strong.size = 500, max.strong.size = nrow(A),n.threads=1,use.strong.rule=F)
beta<-fit.basil$betas[,ncol(fit.basil$betas)]

##variant level
T_0<-as.matrix(abs(beta[grep('_0',id)]))
T_k<-abs(cbind(beta[grep('_1',id)],beta[grep('_2',id)],beta[grep('_3',id)],beta[grep('_4',id)],beta[grep('_5',id)]))
MK.stat<-MK.statistic(T_0,T_k)
variant.kappa.lasso<-MK.stat[,'kappa']
variant.tau.lasso<-MK.stat[,'tau']
variant.W.lasso<-(variant.kappa.lasso==0)*variant.tau.lasso
results.all<-cbind(results.all,
               variant.kappa.lasso,variant.tau.lasso,variant.W.lasso)

##group level
clusters<-results.all[,'clusters']
Get.group.stat<-function(x){tapply(x,clusters,sum,na.rm=T)}
T.group_0<-Get.group.stat(T_0[,1])
T.group_k<-apply(T_k,2,Get.group.stat)
MK.stat<-MK.statistic(T.group_0,T.group_k)
group.kappa.lasso<-MK.stat[,'kappa']
group.tau.lasso<-MK.stat[,'tau']
group.W.lasso<-(group.kappa.lasso==0)*group.tau.lasso
lasso.group.results<-cbind(data.frame(names(T.group_0)),group.kappa.lasso,group.tau.lasso,group.W.lasso)
group.results.all<-cbind(group.results.all,lasso.group.results[match(group.results.all[,1],lasso.group.results[,1]),-1])

results.all<-cbind(results.all,group.results.all[match(results.all[,'clusters'],group.results.all[,1]),-1])

#marginal q-value
group.q.marginal<-MK.q.byStat(group.results.all[,'group.kappa.marginal'],group.results.all[,'group.tau.marginal'],M=5)
group.q.marginal<-group.q.marginal[match(results.all[,'clusters'],group.results.all[,1])]
variant.q.marginal<-MK.q.byStat(results.all[,'variant.kappa.marginal'],results.all[,'variant.tau.marginal'],5,clusters=results.all[,'clusters'])

#susie q-value
group.q.susie<-MK.q.byStat(group.results.all[,'group.kappa.susie'],group.results.all[,'group.tau.susie'],M=5)
group.q.susie<-group.q.susie[match(results.all[,'clusters'],group.results.all[,1])]
variant.q.susie<-MK.q.byStat(results.all[,'variant.kappa.susie'],results.all[,'variant.tau.susie'],5,clusters=results.all[,'clusters'])

#lasso split q-value
group.q.lassosplit<-MK.q.byStat(group.results.all[,'group.kappa.lassosplit'],group.results.all[,'group.tau.lassosplit'],M=5)
group.q.lassosplit<-group.q.lassosplit[match(results.all[,'clusters'],group.results.all[,1])]
variant.q.lassosplit<-MK.q.byStat(results.all[,'variant.kappa.lassosplit'],results.all[,'variant.tau.lassosplit'],5,clusters=results.all[,'clusters'])

#lasso q-value
group.q.lasso<-MK.q.byStat(group.results.all[,'group.kappa.lasso'],group.results.all[,'group.tau.lasso'],M=5)
group.q.lasso<-group.q.lasso[match(results.all[,'clusters'],group.results.all[,1])]
variant.q.lasso<-MK.q.byStat(results.all[,'variant.kappa.lasso'],results.all[,'variant.tau.lasso'],5,clusters=results.all[,'clusters'])

#all results
results.all<-cbind(results.all,
                   variant.q.marginal,variant.q.susie,variant.q.lassosplit,variant.q.lasso,
                   group.q.marginal,group.q.susie,group.q.lassosplit,group.q.lasso)

results.all[results.all[,'variant.q.marginal']<=0.10,]
results.all[results.all[,'variant.q.susie']<=0.10,]
results.all[results.all[,'variant.q.lassosplit']<=0.10,]
results.all[results.all[,'variant.q.lasso']<=0.10,]

results.all[results.all[,'variant.q.lasso']<=0.10,c('KF.Z','p.marginal','variant.q.lasso')]

temp.filename<-paste0(out.dir,'results_AD_Meta_All_',replicate,'.txt')
write.table(results.all,temp.filename,col.names=T,row.names=F,sep='\t',quote=F)

print('done!')



#out.dir<-'/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/Results/'
#temp.filename<-paste0(out.dir,'results_AD_Meta_All_',replicate,'.txt')
#a<-fread(temp.filename)
#a<-data.frame(a)
#a[which(a[,'variant.q.lassosplit']<=0.1),]
#a[which(a[,'variant.q.lasso']<=0.1),]

