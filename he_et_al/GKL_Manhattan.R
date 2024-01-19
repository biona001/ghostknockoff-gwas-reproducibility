# Ben's note: This file is essentially running GKL_Manhattan_GWAS.R first,
# then running GKL_Manhattan.R




#Enable command line arguments
args = commandArgs(TRUE)

replicate = as.numeric(args[1])
zscores_file = args[2]
input_dir = args[3]
ref_gene_file = args[4]
out_dir = args[5]

# for testing
# replicate = 11
# zscores_file = "/home/groups/sabatti/.julia/dev/ghostknockoff-gwas-reproducibility/data/AD_Zscores_Meta.txt"
# input_dir = "/home/groups/sabatti/.julia/dev/ghostknockoff-gwas-reproducibility/Results"
# ref_gene_file = "/home/groups/sabatti/.julia/dev/ghostknockoff-gwas-reproducibility/data/refGene_hg38.txt"
# out_dir = "/home/groups/sabatti/.julia/dev/ghostknockoff-gwas-reproducibility/Results"


library(data.table)
library(plyr)
library(dplyr)
library(CMplot)

#############################################################################
############ FIRST, CREATE MANHATTAN PLOT FOR MARGINAL ANALYSIS #############
############ In particular, this will generate an intermediate .csv file ####
############ needed for creating knockoff manhattan plots ###################
#############################################################################

################################################################################
## Main input
# input_dir<-'/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/Results/'
# out_dir<-'/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/Results/CMplot/'
input_file <- file.path(input_dir, paste0("results_AD_Meta_All_",replicate,".txt"))

study.names<-c('Jansen_etal','Kunkle_etal','Schwartzentruber_etal','Bis_etal',
               'Stanford2021_ADGC','Stanford2021_ADSP_WES','Stanford2021_ADSP_WGS',
               'LeGuan_etal','Huang_etal','Bellenguez_etal','Meta-analysis')

main.text<-paste0(study.names[replicate],'; marginal association test')
memo.text = paste0('MarginalAssociationTest_',study.names[replicate])

############
## read in hg38 gene list with position boundaries

#ref_seq <- read.table("/mnt/scratch/HARM_mar_2021/4d_supporting_files/ncbiRefSeq/ncbiRefSeqCurated.txt")
#colnames(ref_seq) <- c("bin",	"name", "chrom", "strand", "txStart", "txEnd",	"cdsStart",
#                       "cdsEnd", "exonCount",	"exonStarts", "exonEnds",	"score", "name2",	
#                       "cdsStartStat", "cdsEndStat", "exonFrames")

ref_seq<-read.table(ref_gene_file)
colnames(ref_seq)[3]<-'chrom'
colnames(ref_seq)[5:6]<-c('txStart','txEnd')
colnames(ref_seq)[13]<-c('name2')

## read result file
x1<- fread(input_file,header=T)
setnames(x1, c(3,7), c("CHR","BP"))
setnames(x1, 4:5, c("REF","ALT"))
x1 <- x1[, SNP:=paste(CHR,BP,REF,ALT,sep=":")]
x1 <- x1[which(!is.na(x1[,'CHR'])),]
setorderv(x1, c("CHR","BP"))
x1 <- as.data.frame(x1)
x1<-x1[match(unique(x1$SNP),x1$SNP),]

############ selected variants
x1_sug <- x1[x1$p.marginal<=5e-8,]
x1_sug <- x1_sug[!is.na(x1_sug$p.marginal),]

if(nrow(x1_sug)>2){
  
  #get individual z-scores
  Zscores<-as.data.frame(fread(zscores_file))
  temp.Zscores<-Zscores[match(x1_sug[,'rsid'],Zscores[,'rsid']),9:18]
  temp.Zscores<-cbind(temp.Zscores,pnorm(as.matrix(abs(temp.Zscores)),lower.tail=F))
  colnames(temp.Zscores)<-c(paste0('z.',study.names[1:10]),paste0('p.',study.names[1:10]))
  x1_sug<-cbind(x1_sug,temp.Zscores)

  x1_sug$IND <- NA
  x1_sug$TOP <- NA
  x1_sug$GENE <- NA
  x1_sug$RAst <- NA
  x1_sug$RAen <- NA
  
  loc_w <- 1000000
  ra_w <- 500000
  
  # find independent loci
  CHs <- unique(x1_sug$CHR)
  inds <- 0
  for (ch in CHs) {
    BPs <- x1_sug[which(x1_sug$CHR==ch),"BP"]
    BP_first <- BPs[1]
    for (bp in BPs) {
      if (bp==BP_first) {
        inds <- inds + 1
      } else {
        bp_sw <- bp - BP_prv
        if (bp_sw > loc_w) {
          inds <- inds + 1
        }
      }
      x1_sug[which(x1_sug$CHR==ch & x1_sug$BP==bp),"IND"] <- inds
      BP_prv <- bp
    }
  }
  
  # find top SNPs in independent loci
  for (i in x1_sug$IND) {
    x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","SNP","p.marginal")]
    x1_sug_t[which.max(-log10(x1_sug_t$p.marginal)),"TOP"] <- 
      x1_sug_t[which.max(-log10(x1_sug_t$p.marginal)),"SNP"]
    x1_sug[which(x1_sug$IND==i),c("TOP","SNP","p.marginal")] <- x1_sug_t
  }
  
  # set range +/- extension for independent loci
  for (i in x1_sug$IND) {
    x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")]
    x1_sug_t[,"RAst"] <- x1_sug_t[1,"BP"] - ra_w
    x1_sug_t[,"RAen"] <- x1_sug_t[dim(x1_sug_t)[1],"BP"] + ra_w
    x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")] <- x1_sug_t
  }
  
  # find gene next to top SNPs in independent loci
  for (t in x1_sug[which(!is.na(x1_sug$TOP)),"TOP"]) {
    ra <- x1_sug[which(x1_sug$TOP==t),c("CHR","BP","GENE")]
    ref_seq_t <- ref_seq[which(ref_seq$chrom==paste("chr",ra$CHR,sep="")),]
    st_dist <- min(abs(ref_seq_t$txStart - ra$BP))
    st_dist_w <- which.min(abs(ref_seq_t$txStart - ra$BP))
    en_dist <- min(abs(ref_seq_t$txEnd - ra$BP))
    en_dist_w <- which.min(abs(ref_seq_t$txEnd - ra$BP))
    gene <- ""
    if (st_dist<=1000000 | en_dist<=1000000) {
      if (st_dist<=en_dist) {
        gene <- ref_seq_t[st_dist_w,"name2"]
      } else {
        gene <- ref_seq_t[en_dist_w,"name2"]
      }
    }
    x1_sug[which(x1_sug$TOP==t),c("GENE")] <- as.character(gene)
  }
  
  #flag new loci
  new_loc_w <- 1000000
  x1_sug$NEW_hit <- "Y"
  #temp<-cbind(unique(x1_sug$IND),tapply(x1_sug$Schwartzentruber_etal_P,x1_sug$IND,min,na.rm=T),tapply(x1_sug$Schwartzentruber_etal_P,x1_sug$IND,min,na.rm=T)<= 5E-8)
  #x1_sug[which(x1_sug$Schwartzentruber_etal_P <= 5E-8),"NEW_hit"] <- "N"
  AD_loci<-x1_sug
  AD_loci<-AD_loci[AD_loci[,'p.marginal']<=5e-8,]
  temp.new<-c()
  for(i in 1:nrow(x1_sug)){
    temp.chr<-x1_sug[i,'CHR']
    temp.pos<-x1_sug[i,'BP']
    temp.new<-c(temp.new,sum((temp.chr==AD_loci[,'CHR']) & (temp.pos>=AD_loci[,'BP']-new_loc_w) & (temp.pos<=AD_loci[,'BP']+new_loc_w), na.rm=T))
  }
  for(i in 1:nrow(x1_sug)){
    if(sum(temp.new[which(x1_sug[,'IND']==x1_sug[i,'IND'])])>0){x1_sug[i,"NEW_hit"] <- "N"}
  }
  
  # extract signal from all variants for respective independent loci (for plotting)
  x1_sig<-x1_sug[x1_sug[,'p.marginal']<=5e-8,]
  signal <- c()
  for (t in x1_sig[which(!is.na(x1_sig$TOP)),"TOP"]) {
    ra <- x1_sig[which(x1_sig$TOP==t),c("CHR","RAst","RAen","GENE","p.marginal","IND")]
    indt <- ra$IND
    xs <- x1[which(x1$CHR==ra$CHR & (x1$BP >= ra$RAst) & (x1$BP <= ra$RAen) ),
             c("SNP","p.marginal","BP","CHR")]
    #  if (any(xs$BP > (rs429358_pos-1000000) & xs$BP < (rs429358_pos+1000000) & xs$CHR==19)) {
    #    xs <- xs[-which(xs$BP > (rs429358_pos-1000000) & 
    #                      xs$BP < (rs429358_pos+1000000) & xs$CHR==19),]
    #  }
    if (length(xs$SNP)>0) {
      xs <- xs[,c("SNP",'p.marginal')]
      if (ra$p.marginal<=5e-8) {xs$col <- "purple"}
      #if (ra$q<=0.05) {xs$col <- "red"}
      xs$text_col <- "red"; 
      if (any(x1_sig[which(x1_sig$IND==indt),"NEW_hit"] %in% "N")) {xs$text_col <- "red"}
      xs$pch <- 19
      #xs$text <- NA; xs[which.max(xs$W),"text"] <- ra$GENE
      xs$text <- NA; xs[,"text"] <- ra$GENE
      signal <- rbind(signal,xs)
    }
  }
  
  # keep signal with highlights
  signal_high <- signal[which(signal$p.marginal<=5e-8),]
  signal_top <- signal[which(signal$SNP %in% x1_sug[which(!is.na(x1_sug$TOP)),"TOP"]),]
  signal_topp <- signal_top
  
  # remove any empty gene names in signal_topp if applicable
  if (any(which(signal_topp$text %in% ""))) {
    signal_topp <- signal_topp[-which(signal_topp$text %in% ""),]
  }
  
  ############
  ## Save SNPs that pass suggestive significance
  
  # write into files
  write.csv(x1_sug, file.path(out_dir, paste(memo.text,".stats.sug_hits.info.csv",sep="")),
            row.names = F)
  #write.table(x1_sug$SNP, paste(in_fold,out_file,".stats.sug_hits.snpID.txt",sep=""), 
  #            col.names=F, row.names=F, quote=F)
  
  # thresholds
  #T010 <- min(x1[x1$q<=0.10,"W"]); T005 <- x1_sug[which.min(abs(x1_sug$q-0.05)),"W"]
  ths <- c(-log10(5e-8))
  
  ############
  # CMplot
  
  x1t <- x1[,c("SNP","CHR","BP",'p.marginal')]
  #x1t <- x1t[-which(x1t$BP > (rs429358_pos-1000000) & 
  #                    x1t$BP < (rs429358_pos+1000000) & x1t$CHR==19),]
  
  x1t[which(x1t$p.marginal<=1e-50),"p.marginal"] <- 1e-50
  #main.text<-'Meta-analysis of Alzheimer\'s disease via conventional marginal assocation test'
  
  
  setwd(out_dir)
  x1t <- x1t[,c("SNP","CHR","BP","p.marginal")]
  x1t[,4]<--log10(x1t[,4])
  #colnames(x1t) <- c("SNP","chr","pos",paste(sub('.txt','',input_file),'_suggestive',sep=""))
  CMplot(x1t, plot.type="m", LOG10=FALSE, col=c("grey30","grey60"), ylab="-log10(p)",bin.range=c(0,500),
         chr.den.col=c("darkgreen", "yellow", "red"), main=main.text,
         threshold=ths, threshold.lty=c(2), threshold.lwd=c(1), threshold.col=c("black"),
         highlight=signal_topp$SNP, highlight.cex=1, 
         highlight.col=signal_topp$col, highlight.text.col=signal_topp$text_col, highlight.text=signal_topp$text,
         signal.col=c("cornflowerblue"),signal.cex=c(1),
         file="jpg",file.name=memo.text,dpi=300,file.output=TRUE,verbose=TRUE,width=20,height=6)
  
}






############################################################################
############ Next, CREATE MANHATTAN PLOT FOR KNOCKOFF ANALYSIS #############
############################################################################

# library(data.table)
# library(plyr)
# library(dplyr)
# library(CMplot)
# library(foreach)
# library(doParallel)

# source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
# source('/oak/stanford/groups/zihuai/ESGWAS/ESGWAS.R')
# source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')

Get.transfer.Tau.fullmodel.recursive<-function(kappa,tau,U,reveal_prop=0.5){
  revealed_id<-which(tau<=quantile(tau,reveal_prop))
  unrevealed_id<-(1:length(tau))[-revealed_id]
  U<-as.matrix(U)
  X<-cbind(1,tau,U)
  C<-X%*%solve(t(X)%*%X)
  max.ini.tau<-max(tau[revealed_id])
  transfer.tau<-tau
  for(k in 1:length(unrevealed_id)){
    print(k)
    if(k==1){
      temp.y<-rep(0,length(tau))
      temp.y[revealed_id] = as.numeric(kappa[revealed_id]!=0)
      new.pred.y<-rep(NA,length(tau))
      new.pred.y[unrevealed_id]<-C[unrevealed_id,,drop=F]%*%(t(X)%*%temp.y)
    }else{
      if(temp.y[add_id]!=0){
        new.pred.y[unrevealed_id]<-old.pred.y[unrevealed_id]+C[unrevealed_id,,drop=F]%*%t(X[add_id,,drop=F])
      }
    }
    add_id<-unrevealed_id[which.max(new.pred.y[unrevealed_id])]
    transfer.tau[add_id]<-max.ini.tau+k
    temp.y[add_id]<-as.numeric(kappa[add_id]!=0)
    revealed_id<-c(revealed_id,add_id)
    unrevealed_id<-(1:length(tau))[-revealed_id]
    old.pred.y<-new.pred.y
  }
  return(transfer.tau)
}


# SummaryStat.info<-read.table('/oak/stanford/groups/zihuai/ESGWAS_lasso/AD_Analysis/SummaryStatInfo.txt',header=T,sep='\t')
study.names<-c('Jansen_etal','Kunkle_etal','Schwartzentruber_etal','Bis_etal',
               'Stanford2021_ADGC','Stanford2021_ADSP_WES','Stanford2021_ADSP_WGS',
               'LeGuan_etal','Huang_etal','Bellenguez_etal','Meta-analysis')

for(k in 1:2){
  print(k)
  for(index in 1:4){
    print(index)
    ################################################################################
    ## Main input
    #input_dir<-'/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/Results_Before2023July28/'
    #out_dir<-'/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/Results_Before2023July28/CMplot/'
    # input_dir<-'/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/Results/'
    # out_dir<-'/oak/stanford/groups/zihuai/ESGWAS_lasso/Z_scores/AD_meta/Results/CMplot/'
    
    Zscores<-as.data.frame(fread(zscores_file))
    ref_seq<-read.table(ref_gene_file)

    #input_file <- paste0("results_AD_Meta_All_",replicate,"_susieCSadded.txt")
    input_file <- paste0("results_AD_Meta_All_",replicate,".txt")
    x1.all <- fread(file.path(input_dir,input_file),header=T)
    
    if(k==1){type<-'variant'}
    if(k==2){type<-'group'}
    if(k==3){type<-'variantside'}
    
    if(index==1){main.text<-paste0(study.names[replicate],'; conditional independent test CIT-Marginal');memo.text=paste0('GhostKnockoffMarginal_',type,'_',study.names[replicate])}
    if(index==2){main.text<-paste0(study.names[replicate],'; conditional independent test CIT-Lasso');memo.text=paste0('GhostKnockoffLasso_',type,'_',study.names[replicate])}
    if(index==3){main.text<-paste0(study.names[replicate],'; conditional independent test CIT-LassoSplit');memo.text=paste0('GhostKnockoffLassoSplit_',type,'_',study.names[replicate])}
    if(index==4){main.text<-paste0(study.names[replicate],'; conditional independent test CIT-Susie');memo.text=paste0('GhostKnockoffSusie_',type,'_',study.names[replicate])}
    
    ############
    ## read in hg38 gene list with position boundaries
    
    #ref_seq <- read.table("/mnt/scratch/HARM_mar_2021/4d_supporting_files/ncbiRefSeq/ncbiRefSeqCurated.txt")
    #colnames(ref_seq) <- c("bin",	"name", "chrom", "strand", "txStart", "txEnd",	"cdsStart",
    #                       "cdsEnd", "exonCount",	"exonStarts", "exonEnds",	"score", "name2",	
    #                       "cdsStartStat", "cdsEndStat", "exonFrames")
    
    colnames(ref_seq)[3]<-'chrom'
    colnames(ref_seq)[5:6]<-c('txStart','txEnd')
    colnames(ref_seq)[13]<-c('name2')
    
    ############
    
    ## read result file
    x1<-x1.all
    setnames(x1, c(3,7), c("CHR","BP"))
    setnames(x1, 4:5, c("REF","ALT"))
    x1 <- x1[, SNP:=paste(CHR,BP,REF,ALT,sep=":")]
    x1 <- x1[which(!is.na(x1[,'CHR'])),]
    setorderv(x1, c("CHR","BP"))
    x1 <- as.data.frame(x1)
    x1 <- x1 [,colnames(x1)!='p.marginal']
    temp<-c('.marginal','.lasso','.lassosplit','.susie')
    colnames(x1)[match(paste0(c('variant.kappa','variant.tau','variant.W',
                                'group.kappa','group.tau','group.W',
                                'variant.q','group.q'),temp[index]),
                       colnames(x1))]<-c('variant.kappa','variant.tau','variant.W',
                                         'group.kappa','group.tau','group.W',
                                         'variant.q','group.q')
    #x1<-x1[,-c(grep(temp[-index][1],colnames(x1)),grep(temp[-index][2],colnames(x1)))]
    #find unique variants
    x1<-x1[match(unique(x1$SNP),x1$SNP),]
    
    #side info
    if(k==3){
      U<-x1[,c('CADD','DANN','FATHMM-MKL','FunSeq2','Eigen','GenoCanyon','FIRE','LINSIGHT')]
      U[is.na(U)]<--10*log10(0.5) #median of pred score
      variantside.tau<-Get.transfer.Tau.fullmodel.recursive(x1[,'variant.kappa'],x1[,'variant.tau'],U,reveal_prop=0.5)
      variantside.tau<-variantside.tau/max(variantside.tau)*max(x1[,'variant.tau'])
      variantside.kappa<-x1[,'variant.kappa']
      variantside.W<-(variantside.kappa==0)*variantside.tau
      variantside.q<-MK.q.byStat(x1[,'variant.kappa'],variantside.tau,5,clusters=x1[,'clusters'],Rej.Bound=10000)
      x1<-cbind(x1,variantside.kappa,variantside.tau,variantside.W,variantside.q)
    }
    
    ############ selected variants
    x1_sug <- x1[x1[,paste0(type,'.q')]<=0.10,,drop=F]
    x1_sug <- x1_sug[!is.na(x1_sug[,paste0(type,'.q')]),,drop=F]
    print(nrow(x1_sug))
    if(nrow(x1_sug)<=2){next}
    
    # calculate P-vals for cohort data and add flag for at least one genome wide sig
    P<-2*pnorm(as.matrix(abs(x1_sug[,grep('KF.Z',colnames(x1_sug))])),lower.tail = F)
    colnames(P)<-'Marginal_P'
    x1_sug<-cbind(x1_sug,P)
    
    #get individual z-scores
    temp.Zscores<-Zscores[match(x1_sug[,'rsid'],Zscores[,'rsid']),9:18]
    temp.Zscores<-cbind(temp.Zscores,pnorm(as.matrix(abs(temp.Zscores)),lower.tail=F))
    colnames(temp.Zscores)<-c(paste0('z.',study.names[1:10]),paste0('p.',study.names[1:10]))
    x1_sug<-cbind(x1_sug,temp.Zscores)
    
    x1_sug$IND <- NA
    x1_sug$TOP <- NA
    x1_sug$GENE <- NA
    x1_sug$RAst <- NA
    x1_sug$RAen <- NA
    
    loc_w <- 1000000
    ra_w <- 500000
    
    # find independent loci
    CHs <- unique(x1_sug$CHR)
    inds <- 0
    for (ch in CHs) {
      BPs <- x1_sug[which(x1_sug$CHR==ch),"BP"]
      BP_first <- BPs[1]
      for (bp in BPs) {
        if (bp==BP_first) {
          inds <- inds + 1
        } else {
          bp_sw <- bp - BP_prv
          if (bp_sw > loc_w) {
            inds <- inds + 1
          }
        }
        x1_sug[which(x1_sug$CHR==ch & x1_sug$BP==bp),"IND"] <- inds
        BP_prv <- bp
      }
    }
    
    # find top SNPs in independent loci
    for (i in x1_sug$IND) {
      x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","SNP",paste0('variant','.W'))]
      x1_sug_t[which.max(x1_sug_t[,paste0('variant','.W')]),"TOP"] <- 
        x1_sug_t[which.max(x1_sug_t[,paste0('variant','.W')]),"SNP"]
      x1_sug[which(x1_sug$IND==i),c("TOP","SNP",paste0('variant','.W'))] <- x1_sug_t
    }
    
    # set range +/- extension for independent loci
    for (i in x1_sug$IND) {
      x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")]
      x1_sug_t[,"RAst"] <- x1_sug_t[1,"BP"] - ra_w
      x1_sug_t[,"RAen"] <- x1_sug_t[dim(x1_sug_t)[1],"BP"] + ra_w
      x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")] <- x1_sug_t
    }
    
    # find gene next to top SNPs in independent loci
    for (t in x1_sug[which(!is.na(x1_sug$TOP)),"TOP"]) {
      ra <- x1_sug[which(x1_sug$TOP==t),c("CHR","BP","GENE")]
      ref_seq_t <- ref_seq[which(ref_seq$chrom==paste("chr",ra$CHR,sep="")),]
      st_dist <- min(abs(ref_seq_t$txStart - ra$BP))
      st_dist_w <- which.min(abs(ref_seq_t$txStart - ra$BP))
      en_dist <- min(abs(ref_seq_t$txEnd - ra$BP))
      en_dist_w <- which.min(abs(ref_seq_t$txEnd - ra$BP))
      gene <- ""
      if (st_dist<=1000000 | en_dist<=1000000) {
        if (st_dist<=en_dist) {
          gene <- ref_seq_t[st_dist_w,"name2"]
        } else {
          gene <- ref_seq_t[en_dist_w,"name2"]
        }
      }
      x1_sug[which(x1_sug$TOP==t),c("GENE")] <- as.character(gene)
    }
    
    #flag new loci
    new_loc_w <- 1000000
    x1_sug$NEW_hit <- "Y"
    #x1_sug[x1_sug[,'q']<=0.10,'NEW_hit']<-'N'
    #temp<-cbind(unique(x1_sug$IND),tapply(x1_sug$Schwartzentruber_etal_P,x1_sug$IND,min,na.rm=T),tapply(x1_sug$Schwartzentruber_etal_P,x1_sug$IND,min,na.rm=T)<= 5E-8)
    #x1_sug[which(x1_sug$Schwartzentruber_etal_P <= 5E-8),"NEW_hit"] <- "N"
    
    AD_loci<-read.csv(file.path(out_dir, paste0('MarginalAssociationTest_',study.names[replicate],".stats.sug_hits.info.csv")))
    #AD_loci<-x1_sug
    #AD_loci<-AD_loci[AD_loci[,paste0(type,'.q')]<=0.10,]
    temp.new<-c()
    for(i in 1:nrow(x1_sug)){
      temp.chr<-x1_sug[i,'CHR']
      temp.pos<-x1_sug[i,'BP']
      temp.new<-c(temp.new,sum((temp.chr==AD_loci[,'CHR']) & (temp.pos>=AD_loci[,'BP']-new_loc_w) & (temp.pos<=AD_loci[,'BP']+new_loc_w), na.rm=T))
    }
    for(i in 1:nrow(x1_sug)){
      if(sum(temp.new[which(x1_sug[,'IND']==x1_sug[i,'IND'])])>0){x1_sug[i,"NEW_hit"] <- "N"}
    }
    
    # extract signal from all variants for respective independent loci (for plotting)
    x1_sig<-x1_sug[x1_sug[,paste0(type,'.q')]<=0.20,]
    signal <- c()
    for (t in x1_sig[which(!is.na(x1_sig$TOP)),"TOP"]) {
      ra <- x1_sig[which(x1_sig$TOP==t),c("CHR","RAst","RAen","GENE",paste0(type,'.q'),paste0(type,'.W'),"IND")]
      indt <- ra$IND
      xs <- x1[which(x1$CHR==ra$CHR & (x1$BP >= ra$RAst) & (x1$BP <= ra$RAen) ),
               c("SNP",paste0(type,'.q'),paste0(type,'.W'),"BP","CHR")]
      #  if (any(xs$BP > (rs429358_pos-1000000) & xs$BP < (rs429358_pos+1000000) & xs$CHR==19)) {
      #    xs <- xs[-which(xs$BP > (rs429358_pos-1000000) & 
      #                      xs$BP < (rs429358_pos+1000000) & xs$CHR==19),]
      #  }
      if (length(xs$SNP)>0) {
        xs <- xs[,c("SNP",paste0(type,'.W'),paste0(type,'.q'))]
        if (ra[,paste0(type,'.q')]<=0.20) {xs$col <- "purple"}
        #if (ra$q<=0.05) {xs$col <- "red"}
        xs$text_col <- "blue";
        if (any(x1_sig[which(x1_sig$IND==indt),"NEW_hit"] %in% "N")) {xs$text_col <- "red"}
        xs$pch <- 19
        #xs$text <- NA; xs[which.max(xs$W),"text"] <- ra$GENE
        xs$text <- NA; xs[,"text"] <- ra$GENE
        signal <- rbind(signal,xs)
      }
    }
    
    # keep signal with highlights
    signal_high <- signal[which(signal[,paste0(type,'.q')]<=0.10),]
    signal_top <- signal[which(signal$SNP %in% x1_sug[which(!is.na(x1_sug$TOP)),"TOP"]),]
    signal_topp <- signal_top
    
    # remove any empty gene names in signal_topp if applicable
    if (any(which(signal_topp$text %in% ""))) {
      signal_topp <- signal_topp[-which(signal_topp$text %in% ""),]
    }
    
    ############
    ## Save SNPs that pass suggestive significance
    
    # write into files
    colnames(x1_sug)[ncol(x1_sug)]<-'NewHits'
    write.csv(x1_sug, file.path(out_dir, paste(memo.text,".stats.sug_hits.info.csv",sep="")),
              row.names = F)
    
    #write.table(x1_sug$SNP, paste(in_fold,out_file,".stats.sug_hits.snpID.txt",sep=""), 
    #            col.names=F, row.names=F, quote=F)
    
    # thresholds
    T020 <- min(x1[x1[,paste0(type,'.q')]<=0.20,paste0(type,'.W')]); T010 <- min(x1[x1[,paste0(type,'.q')]<=0.10,paste0(type,'.W')])#x1_sug[which.min(abs(x1_sug$q-0.10)),"W"]
    #T020 <- min(x1[x1$q<=0.05,"W"])
    #T010 <- -log10(0.10); T005 <- -log10(0.05)
    #ths <- c(T010,T005)
    
    ############
    # CMplot
    
    x1t <- x1[,c("SNP","CHR","BP",paste0(type,'.W'),paste0(type,'.q'))] 
    #x1t <- x1t[-which(x1t$BP > (rs429358_pos-1000000) & 
    #                    x1t$BP < (rs429358_pos+1000000) & x1t$CHR==19),]
    x1t[which(x1t[,paste0(type,'.W')]>T010*8),paste0(type,'.W')]<- T010*8;ylim=c(0,T010*8)
    #if(index==1){
    #  x1t[which(x1t[,paste0(type,'.W')]>200),paste0(type,'.W')] <- 200;ylim=NULL
    #}
    #if(index==2){
    #  x1t[which(x1t[,paste0(type,'.W')]>0.02),paste0(type,'.W')]<- 0.02;ylim=c(0,0.02)
    #}
    #if(index==3){
    #  x1t[which(x1t[,paste0(type,'.W')]>0.02),paste0(type,'.W')]<- 0.02;ylim=c(0,0.02)
    #}
    #if(index==4){
    #  x1t[which(x1t[,paste0(type,'.W')]>1),paste0(type,'.W')] <- 1;ylim=c(0,1.5)
    #}
    if(T010!=Inf){T010<-min(T010,max(x1t[,paste0(type,'.W')]))}
    if(T020!=Inf){T020<-min(T020,max(x1t[,paste0(type,'.W')]))}
    ths <- c(T010)#,T020)
    
    #x1t<-cbind(x1t,-log10(x1t[,'q']))
    #colnames(x1t)[ncol(x1t)]<-'-log10(q)'
    signal_topp<-signal_topp[signal_topp[,paste0(type,'.q')]<=0.1,]
    
    library(CMplot)
    setwd(out_dir)
    x1t <- x1t[,c("SNP","CHR","BP",paste0(type,'.W'))]
    colnames(x1t)[4]<-'Test statistic'
    #colnames(x1t) <- c("SNP","chr","pos",paste(sub('.txt','',input_file),'_suggestive',sep=""))
    CMplot(x1t, plot.type="m", LOG10=FALSE, col=c("grey30","grey60"), ylab="Test statistic", ylim=ylim,bin.range=c(0,500),
           chr.den.col=c("darkgreen", "yellow", "red"), main=main.text,
           threshold=ths, threshold.lty=c(2), threshold.lwd=c(1), threshold.col=c("red"),
           highlight=signal_topp$SNP, highlight.cex=1, 
           highlight.col=signal_topp$col, highlight.text.col=signal_topp$text_col, highlight.text=signal_topp$text,
           signal.col=c("cornflowerblue"),signal.cex=c(1),
           file="jpg",file.name=memo.text,dpi=300,file.output=TRUE,verbose=TRUE,width=20,height=6)
    
  }
}




