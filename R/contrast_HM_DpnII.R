#source(file.path(Sys.getenv("HOME"), "R", "source","myR.R"))
#source(file.path(Sys.getenv("HOME"), "R", "source","HiClib.R"))


options(stringsAsFactors = FALSE)
#library(ggplot2)
library(GenomicFeatures)
library(Matrix)

expand_hm_to_distance_matrix<-function(hm,nfrag,band_width){
  int_ud<-matrix(NA,nrow=nfrag,ncol=(2*nfrag-1))  #number of int up/down-stream
  for(i in 1:nfrag){
    start<-nfrag+1-i; end<-start+nfrag-1
    int_ud[i,start:end]<-hm[i,]
  }
  if(band_width==nfrag){
    return(int_ud)
  }else{
    NA_cols<-c( 1:(nfrag-band_width),(nfrag+band_width):(2*nfrag-1) )
                            #if(length(NA_cols)==2) NA_cols<-NA #dirty,but for consistent when band width=nfrag
                            #int_ud[,NA_cols]<-NA   #trim regions(columns) outside band width
    return(int_ud[,-NA_cols])
  }
}
clip_distance_matrix_to_hm<-function(int_ud,nfrag){
  hm<-matrix(NA,nrow=nfrag,ncol=nfrag)
  for(i in 1:nfrag){
    start<-nfrag+1-i; end<-start+nfrag-1
    hm[i,]<-int_ud[i,start:end]
  }
  return(hm)
}

#Samples<-c("AllCTR","AllCAPH2","AllSlimB","AllFlv","AllTrp")
Samples<-c("All4CTR")
resl<-"Frag"
Frag_count_cut<-60L#filter out fragments with too few reads;change with sample

ns<-25L#block summation size
#nc_inter<-ns*ns-1#number of valid entries for inter-region
nc_inter<-(ns-1)*(ns-1)#number of valid entries for inter-region
nc_intra<-ns*ns-ns#number of valid entries for intra-region
bw<-1000L#range to calculate raduie of gyration
setwd("/home/lli30/results/HiC/Domain/new_method")

if(resl=="Frag") {
  res_sit<-ReadRes(chrmajor,chrmajorlength,enz="DpnII_r603")
}

for(sample in Samples){

#sample<-Samples[1]

path<-paste("/compbio/projects/lli30/HiC/data",sample,sep="/")

vertical_ptg<-list()#percentage at vertical

nfrags<-rep(0L,length(chrmajor))

frag_sums<-list()

#======main===========
cntr_all<-interc_all<-intrac_all<-ncount_all<-list()

rg_all<-list()

for (ichr in seq_along(chrmajor)){
  #ichr=1#test only
  chrn<-chrmajor[ichr]
  cat(sample,chrn,"\n")
  FreqFile<-paste(path,"/",chrmajor[ichr],"_freq_table.txt",sep="") #use chr freqneucy table
  ft<-RT_easy(FreqFile)
  ft<-ft[ft[,2]>=1L & ft[,3]>=1L,]
  ft<-ft[abs(ft[,2]-ft[,3])>1L,]#remove adjacent cells
  ncount_all[[chrn]]<-sum(ft[,4])
  #ft_rev<-ft[,c(1,3,2,4)]
  #ft<-rbind(ft,ft_rev)
  nfrag<-length(res_sit[[chrn]])-1L
  hm<-sparseMatrix(c(ft[,2],ft[,3]), c(ft[,3],ft[,2]), x = c(ft[,4],ft[,4]),dims=c(nfrag,nfrag))
  #hist(log2(colSums(hm)+1),breaks=100)
  #hm<-hm[31:530,31:530]#test only
  #d<-nrow(hm)
  #rows<-1:(d-1)#remove adjacent cells
  #ind<-cbind(c(rows,rows+1), c(rows+1,rows))
  #hm[ind]<-0#remove adjacent ligation counts
  #diag(hm)<-0L
  frag_sum<-rowSums(hm,na.rm=T)
  hist(log2(frag_sum+1),100)
  frag_sums[[chrn]]<-c(frag_sum,NA)
  frag_sums[[chrn]]<-frag_sums[[chrn]]/(ncount_all[[chrn]]/1e7)
  ind_filter<-which(frag_sum<=Frag_count_cut)
  frag_sum[ind_filter]<-NA
  cind<-(-bw:bw)
  rg<-rep(NA,nfrag)
  for(ic in 1:nfrag){
    if(! ic %% 100 ) cat(ic," ")
    if(ic>bw & ic<(nfrag-bw)){
      if(ic %% 5000L ==1001L){#read sub-heatmap
        up<-min(nfrag,ic+5000L+1000L)
        hm_sub<-as.matrix(hm[(ic-1000L):up,(ic-1000L):up])
      }
      ind_local<-((ic-1001L) %% 5000L)+1001L
      rg[ic]<-weighted.mean(abs(cind),hm_sub[ind_local,(ind_local-bw):(ind_local+bw)]
                            ,na.rm=T)
    }else{#near two ends
      ind<-which((cind+ic)>0 & (cind+ic)<nfrag)
      rg[ic]<-weighted.mean(abs(cind[ind]),hm[ic,ind-bw+ic-1L],na.rm=T)
    }
  }
  #hm[ind_filter,]<-NA;
  #hm[,ind_filter]<-NA
  #hm_na_0<-hm
  #hm_na_0[ind_filter,]<-hm_na_0[,ind_filter]<-0
  #int_ud<-expand_hm_to_distance_matrix(hm_na_0,nfrag,bw+1)
  #int_ud[is.na(int_ud)]<-0
  #int_ud[ind_filter]<-0
  #rg<-int_ud %*% abs(-bw:bw) / (rowSums(int_ud)+0.1)
  #rg<-0
  rg[ind_filter]<-NA
  #hm_fr<-sweep(hm,1,frag_sum,"/")
  #hm_fc<-sweep(hm,2,frag_sum,"/")
  hm_f<-hm#hm or hm_fr or hm_fc or sqrt(hm_fr*hm_fc)
  cntr<-rep(NA,nfrag+1L)# number of terminal=nfrag+1
  interc<-intrac<-rep(NA,nfrag+1L)
  for(ic in (ns+1L):(nfrag-ns+1L)){
    #ic<-21L#test only
    if(! ic %% 500 ) cat(ic," ")
    ind_r<-ind_c<-(ic-ns):(ic+ns-1L)
    hms<-as.matrix(hm[ind_r,ind_c])
    fs<-frag_sum[ind_r]
    ind_na<-which(is.na(fs))
    hms[ind_na,]<-hms[,ind_na]<-NA
    region_inter<-c(hms[1:(ns-5),-(1:(ns+5))])#skip 5 frag around the end for inter region
    region_intra<-c(hms[1:ns,1:ns],hms[-(1:ns),-(1:ns)])
    #sum_inter<-sum(hms[1:ns,-(1:ns)],na.rm=T)
    sum_inter<-sum(region_inter,na.rm=T)
    sum_intra<-sum(region_intra,na.rm=T)/2L#on both sides
    #nna_inter<-length(which(is.na(hms[1:ns,-(1:ns)])))
    nna_inter<-length(which(is.na(region_inter)))
    nna_intra<-(length(which(is.na(region_intra))))/2L
    cs_top_inter<-cumsum(sort(region_inter,decreasing=T))
    cs_top_intra<-cumsum(sort(region_intra,decreasing=T))/2L
    outlier_inter<-cs_top_inter[1]>sum_inter*0.4#proportion of outlier
    outlier_intra<-cs_top_intra[2]>sum_intra*0.4 || cs_top_intra[4]>sum_intra*0.5
    if(nna_inter<=nc_inter/2 && nna_intra <=nc_intra/2
       &&(!outlier_inter) &&(!outlier_intra)){#enough non-0 entries, also no single/double outlier
      ratio<-sum_intra/sum_inter*((nc_inter-nna_inter)/nc_inter)*(nc_intra/(nc_intra-nna_intra))
      cntr[ic]<-ifelse(sum_intra >=Frag_count_cut,ratio,NA)
      interc[ic]<-sum_inter
      intrac[ic]<-sum_intra
    }
  }
  cntr_all[[chrn]]<-cntr;interc_all[[chrn]]<-interc;intrac_all[[chrn]]<-intrac
  rg_all[[chrn]]<-c(rg,NA)
}

###output in fragment unit
chr_f<-rep(names(cntr_all),unlist(lapply(cntr_all,length)))
bin_ind<-unlist(lapply(cntr_all,function(x) 1:length(x)))
cntr_v<-unlist(cntr_all)
cntr_v<-cntr_v/median(cntr_v,na.rm=T)#normalize
outd<-cbind(chr_f,bin_ind,unlist(res_sit)+1L,round(cntr_v,3),unlist(intrac_all),unlist(interc_all)
            ,round(unlist(rg_all),2),round(unlist(frag_sums),2))
colnames(outd)<-c("chr","ind","pos","cntr","intrac","interc","rg","frag_count")
outfn<-paste("contrastnew_adjremoved_outlieremoved_",sample,"_ns_",ns,"_intertrim5_bw_frag_withrg.txt",sep="")
WT_easy(outd,outfn,col.names=TRUE)

}#end loop over samples
