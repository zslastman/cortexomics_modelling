source(here::here('src/R/Rprofile.R'))
#rmarkdown::render(here('src/R/Modeling/run_degmodel_dropout.Rmd'))
dp_stanfile = here('src/Stan/mod_proDD.stan')%T>%{stopifnot(file.exists(.))}
dp_model = rstan::stan_model(dp_stanfile)



# library(doParallel)
# registerDoParallel(cores=2)

# library(doParallel)
# cl <- makeCluster(20)
# registerDoParallel(cl)
# foreach(i=1:3) %dopar% sqrt(i)

#load the metadata on all our genes and the many IDs involved
metainfo<-suppressMessages({read_tsv(here('data/metainfo.tsv'))})
#Pull out the gene names we want to analyze
uids4stan <- metainfo%>%
  filter(isbest)%>%#these are thee final pairs of gene/cds/mass spec ids that we use for modeling MS
  filter(sig_MS_change)%>%
  filter(n_stagemissing<=2)%>%#now filter them further - those with too many missing values will crash Rstan
  .$uprotein_id

#Make sure these interesting genes come first
testgenes <- metainfo%>%filter(gene_name%in%c('Flna','Satb2'))%>%filter(isbest)%>%.$uprotein_id%>%unique
uids4stan <- union(testgenes,uids4stan)

#
best_uprotein_ids <- uids4stan


c(mscountvoom) %<-% with(new.env(),{
	load(here('../cortexomics/data/integrate_exprdata2.Rdata'))
	list(mscountvoom)	
})

matchedms_mat<-readRDS(here('../cortexomics/data/matched_ms_mat.rds'))

inclusiontable(best_protein_ids,rownames(countvoom$E))

best_uprotein_ids <- metainfo%>%filter(isbest)%>%.$uprotein_id%>%unique
best_protein_ids <- best_uprotein_ids%>%str_replace('_\\d+$','')
rnamat <- mscountvoom$E%>%.[,str_subset(colnames(.),'total')]
rnamat %>% saveRDS('../cortexomics/data/rnamat.rds')
matchedms_mat %>% saveRDS('../cortexomics/data/matchedms_mat.rds')
rna_sigma <- countvoom$weights%>%{1/.}%>%
  set_rownames(rownames(countvoom$E))%>%
  set_colnames(colnames(countvoom$E))%>%
  .[best_protein_ids,str_subset(colnames(.),'total'),drop=F]%>%
  set_rownames(best_uprotein_ids)
rna_sigma %>% saveRDS('../cortexomics/data/rna_sigma.rds')
rnamat %>%saveRDS('../cortexomics/data/rnamat.rds')
# matchedms_mat%>%colMedians(na.rm=T)%>%txtplot


#this is the count data
rnamat <- (readRDS)(here('../cortexomics/data/rnamat.rds'))
#this contains the variance of the log count data
rna_sigma <- (readRDS)(here('../cortexomics/data/rna_sigma.rds'))

################################################################################
########Get prodd params 
################################################################################
pids4stan <- uids4stan%>%str_replace('_\\d+$','')
best_uprotein_ids <- uids4stan
best_protein_ids<-best_uprotein_ids%>%str_extract('[^_]+')

#this is the ms data
matchedms_mat <- (readRDS)(here('../cortexomics/data/matched_ms_mat.rds'))    

if(!file.exists(here('../cortexomics/data/rnaproddparams.Rdata'))){
  library(proDD)
  #get matrices of rnaseq and ms data, median norm them, then subtract again so the values center on zero
  #(the former is necessary, the second is simply for numeric reasons - stan will freak out if values are high)
    #ids of the ms data for these
  best_ms_ids <- metainfo%>%{.$ms_id[match(best_uprotein_ids,.$uprotein_id)]}
  #and of the rna data
  best_protein_ids <- metainfo%>%{.$protein_id[match(best_uprotein_ids,.$uprotein_id)]}

  matchedms_mat_rscl <- matchedms_mat[best_ms_ids,]
  matchedms_mat_rscl <- matchedms_mat_rscl%>%{proDD::median_normalization(.)}
  msmed <- matchedms_mat_rscl%>%apply(2,median,na.rm=T)%>%median(.,na.rm=T)
  matchedms_mat_rscl %<>% subtract(msmed)
  matchedms_mat_rscl%<>% set_rownames(best_uprotein_ids)
  
  tpvect<-matchedms_mat_rscl%>%colnames%>%str_extract('[^_]+')

  #Now fit our hyperparameters
  proddparam <- proDD::fit_parameters(matchedms_mat_rscl,tpvect)
  save(matchedms_mat_rscl,rnamatrscl,rnamed,msmed,proddparam,
       tpvect,rna_sigma,file='../cortexomics/data/rnaproddparams.Rdata')
}else{
  library(zeallot)
  # c(matchedms_mat_rscl,rnamatrscl,rnamed,msmed,proddparam,tpvect) %<-% get()
  # file.remove(file='data/proddparams.Rdata')
  load(here('../cortexomics/data/rnaproddparams.Rdata'))
}

matchedms_mat_rscl<-matchedms_mat_rscl[uids4stan,]

#align it all
rnamatrscl<-rnamatrscl[uids4stan,]
rna_sigma<-as.matrix(rna_sigma)%>%.[uids4stan,]%>%set_rownames(uids4stan)


rnamed <- rnamat[]%>%median(.,na.rm=T)
rnamatrscl <- rnamat[best_uprotein_ids,]%>%
  set_rownames(best_uprotein_ids)
rnamatrscl <- rnamatrscl%>%{proDD::median_normalization(.)}
rnamed <- rnamatrscl%>%apply(2,median,na.rm=T)%>%median(.,na.rm=T)
rnamatrscl %<>% subtract(rnamed)


################################################################################
########Function to serve data to stan
################################################################################
  

#function that pulls out the stan data we need in the right shape, and pairs it with proDD parameters etc. 
get_dp_standata <- function(sel_uprotein_id,
                            matchedms_mat_rscl,ribomatrscl,ribo_sigma,
                            params=proddparam
){
  stopifnot(rownames(matchedms_mat_rscl)==rownames(ribomatrscl))
  stopifnot(all(sel_uprotein_id%in%rownames(matchedms_mat_rscl)))
  # msids <- ms_id2protein_id%>%filter(uprotein_id %in% sel_uprotein_id)%>%.$ms_id
  msdata <- matchedms_mat_rscl[sel_uprotein_id,,drop=F]
  # msdata <- sizefactnorm(msdata)
  # MSmed = 0
  # MSmed = median(msdata,na.rm=T)
  # msdata = msdata - MSmed
  msdata %<>% replace_na(-Inf)
  n_missing_ms <- sum(!is.finite(msdata))
  #
  tpvect <- colnames(msdata)%>%str_extract('[^_]+')%>%as.factor%>%as.numeric
  #clear unecessary columns from this object
  ribotpvect <- colnames(ribomatrscl)%>%str_extract('[^_]+')%>%as.factor%>%as.numeric
  ribomat <- ribomatrscl[sel_uprotein_id,,drop=F]
  stopifnot(ncol(ribomat)==length(ribotpvect))
  # ribomed = median(ribomat,na.rm=T)
  # ribomat <- ribomat - ribomed
  #
  prodd_params <- params
  #
  require(splines2)
  time=1:n_distinct(tpvect)
  timeknots <- time[c(-1,-length(time))]
  mybs <- cbind(1,ibs(time, knots = timeknots,degree = 1, intercept = TRUE))
  mydbs = bs(time, knots = timeknots,degree = 1, intercept = TRUE)
  #
  standata = list(
    nsamples=ncol(msdata),#number of samples
    nribosamples=ncol(ribomat),#number of samples
    G=nrow(msdata),#number of proteins
    T=n_distinct(tpvect),#info on the number of conditions
    totalmissing=n_missing_ms,#info on th total amount of missing data
    lMS=msdata,#data
    lribo=ribomat,#data
    experimental_design=tpvect, #indicator variable matching the samples to conditions
    experimental_design_r=ribotpvect,
    zeta=prodd_params$hyper_params$zeta,#the spread of the dropout point for a library, gets combined with the varianc per protein
    rho=prodd_params$hyper_params$rho,#rho the location of the dropout point for a given library
    mu0=prodd_params$hyper_params$mu0,#fit by proDD the mean of means
    sigma20=prodd_params$hyper_params$sigma20,#fit by proDD - the variance in means
    eta=prodd_params$hyper_params$eta,#fit by proDD - th evariance in protein variances
    nu=prodd_params$hyper_params$nu,#fit by proDD - the mean of protein variances
    ribo_sigma=ribo_sigma[sel_uprotein_id,,drop=F],
    mybs = mybs,
    mydbs = mydbs
  )
  standata
  invisible(standata)
}

get_dp_standata <- partial(get_dp_standata,
                           matchedms_mat_rscl=matchedms_mat_rscl,
                           params=proddparam)

gnm2uid<-metainfo%>%filter(!is.na(uprotein_id),(isbest))%>%distinct(gene_name,uprotein_id)%>%{safe_hashmap(.[[1]],.[[2]])}

edist <- dist(rnamatrscl)

satb2like_uids <- as.matrix(edist)[gnm2uid[['Satb2']],]%>%sort%>%head(10)%>%names


get_dp_standata_withpriors <- function(id,...){
  c(get_dp_standata(id,...),list(
    l_pihalf_priorsd=3,
    l_pihalf_priormu=0.5,
    l_st_priorsd = 3,
    l_ribo_priorsd = 3,
    l_ribo_priornu = 2.7
    )
  )
}

stopifnot('l_ribo_priornu' %in% names(get_dp_standata_withpriors('ENSMUSP00000110057_4528',ribomatrscl=rnamatrscl,ribo_sigma=rna_sigma)))
stopifnot(get_dp_standata_withpriors(rownames(rnamatrscl)[1:3],ribomatrscl=rnamatrscl,ribo_sigma=rna_sigma)$lribo%>%dim%>%identical(c(3L,10L)))
stopifnot(rna_sigma[satb2like_uids,]%>%dim%>%identical(10L,10L))

get_dp_standata_withpriors(rownames(rnamatrscl)[1:3],ribomatrscl=rnamatrscl,ribo_sigma=rna_sigma)
get_dp_standata_withpriors(satb2like_uids,ribomatrscl=rnamatrscl,ribo_sigma=rna_sigma)


################################################################################
########Get initial values with opt
################################################################################
  
#Now optimize with many genes at once
# satb2like_uids <- as.matrix(edist)[gnm2uid[['Satb2']],]%>%sort%>%head(100)%>%names
#Random initializations don't work well, so we'lll 
satb2fits<-replicate(simplify=F,5,{
  pre = Sys.time()
  optres<-safely(rstan::optimizing)(dp_model,data=get_dp_standata_withpriors(gnm2uid[['Satb2']],ribomatrscl=rnamatrscl,rna_sigma),algorithm='Newton',as_vector=FALSE,hessian=TRUE)
 # optres<-(rstan::optimizing)(dp_model,data=get_dp_standata_withpriors(satb2like_uids,ribomatrscl=rnamatrscl,rna_sigma),algorithm='Newton',as_vector=FALSE,hessian=TRUE)
  Sys.time() - pre
  optres
})

#get the ones that didn't crash
satb2fits%<>%map('result')
message(str_interp("of ${length(satb2fits)} runs,  ${sum(satb2fits%>%map_lgl(is.null))} failed with an error"))
satb2fits%<>% .[!map_lgl(satb2fits,is.null)]
#get best fit results
bestoptes <- satb2fits[order(satb2fits%>%map_dbl('value'))%>%rev]%>%.[[1]]
#get initial values for single genes - set all to zero except the variance parameter
zeroinitvals <- bestoptes$par
for(p in names(zeroinitvals))zeroinitvals[[p]] <- zeroinitvals[[p]] - zeroinitvals[[p]]
zeroinitvals[['sigma2']] = zeroinitvals[['sigma2']]+0.2


counttype='ribo'
for(counttype in c('rna','ribo')){
  optfitlistfile <- str_interp('data/${counttype}_geneoptfits.rds')
  exprmat = if(counttype=='rna') rnamatrscl else ribomatrscl
  exprmat = if(counttype=='rna') rna_sigma else ribo_sigma
  if(!file.exists(optfitlistfile)){
  #
  require('R.utils')
  genefits <- list()
  #       
  genefits <- mclapply(mc.cores=20,uids4stan%>%setNames(.,.),(function(selgene){
    cat(paste0('.',which(uids4stan==selgene),'..'))
      #
      selgenefits<-lapply(rep(selgene,5),safely(function(selgene){
        withTimeout(timeout=15,{
          # message(selgene)
          optdat <- get_dp_standata_withpriors(selgene,rnamatrscl,rna_sigma)
          if(sum(is.finite(optdat$lMS))<2) return(NULL)
          #apparently we need actual hessian - Netwon succeeds for individual genes where
          #L-BFGS does not. 
          optimizing(dp_model,data=optdat,init=zeroinitvals,algorithm='Newton',as_vector=F,hessian=TRUE) 
          optres<-identity(rstan::optimizing)(dp_model,data=optdat,algorithm='Newton',as_vector=FALSE,hessian=TRUE)
        })
      }))
    selgenefits
  }))
  saveRDS(genefits,file=optfitlistfile)
}else{
  genefits<- readRDS((file=optfitlistfile))
}
}

genefits <- readRDS('data/rna_geneoptfits.rds')
genefits%>%setNames(uids4stan)%>%saveRDS('data/rna_geneoptfits.rds')



#Now optimize with many genes at once
satb2like_uids <- as.matrix(edist)[gnm2uid[['Satb2']],]%>%sort%>%head(100)%>%names
#Random initializations don't work well, so we'lll 
satb2_likefits<-replicate(simplify=F,5,{
  pre = Sys.time()
  optres<-safely(rstan::optimizing)(dp_model,data=get_dp_standata_withpriors(satb2like_uids,ribomatrscl=rnamatrscl,rna_sigma),algorithm='Newton',as_vector=FALSE,hessian=TRUE)
 # optres<-(rstan::optimizing)(dp_model,data=get_dp_standata_withpriors(satb2like_uids,ribomatrscl=rnamatrscl,rna_sigma),algorithm='Newton',as_vector=FALSE,hessian=TRUE)
  Sys.time() - pre
  optres
})


satb2like_uids%in%rownames(matchedms_mat_rscl)



stop()

mgsampling<-rstan::sampling(dp_model,data=get_dp_standata_withpriors(satb2like_uids,ribomatrscl=rnamatrscl,ribo_sigma=rna_sigma),chains=4,
                            init=function(){bestoptes$par},control=list(adapt_delta=.98,max_treedepth=20),iter=2e3)





stop()

dir.create(here('data/rnadpfits/'))
mclapply(mc.cores=20,rownames(matchedms_mat_rscl),function(uid){
	library(tidyverse)
	rdsfile <- str_interp(here('data/rnadpfits/${uid}_dpfit.rds'))
	if(!file.exists(rdsfile)){
	 	sampling<-rstan::sampling(dp_model,data=get_dp_standata_withpriors(uid),chains=4,
	                          init=function(){bestoptes$par},control=list(adapt_delta=.98,max_treedepth=15),iter=2e3)

	 	sampling%>%saveRDS(rdsfile)
	}
})

{
  require(BiocParallel)
#  BiocManager::install('batchtools')
  bpparam<- BatchtoolsParam(workers=1, cluster="slurm", resources=list(queue='all'))
  # sge_template <- '~/tmp.tmpl' 
  bpparam$template%>%readLines
  bpparam<- BatchtoolsParam(workers=1, cluster="slurm", resources=list(walltime = 4,queue='all',ncpus=4))
  #bplapply(BPPARAM=param
  #library(GenomicRanges)
  #result <- bplapply(BPPARAM=bpparam, c('1:1','1:3'), function(x)GenomicRanges::GRanges(x) )
}



#estimatesbak<-estimates
library(rstan)




estimates <- bplapply(BPPARAM=bpparam,uids4stan[1:2]%>%setNames(.,.),
  arglist=c('get_dp_standata_withpriors','rna_sigma','get_dp_standata','rnamatrscl','dp_model','matchedms_mat_rscl','proddparam')%>%setNames(.,.)%>%lapply(get),
  function(uid,arglist){
  attach(arglist)
  require(tidyverse)
  require(rstan)
  library(splines)
  library(purrr)
  library(magrittr)
  library(tidyverse)
  rdsfile <- str_interp(here('data/rnadpfits/${uid}_dpfit.rds'))
  if(!file.exists(rdsfile)){
    sampling<-rstan::sampling(dp_model,data=get_dp_standata_withpriors(uid),chains=4,
                            init=function(){bestoptes$par},control=list(adapt_delta=.98,max_treedepth=15),iter=2e3)

    
    if(!exists(here('data/samplingfile.rds'))){
      samplingfile <- create_samplingfile(dds, blind = TRUE)
      saveRDS(samplingfile,here('data/samplingfile.rds'))
    }else{
      samplingfile<-readRDS(here('data/samplingfile.rds'))
    }
})
estimates


allfitrds <- Sys.glob('data/dpfits/*')

lKDmeans <- allfitrds%>%mclapply(.%>%readRDS%>%summary%>%.[[1]]%>%.[,'50%']%>%.['lKd[1]'])
lKDmeans <- allfitrds%>%mclapply(.%>%readRDS%>%summary%>%.[[1]]%>%.[,'50%']%>%.['lKd[1]'])
lKDmeans%<>%setNames(allfitrds%>%basename%>%str_replace('_dpfit.rds',''))
lKDmeans%<>%keep(Negate(is.error))

lkdmeandf<-lKDmeans%>%map_dbl(1)%>%enframe('uprotein_id','lKd')

allcis <- allfitrds%>%mclapply(.%>%readRDS%>%summary%>%.[[1]]%>%.['lKd[1]',]%>%{.[c('2.5%','97.5%')]})
allcis%<>%setNames(allfitrds%>%basename%>%str_replace('_dpfit.rds',''))
allcis%<>%keep(Negate(is.error))
allcis%>%map_df(as.data.frame)

confuproteinids<-allcis%>%map_df(.id='uprotein_id',.%>%t%>%as.data.frame)%>%set_colnames(c('uprotein_id','lower','upper'))%>%mutate(width=upper-lower)%>%
	filter(width<3)%>%
	.$uprotein_id

mcshanedf<-fread('ext_data/mcshane_etal_2016_S1.csv')
#
mcshanethalfs<-mcshanedf%>%select(2,38,41)%>%set_colnames(c('gene_name','half_life','McShane_deg_cat'))
#
mcshanethalfs$half_life%<>%str_replace('> 300','300')%>%as.numeric
mcshanethalfs%<>%filter(half_life<299)
mcshanethalfs$half_life %<>% {./24}

mcshanethalfs%>%head
#now plot
plotfile<- here(paste0('plots/','rnalKD_indiv_v_mcshane_pihalf2','.pdf'))
dir.create(dirname(plotfile))
pdf(plotfile)
p = lkdmeandf%>%inner_join(metainfo%>%distinct(gene_name,uprotein_id))%>%inner_join(mcshanethalfs)%>%
	filter(McShane_deg_cat=='ED')%>%
	# filter(between(log2(half_life),-3,3))%>%
	# filter(uprotein_id%in%confuproteinids)%>%
	ggplot(.,aes(log2(half_life),lKd,color=McShane_deg_cat))+
	geom_point()+
	scale_color_discrete(name='McShane_deg_cat')+
	scale_x_continuous(paste0('log2(Half Life) (McShane et al)'))+
	scale_y_continuous(paste0('Estimated lKd'))+
	ggtitle(paste0('Measured Half Lives vs Estimated'))+
	geom_smooth(method='lm')+
	theme_bw()
p
dev.off()
normalizePath(plotfile)
cor.test(p$data$half_life,p$data$lKd,method='spearman')
cor.test(p$data$half_life,p$data$lKd)

################################################################################
########Now with ks
################################################################################


getwd()

fitlist <- list(
	'rna' = Sys.glob('data/rnadpfits/*'),
	'ribo' = Sys.glob('data/dpfits/*')
)


# fitlist%>%lapply(head,2)%>%lapply(.%>%setNames(.,basename(.)%>%str_replace('_(rna)?dpfit.rds','')))%>%map_df(.id='source',.%>%mclapply(.%>%readRDS%>%summary)
# stansumdata%<>%lapply(.%>%setNames(rnaallfitrds%>%basename%>%))
# stansumdata%<>%lapply(keep,Negate(is.error))

#read in the data from stan, summarising it as we go
if(!file.exists('data/stansumdata.rds')){
	stansumdata <- fitlist%>%
		lapply(.%>%setNames(.,basename(.)%>%str_replace('_(rna)?dpfit.rds','')))%>%
		# lapply(head)%>%
		lapply(.%>%mclapply(.%>%readRDS%>%summary%>%.[[1]]%>%as.data.frame%>%rownames_to_column('par'))%>%keep(Negate(is.error)))%>%
		map_df(.id='source',.%>%bind_rows(.id='uprotein_id'))
	stansumdata %<>% rename(estimate:='50%')
	stansumdata %<>% rename('lower':='2.5%')
	stansumdata %<>% rename('upper':='97.5%')

	stansumdata %>% saveRDS('data/stansumdata.rds')
}else{
	stansumdata <- readRDS('data/stansumdata.rds')
}
rename <- dplyr::rename

confuproteinids <- stansumdata%>%mutate(width=upper-lower)%>%
	filter(width<3)%>%
	.$uprotein_id%>%unique


################################################################################
########Get real data for comparison
################################################################################
	
countlinearTEs <- get_contrast_cis(
	bestonlycountebayes,
	t(countonly_pred_te_design['TE',,drop=F])
)%>%select(protein_id = uprotein_id,logFC,CI.L,CI.R,adj.P.Val)

countlinearTEs%<>%safe_left_join(metainfo%>%distinct(protein_id,gene_name),by='protein_id')

mcshanedf<-fread('ext_data/mcshane_etal_2016_S1.csv')
#
mcshanethalfs<-mcshanedf%>%select(2,38,41)%>%set_colnames(c('gene_name','half_life','McShane_deg_cat'))
#
mcshanethalfs$half_life%<>%str_replace('> 300','300')%>%as.numeric
mcshanethalfs%<>%filter(half_life<299)
mcshanethalfs$half_life %<>% {./24}

mcshanethalfs%>%head

library(txtplot)

isource='rna'
ipar='l_pihalf'

for(isource in names(fitlist)){
for(ipar in ipar ){

#now plot
plotfile<- here(paste0('plots/',isource,'_',ipar,'_indiv_v_mcshane_pihalf2','.pdf'))
dir.create(dirname(plotfile))
pdf(plotfile)
p = stansumdata%>%
	filter(source==isource)%>%
	filter(par%>%str_detect(ipar))%>%inner_join(metainfo%>%distinct(gene_name,uprotein_id))%>%inner_join(mcshanethalfs)%>%
	filter(McShane_deg_cat=='ED')%>%
	# filter(between(log2(half_life),-3,3))%>%
	# filter(uprotein_id%in%confuproteinids)%>%
	ggplot(.,aes(log2(half_life),estimate,color=McShane_deg_cat))+
	geom_point()+
	scale_color_discrete(name='McShane_deg_cat')+
	scale_x_continuous(paste0('log2(Half Life) (McShane et al)'))+
	scale_y_continuous(paste0('Estimated ',ipar,' individual model fits'))+
	ggtitle(paste0('Measured Half Lives vs Estimated'))+
	geom_smooth(method='lm')+
	theme_bw()
p
dev.off()
normalizePath(plotfile)
p$data%>%head
cor.test(p$data$half_life,p$data[,'estimate'],method='spearman')
cor.test(p$data$half_life,p$data[,'estimate'])

txtplot(p$data$half_life,p$data[,'estimate'],width=100)

}}
