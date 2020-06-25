source(here::here('src/R/Rprofile.R'))

#rmarkdown::render(here('src/R/Modeling/run_degmodel_dropout.Rmd'))
dp_stanfile = here('src/Stan/mod_proDD.stan')%T>%{stopifnot(file.exists(.))}
dp_model = rstan::stan_model(dp_stanfile)


require(tidyverse)
require(rstan)
library(splines)
library(purrr)
library(magrittr)
library(here)

################################################################################
########Functions to serve data to stan
################################################################################
#load the metadata on all our genes and the many IDs involved
metainfo<-suppressMessages({read_tsv(here('data/metainfo.tsv'))})
#Pull out the gene names we want to analyze
uids4stan <- metainfo%>%
  filter(isbest)%>%#these are thee final pairs of gene/cds/mass spec ids that we use for modeling MS
#  filter(sig_MS_change)%>%
  filter(n_stagemissing<=2)%>%#now filter them further - those with too many missing values will crash Rstan
  .$uprotein_id

gnm2uid<-metainfo%>%filter(uprotein_id %in% uids4stan)%>%distinct(gene_name,uprotein_id)%>%{safe_hashmap(.[[1]],.[[2]])}



get_count_mat_sigma <- function(uids4stan,countvoom,isribo=F){
	assaystring = if(isribo)'ribo' else 'total'

	mat <- rnamat <- countvoom$E%>%.[,str_subset(colnames(.),assaystring)]
	stopifnot(ncol(mat)==10)
	sigma <- countvoom$weights%>%{1/.}%>%
	  set_rownames(rownames(countvoom$E))%>%
	  set_colnames(colnames(countvoom$E))%>%
	  .[uids4stan,str_subset(colnames(.),assaystring),drop=F]%>%
	  set_rownames(uids4stan)
	list(mat[uids4stan,],sigma[uids4stan,])
}

scale_countmat <- function(uids4stan,rnamat){
	require(magrittr)
	require(proDD)
	rnamatrscl <- rnamat[uids4stan,]%>%
	  set_rownames(uids4stan)
	rnamatrscl <- rnamatrscl%>%{proDD::median_normalization(.)}
	rnamed <- rnamatrscl%>%apply(2,median,na.rm=T)%>%median(.,na.rm=T)
	rnamatrscl %<>% subtract(rnamed)%>%set_rownames(uids4stan)

}

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

matchedms_mat <- (readRDS)(here('data/matched_ms_mat.rds'))
mscountvoom <- (readRDS)(here('data/mscountvoom.rds'))

# c(mscountvoom$E) %<-% with(new.env(),{
# 	load('data/integrate_exprdata2.Rdata')
# 	list(mscountvoom$E)	
# })

#The below seems to be rna specific but isn't actually
#file.copy('../cortexomics/data/rnaproddparams.Rdata',here('data/proddparams.rda'))
if(!file.exists(here('data/proddparams.rda'))){
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
       tpvect,rna_sigma,file='data/proddparams.rda')
}else{
  library(zeallot)
  # c(matchedms_mat_rscl,rnamatrscl,rnamed,msmed,proddparam,tpvect) %<-% get()
  # file.remove(file='data/proddparams.Rdata')
  load(here('data/proddparams.rda'))
}

matchedms_mat_rscl<-matchedms_mat_rscl[uids4stan,]

get_dp_standata <- partial(get_dp_standata,
                           matchedms_mat_rscl=matchedms_mat_rscl[uids4stan,],
                           params=proddparam)



#add in prior data
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


#add in prior linear
get_dp_standata_withpriors_linear <- function(id,...){
  c(get_dp_standata(id,...),list(
    l_pihalf_priorsd= 0.001,
    l_pihalf_priormu= -4,
    l_st_priorsd = 3,
    l_ribo_priorsd = 3,
    l_ribo_priornu = 2.7
  )
  )
}


countmats <- list(
	rna = get_count_mat_sigma(uids4stan,mscountvoom,isribo=F)[[1]],
	ribo = get_count_mat_sigma(uids4stan,mscountvoom,isribo=T)[[1]]
)
sigmas <- list(
	rna = get_count_mat_sigma(uids4stan,mscountvoom,isribo=F)[[2]],
	ribo = get_count_mat_sigma(uids4stan,mscountvoom,isribo=T)[[2]]
)
sigmas%>%str
all((uids4stan %in% rownames(matchedms_mat_rscl)))


stopifnot('l_ribo_priornu' %in% names(get_dp_standata_withpriors(uids4stan[[1]],ribomatrscl=countmats[[1]],ribo_sigma=sigmas[[1]])))
stopifnot(get_dp_standata_withpriors(rownames(countmats[[1]])[1:3],ribomatrscl=countmats[[1]],ribo_sigma=sigmas[[1]])$lribo%>%dim%>%identical(c(3L,10L)))
stopifnot(sigmas[[2]][uids4stan[1:10],]%>%dim%>%identical(c(10L,10L)))


################################################################################
########Do optimization to get initial values
################################################################################

#Now optimize with many genes at once
#satb2like_uids <- as.matrix(edist)[gnm2uid[['Satb2']],]%>%sort%>%head(100)%>%names
#Random initializations don't work well, so we'lll 
satb2fits<-replicate(simplify=F,5,{
  pre = Sys.time()
  optres<-safely(rstan::optimizing)(dp_model,data=
  	get_dp_standata_withpriors(gnm2uid[['Satb2']],ribomatrscl=countmats[[1]],ribo_sigma=sigmas[[1]]),
  	algorithm='Newton',as_vector=FALSE,hessian=TRUE)
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
initvals <- bestoptes$par
for(p in names(initvals))initvals[[p]] <- initvals[[p]] - initvals[[p]]
initvals[['sigma2']] = initvals[['sigma2']]+0.2

counttype='ribo'
selgene=gnm2uid[['Satb2']]

fittype = 'ribo_linear'
fittype = 'ribo'
fittypes <- c('rna','ribo','rna_linear','ribo_linear') 
for(fittype in fittypes){
  message(fittype)
  optfitlistfile <- str_interp('data/${fittype}_geneoptfits.rds')
  counttype = fittype%>%str_extract('^[^_]+')
  paramtype = fittype%>%str_extract('(?<=_)[^_]+$')
  
  exprmat = countmats[[counttype]]
  sigma =  sigmas[[counttype]]
  if(!file.exists(optfitlistfile)){
  #
  require('R.utils')
  genefits <- list()
  #       
  genefits <- mclapply(mc.cores=20,uids4stan%>%setNames(.,.),(function(selgene){
    message(paste0('.',which(uids4stan==selgene),'..'))
      #
      selgenefits<-lapply(rep(selgene,5),safely(function(selgene){
        withTimeout(timeout=15,{
          # message(selgene)
          datafun = if('linear'%in%paramtype) get_dp_standata_withpriors_linear else get_dp_standata_withpriors
          optdat <- datafun(selgene,exprmat,sigma)
          if(sum(is.finite(optdat$lMS))<2) return(NULL)
          #apparently we need actual hessian - Netwon succeeds for individual genes where
          #L-BFGS does not. 
          #optimizing(dp_model,data=optdat,init=initvals,algorithm='Newton',as_vector=F,hessian=TRUE,verbose=TRUE) 
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

Sys.glob(str_interp('data/*_geneoptfits.rds'))
# file.remove(Sys.glob(str_interp('data/*_geneoptfits.rds')))
counttype='rna'
bestfits<-lapply(fittypes%>%setNames(.,.),function(counttype){
  optfitlistfile <- str_interp('data/${counttype}_geneoptfits.rds')
  genefits <- readRDS(optfitlistfile)

  genefitsres <- genefits%>%map(.%>%map('result')%>%keep(Negate(is.null)))
  genefitsres%>%map_dbl(length)%>%table

  genefitsres <- genefitsres%>%map(.%>%.[ map(.,'return_code')%>%map_lgl(~.==0)])
  genefitsres%>%map_dbl(length)%>%table

  #get the best fit object for each gene
  bestfits <- lapply(genefitsres,function(selgenefits){
	bestfit<-selgenefits%>%map_dbl('value')%>%which.max
    selgenefits[[bestfit]]
  })
  bestfits

})
bestfits%>%map(length)
bestfitinits<-bestfits%>%lapply(.%>%map('par'))

################################################################################
########Collect sampling on the cluster
################################################################################

{
  require(BiocParallel)
  #BiocManager::install('batchtools')
  # bpparam<- BatchtoolsParam(workers=20, cluster="slurm", resources=list(queue='all'))
  # sge_template <- '~/tmp.tmpl' 
  # bpparam$template%>%readLines
  bpparam<- BatchtoolsParam(
  	workers=20, cluster="slurm", resources=list(walltime = 4*(60)*(60),queue='all',ncpus=4),
  	saveregistry = TRUE,registryargs = batchtoolsRegistryargs(file.dir='/fast/work/groups/ag_ohler/dharnet_m/cortexomics_modelling/bplapplydebug')
)
  #bplapply(BPPARAM=param
  #library(GenomicRanges)
  #result <- bplapply(BPPARAM=bpparam, c('1:1','1:3'), function(x)GenomicRanges::GRanges(x) )
}


estimates <- bplapply(BPPARAM=bpparam,uids4stan%>%setNames(.,.),
# estimates <- lapply(uids4stan%>%setNames(.,.),
  arglist=c('countmats','sigmas','uids4stan','get_dp_standata_withpriors','bestfitinits','rna_sigma','get_dp_standata','rnamatrscl','dp_model','matchedms_mat_rscl','proddparam')%>%setNames(.,.)%>%lapply(get),
  function(uid,arglist){
  attach(arglist)
  require(tidyverse)
  require(rstan)
  library(splines)
  library(purrr)
  library(magrittr)
  library(tidyverse)
  library(here)

  outfiles = list()
  for(runtype in c('rna','ribo')){
  	countmat = countmats[[runtype]]
  	countsigma = sigmas[[runtype]]
 	
 	 rdsfile <- str_interp(here('data/${runtype}_fits/${uid}_dpfit.rds'))
  	rdsfile%>%dirname%>%dir.create(showWarnings=F)
  	if(!file.exists(rdsfile)){
  	#gnm2uid[['Satb2']]
  	modeldata <- get_dp_standata_withpriors(uid,ribomatrscl=countmat,countsigma)
  	init <- bestfitinits[[runtype]][[uid]]
  	attempts = 0
  	notdone=TRUE
  	while(notdone){

	    sampling<-safely(rstan::sampling)(dp_model,data=modeldata,chains=4,
	                            init=function(){init},control=list(adapt_delta=.98,max_treedepth=15),iter=2e3)
	    attempts = attempts + 1
	    
	    if(!is.null(sampling$result)){
	    	sampling = sampling$result;
	    	notdone=FALSE
	    }else{
	    	if(attempts>=1){
	    		sampling = sampling$error
	    		notdone = FALSE
	    	}
	    }

  	}
	saveRDS(sampling,rdsfile)
    }else{
      sampling<-readRDS(rdsfile)
    }
    outfiles = append(outfiles,rdsfile)
  }
	outfiles

})


Sys.glob(here('data/*_fits/*_dpfit.rds'))


################################################################################
########Plot Riboseq vs RNA difference.
################################################################################
##we now have the modes for each gene, let's combine these into a list of parameters for the joint model
get_comb_initvals <- function(bestfitinits){
  combinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
    bestfitinits%>%
      map(argind)%>%
      setNames(.,seq_along(.))%>%
      do.call(what=partial(abind::abind,along=1))
  })
  combinitvals	
}
combinitvals <- bestfitinits%>%map(~ .[uids4stan])%>%lapply(get_comb_initvals)

protmatchvect = matchedms_mat_rscl%>%colnames%>%str_replace('_\\d+$','')%>%as.factor%>%as.integer

matchedms_mat_rscl - combinitvals[['ribo']]$prot[,protmatchvect]
prodadf <- matchedms_mat_rscl%>%colnames%>%data.frame(sample=.)%>%separate(sample,c('time','assay','rep'))




library(proDA)
library(txtplot)

uids4errormeas = metainfo%>%filter(isbest)%>%.$uprotein_id%>%unique%>%intersect(uids4stan)

riboerror = rowSums((matchedms_mat_rscl[uids4errormeas,] - combinitvals[['ribo']]$prot[match(uids4errormeas,uids4stan),protmatchvect])^2,na.rm=TRUE)
rnaerror = rowSums((matchedms_mat_rscl[uids4errormeas,] - combinitvals[['rna']]$prot[match(uids4errormeas,uids4stan),protmatchvect])^2,na.rm=TRUE)
riboerror = rowSums((matchedms_mat_rscl[uids4errormeas,] - combinitvals[['ribo_linear']]$prot[match(uids4errormeas,uids4stan),protmatchvect])^2,na.rm=TRUE)
rnaerror = rowSums((matchedms_mat_rscl[uids4errormeas,] - combinitvals[['rna_linear']]$prot[match(uids4errormeas,uids4stan),protmatchvect])^2,na.rm=TRUE)


countmatchvect = countmats[[1]]%>%colnames%>%str_replace('_\\d+$','')%>%as.factor%>%as.integer
cmatmeans = lapply(countmats,function(cmat){
  lapply(unique(protmatchvect),function(i) rowMeans(cmat[,countmatchvect==i],na.rm=TRUE))%>%bind_cols
})


dteuproteinids<-metainfo%>%filter(dTE)%>%.$uprotein_id%>%unique

plot_riboerr = log2(riboerror/rnaerror)%>%enframe('uprotein_id','log2_ribo_err_v_rna_err')%>%
    mutate(dTE=uprotein_id%in%dteuproteinids)%>%
    filter(log2_ribo_err_v_rna_err%>%{between(.,quantile(.,0.02),quantile(.,1-0.02))})%>%
    qplot(data=.,geom='histogram',fill=dTE,x=log2_ribo_err_v_rna_err,alpha=I(1))+facet_grid(dTE~.,scales = 'free_y')+
    lims(x=c(-3,3))+
    ggtitle('Protein Variance Explained by Riboseq vs RNA\n(upper/lower 2% quantiles trimmed)')


plot_riboerr_lin = log2(riboerror_lin/rnaerror_lin)%>%enframe('uprotein_id','log2_ribo_err_v_rna_err')%>%
  mutate(dTE=uprotein_id%in%dteuproteinids)%>%
  filter(log2_ribo_err_v_rna_err%>%{between(.,quantile(.,0.02),quantile(.,1-0.02))})%>%
  qplot(data=.,geom='histogram',fill=dTE,x=log2_ribo_err_v_rna_err,alpha=I(1))+facet_grid(dTE~.,scales = 'free_y')+
  lims(x=c(-3,3))+
  ggtitle('Protein Variance Explained by Riboseq vs RNA\nlinear_model\n(upper/lower 2% quantiles trimmed)')

ggarrange(ncol=2,plotlist=list(plot_riboerr,plot_riboerr_lin))

log2(riboerror/riboerror_lin)%>%enframe('uprotein_id','log2_ribo_err_v_rna_err')%>%
  mutate(dTE=uprotein_id%in%dteuproteinids)%>%
  filter(log2_ribo_err_v_rna_err%>%{between(.,quantile(.,0.02),quantile(.,1-0.02))})%>%
  qplot(data=.,geom='histogram',fill=dTE,x=log2_ribo_err_v_rna_err,alpha=I(1))+facet_grid(dTE~.,scales = 'free_y')+
  ggtitle('Protein Variance Explained by Riboseq vs RNA\nlinear_model\n(upper/lower 2% quantiles trimmed)')



table(log2(riboerror/riboerror_lin) < 0)

table(riboerror>riboerror_lin)

table(log2(riboerror/riboerror_lin) < 0)

table(log2(riboerror/rnaerror) < 0)
table(log2(riboerror_lin/rnaerror_lin) < 0)

table(log2(riboerror/rnaerror) < 0)


save.image()

table(log2(riboerror/rnaerror) < 0)


trimquants <- function(x,q=0.02){
	q1 = quantile(x,q)
	q2 = quantile(x,1-q)
	x[ (x > q1)&(x<q2)]
}



fisher.test(table(riboerror < rnaerror,
	uids4stan %in% (metainfo%>%filter(dTE)%>%.$uprotein_id)
))



riboerror - rnaerror


library(txtplot)
#or with coorealtions
map_dbl(seq_along(uids4errormeas),function(i){
	ms = matchedms_mat_rscl[uids4errormeas,][i,]
	pred = combinitvals[['ribo']]$prot[match(uids4errormeas,uids4stan),protmatchvect][i,]
	cor(ms,pred,use='complete')
})%>%txtdensity

map_dbl(seq_along(uids4errormeas),function(i){
	ms = matchedms_mat_rscl[uids4errormeas,][i,]
	pred = combinitvals[['rna']]$prot[match(uids4errormeas,uids4stan),protmatchvect][i,]
	cor(ms,pred,use='complete')
})%>%txtdensity



signaldf <- countmats[1]%>%map_df(.%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(sample,signal,-uprotein_id)%>%separate(sample,into=c('time','assay','rep'))%>%group_by(time,assay,uprotein_id)%>%summarise(signal = mean(log2(2^signal))))%>%
	spread(assay,signal)

sample_info_df%>%safe_left_join(signaldf)


countmats%>%lapply(head)

countmats%>%head

#Or use proDA for this
fit <- proDA(matchedms_mat_rscl[uids4stan,], design = ~ time, 
             col_data = prodadf[,], reference_level = "E13")
prodameanests = (proDA::coefficients(fit)%>%{cbind(.[,1],.[,-1]+.[,1])} )
proDA::coefficients(fit)%>%dim
length(uids4stan)




devmats = list(
	ribodevmat = matchedms_mat_rscl[uids4stan,] - combinitvals[['ribo']]$prot[match(uids4errormeas,uids4stan),protmatchvect],
	rnadevmat = matchedms_mat_rscl[uids4stan,] - combinitvals[['rna']]$prot[match(uids4errormeas,uids4stan),protmatchvect],
	msmeandevmat = matchedms_mat_rscl[uids4stan,] - prodameanests[,protmatchvect],
	ribolindevmat = matchedms_mat_rscl[uids4stan,] - cmatmeans$ribo[,protmatchvect],
	rnalindevmat = matchedms_mat_rscl[uids4stan,] - cmatmeans$rna[,protmatchvect]
)

# #Or use proDA for this
# devfits <- devmats%>%mclapply(mc.cores=10,function(x){
# 	proDA(x, design = ~ time, col_data = prodadf[,], reference_level = "E13")
# })
# d
# timevars = devfits%>%map(.%>%proDA::coefficients(.)%>%as.data.frame%>%{cbind(.[,1],.[,-1]+.[,1])}%>%apply(1,var))


predmat = cmatmeans[[1]]
expr_pred_mat<-cbind(matchedms_mat_rscl[uids4stan,],predmat)%>%as.matrix
prodadfpred= rbind(prodadf,data.frame(time=unique(prodadf$time),assay='pred',rep=1))
compfit = proDA(expr_pred_mat,design = ~assay*time,col_data=prodadfpred,reference_level=c('E13'))


prodameanesttvars = prodameanests%>%as.data.frame%>%{cbind(.[,1],.[,-1]+.[,1])}%>%apply(1,var)

txtdensity(log10(trimquants(timevars[[1]] / prodameanesttvars,.00)))
txtdensity(log10(trimquants(timevars[[2]] / prodameanesttvars,.00)))
txtdensity(log10(trimquants(timevars[[3]] / prodameanesttvars,.00)))

riboexplainsmore = (timevars[[1]] / prodameanesttvars) > (timevars[[2]] / prodameanesttvars)


txtdensity( (timevars[[1]] / prodameanesttvars) - (timevars[[2]] / prodameanesttvars) )

riboexplainsmore%>%mean


devfits[[2]]%>%proDA::coefficients(.)%>%as.data.frame%>%{cbind(.[,1],.[,-1]+.[,1])}%>%apply(1,var)

proDA::coefficients()%>%as.data.frame%>%select(-Intercept)%>%apply(1,var)%>%sum


rind <- match(uids4errormeas,uids4stan)

sum(rowVars(prodameanests[rind,] - combinitvals[['rna']]$prot[rind,]))/1e6
sum(rowVars(prodameanests[rind,] - combinitvals[['ribo']]$prot[rind,]))/1e6

combinitvals[['ribo']]$prot[match(uids4errormeas,uids4stan),]%>%dim



map(proDA::coefficient_variance_matrices(fit),~diag(.)[2:5])


################################################################################
########Now plot
################################################################################




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
plotfile<- here(paste0('plots/','rnalKD_indiv_v_mcshane_pihalf','.pdf'))
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


tibble(uprotein_id=uids4stan,
opt_lpihalf = combinitvals$rna$l_pihalf)%>%
inner_join(metainfo%>%distinct(gene_name,uprotein_id))%>%inner_join(mcshanethalfs)%>%
{cor.test(.$opt_lpihalf,log2(.$half_life),method='spearman')}


tibble(uprotein_id=uids4stan,
opt_lpihalf = combinitvals$ribo$l_pihalf)%>%
inner_join(metainfo%>%distinct(gene_name,uprotein_id))%>%inner_join(mcshanethalfs)%>%
{cor.test(.$opt_lpihalf,log2(.$half_life),method='spearman')}

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


isource='ribo'
ipar='l_pihalf'

for(isource in names(fitlist)){
for(ipar in ipar ){

#now plot
plotfile<- here(paste0('plots/',source,'_',ipar,'_indiv_v_mcshane_pihalf','.pdf'))
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

chrmcounts = Sys.glob('../cortexomics/pipeline/star/data/*/*.bam') %>% str_subset(neg=T,'transcript')%>%setNames(.,basename(.))%>% map_dbl(.%>%{x=.;system(str_interp('samtools view -c ${.} chrM'),intern=TRUE)}%>%as.numeric)

chrmcounts%>%enframe%>%mutate(name=str_replace(name,'\\.bam',''))%>%filter


