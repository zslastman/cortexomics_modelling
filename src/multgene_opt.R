if(!exists('get_dp_standata_withpriors')) source('src/R/Modeling/allgene_sampling_clust.R')
source('src/R/Functions/rstan_functions.R')

dp_model_ribnorm = rstan::stan_model('src/Stan/mod_proDD_ribnorm.stan')

fix_param(dp_model_ribnorm,)

################################################################################
########Attempt sampling with multiple genes
################################################################################
edist <- dist(matchedms_mat_rscl)
satb2like_uids <- as.matrix(edist)[gnm2uid[['Satb2']],]%>%sort%>%head(10)%>%names
flnalike_uids <- as.matrix(edist)[gnm2uid[['Flna']],]%>%sort%>%head(10)%>%names



#use the individual estimates to get joint estimates
mgids <- c(satb2like_uids,flnalike_uids)
mginits <- bestfitinits%>%map(~ .[mgids])%>%lapply(get_comb_initvals)


mgsampling<-rstan::optimizing(dp_model_ribnorm,data=get_dp_standata_withpriors(c(mgids),ribomatrscl=countmats$ribo,ribo_sigma=sigmas$ribo),
                              init=mginits$ribo,iter=10)

mgsampling<-rstan::sampling(dp_model_ribnorm,data=get_dp_standata_withpriors(c(mgids),ribomatrscl=countmats$ribo,ribo_sigma=sigmas$ribo),chains=4,
                            init=function(){mginits$ribo},control=list(adapt_delta=.98,max_treedepth=15),iter=2e3)
