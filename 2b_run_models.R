# ---
# title: "Run spatiotemporal models for fisheries bycatch prediction"
# author: "Brian Stock"
# date: "Sept 29, 2017"
# output: html_vignette
# vignette: >
#   %\VignetteEngine{knitr::rmarkdown}
#   \usepackage[utf8]{inputenc} 
# ---

# source("/home/brian/Dropbox/bycatch/manuscript/spatial-bycatch/2b_run_models.R")

# load the data from 2a_process_survey
load("/home/brian/Dropbox/bycatch/manuscript/spatial-bycatch/wcann_processed.RData")

library(randomForest)
library(ROCR)
library(DMwR)
library(caret)
library(mgcv)
library(INLA)
INLA:::inla.dynload.workaround() 

# Use caret package to do stratified k-fold cross validation
k = 10
ind <- list()
species.bin <- c("DBRK","PHLB","YEYE")
for(sp in 1:length(species.bin)){
  sp.col <- paste0(species.bin[sp],"_01")
  ind[[sp]] <- createFolds(dat[,sp.col], k=k)
}

calc_AUC <- function(pred, obs){
  obs <- as.numeric(as.character(obs))
  predict = prediction(as.vector(pred), obs)
  AUC <- round(unlist(slot(performance(predict,"auc"),"y.values")),3)
  return(AUC)
}
calc_RMSE <- function(pred, obs){
  RMSE <- round(sqrt(mean((pred-obs)^2)),3)
  return(RMSE)
}

fit_GLM <- function(dat, sp.ind, covar, modeltype, fit.id, test.id){
  btime <- Sys.time()
  dat$z <- dat[,sp.ind]
	formula.glm = as.formula(paste0("z ~ -1 + YEAR + ", paste(covar, collapse=" + ")))
  if(modeltype=="binomial") fit <- gam(formula.glm, family="binomial", data=dat[fit.id,])
	if(modeltype=="positive") fit <- gam(formula.glm, family=Gamma(link="log"), data=dat[fit.id,])
	
	# calculate and return performance
	obs <- dat[test.id, sp.ind] # observations at test locations
	pred <- predict(fit, newdata=dat[test.id,], type='response')
	etime <- Sys.time()
	rtime <- etime - btime
	if(modeltype=="binomial"){
	  AUC <- calc_AUC(pred, obs)
	  return(list(AUC, fit, pred, obs, rtime))
	}
	if(modeltype=="positive"){
	  RMSE <- calc_RMSE(pred, obs)
	  return(list(RMSE, fit, pred, obs, rtime))
	}
}

fit_GAM_CONSTANT <- function(dat, sp.ind, covar, modeltype, fit.id, test.id){
  btime <- Sys.time()
  dat$z <- dat[,sp.ind]
	formula.gam.const = as.formula(paste0("z ~ -1 + YEAR + s(LON,LAT,k=100) + ",
	                                      paste(covar, collapse=" + ")))
  if(modeltype=="binomial") fit <- gam(formula.gam.const, family="binomial", data=dat[fit.id,])
	if(modeltype=="positive") fit <- gam(formula.gam.const, family=Gamma(link="log"), data=dat[fit.id,])

	# calculate and return performance
	obs <- dat[test.id, sp.ind] # observations at test locations
	pred <- predict(fit, newdata=dat[test.id,], type='response')
	etime <- Sys.time()
	rtime <- etime - btime
	if(modeltype=="binomial"){
	  AUC <- calc_AUC(pred, obs)
	  return(list(AUC, fit, pred, obs, rtime))
	}
	if(modeltype=="positive"){
	  RMSE <- calc_RMSE(pred, obs)
	  return(list(RMSE, fit, pred, obs, rtime))
	}
}

fit_GAM_IID <- function(dat, sp.ind, covar, modeltype, fit.id, test.id){
  btime <- Sys.time()
  dat$z <- dat[,sp.ind]
	formula.gam.iid = as.formula(paste0("z ~ -1 + s(LON,LAT,k=100,by=YEAR) + ",  paste(covar, collapse=" + ")))
  if(modeltype=="binomial") fit <- gam(formula.gam.iid, family="binomial", data=dat[fit.id,])
	if(modeltype=="positive") fit <- gam(formula.gam.iid, family=Gamma(link="log"), data=dat[fit.id,])
	
	# calculate and return performance
	obs <- dat[test.id, sp.ind] # observations at test locations
	pred <- predict(fit, newdata=dat[test.id,], type='response')
	etime <- Sys.time()
	rtime <- etime - btime
	if(modeltype=="binomial"){
	  AUC <- calc_AUC(pred, obs)
	  return(list(AUC, fit, pred, obs, rtime))
	}
	if(modeltype=="positive"){
	  RMSE <- calc_RMSE(pred, obs)
	  return(list(RMSE, fit, pred, obs, rtime))
	}
}

fit_GMRF <- function(dat, sp.ind, covar, modeltype, modeltype.GMRF, fit.id, test.id){
  btime <- Sys.time()
  # response needs to be numeric, not factor
  dat[,sp.ind] <- as.numeric(as.character(dat[,sp.ind]))
  # inRCA needs to be numeric, not factor
  dat$inRCA <- as.numeric(as.character(dat$inRCA))
  
  # turn 2012 to NA so YEAR=2012 will be intercept, other years will be offsets
  if(modeltype.GMRF=="CONSTANT"){
    dat$YEAR <- as.numeric(as.character(dat$YEAR))
    dat$YEAR[which(dat$YEAR==2012)] <- NA
    dat$YEAR <- as.factor(dat$YEAR)
  }
  # record keeping
  n.years <- length(levels(dat$YEAR))
  yr.labs <- levels(dat$YEAR)
  levels(dat$YEAR) <- 1:n.years
  n.sites <- dim(dat)[1]
  
  # CONSTANT includes offset terms for each year
  # EXCHANGEABLE includes separate GMRF for each year, so don't want YEAR terms
  if(modeltype.GMRF=="CONSTANT") covar = c(covar,"YEAR")
  n.covar <- length(covar)
  
  # difference between binomial and positive delta model components
  if(modeltype=="binomial") family.inla = "binomial" 
  if(modeltype=="positive") family.inla = "gamma" 
  n.fit <- length(fit.id)
  n.test <- length(test.id)
  
  # Set up mesh
  coords.fit = cbind(dat$LON[fit.id], dat$LAT[fit.id])
  coords.test = cbind(dat$LON[test.id], dat$LAT[test.id])
  bnd = inla.nonconvex.hull(coords.fit, convex=-0.05, concave=-0.2)
  mesh1 = inla.mesh.2d(loc=coords.fit, boundary=bnd, offset=c(-0.01,-0.03), 
                       cutoff=0.3, max.edge=c(.6,1.3))
  spde = inla.spde2.matern(mesh1, alpha=2) # alpha = matern parameter
  
  # Make index for spatial field
  #   CONSTANT: fit one spatial field
	#   EXCHANGEABLE: fit different spatial fields for each year
  if(modeltype.GMRF=="CONSTANT") iset <- inla.spde.make.index("i", n.spde=mesh1$n)
	if(modeltype.GMRF=="EXCHANGEABLE") iset <- inla.spde.make.index('i', n.spde=spde$n.spde,
	                                                                n.group=n.years)

  # Make A matrix (projects from mesh nodes to data locations)
  if(modeltype.GMRF=="CONSTANT"){
	  A <- inla.spde.make.A(mesh=mesh1, loc=coords.fit)
	  A.test <- inla.spde.make.A(mesh=mesh1, loc=coords.test)
  }
  if(modeltype.GMRF=="EXCHANGEABLE"){
	  A <- inla.spde.make.A(mesh=mesh1, loc=coords.fit, 
	                        group=as.numeric(dat$YEAR)[fit.id])
	  A.test <- inla.spde.make.A(mesh=mesh1, loc=coords.test, 
	                             group=as.numeric(dat$YEAR)[test.id])
  }
  A.list = list(); A.list[[1]] = A;
  for(i in 1:n.covar) A.list[[i+1]] <- 1;
  A.list.test = list(); A.list.test[[1]] = A.test;
  for (i in 1:n.covar) A.list.test[[i+1]] <- 1;

  # Make list of covariates including 'iset', the GMRF
  effect.list = list(i = iset, 
					sst = dat[fit.id,"sst"], 
					sst2 = dat[fit.id,"sst2"], 
					logDEPTH = dat[fit.id,"logDEPTH"], 
					logDEPTH2 = dat[fit.id,"logDEPTH2"], 
					inRCA = dat[fit.id,"inRCA"], 
					DAY = dat[fit.id,"DAY"])
  effect.list.test = list(i = iset, 
					sst = dat[test.id,"sst"], 
					sst2 = dat[test.id,"sst2"], 
					logDEPTH = dat[test.id,"logDEPTH"], 
					logDEPTH2 = dat[test.id,"logDEPTH2"], 
					inRCA = dat[test.id,"inRCA"], 
					DAY = dat[test.id,"DAY"])
  if(modeltype.GMRF=="CONSTANT"){ # CONSTANT model has YEAR terms in it
	  effect.list$YEAR <- dat[fit.id,"YEAR"]
	  effect.list.test$YEAR <- dat[test.id,"YEAR"]
  }
  
  # 'stack' data together
  # sdat.fit: data for GMRF model to fit
  sdat.fit <- inla.stack(tag='sdat.fit', data=list(z=dat[fit.id, sp.ind]), 
                         A=A.list, effects=effect.list)
  # sdat.test: data for GMRF model to predict
  # set z=NA to tell R-INLA to predict for these locations
  sdat.test <- inla.stack(tag='sdat.test', data=list(z=rep(NA, n.test)), 
                          A=A.list.test, effects=effect.list.test)
  sdat.full <- inla.stack(sdat.fit, sdat.test)

  # Make INLA formula
  if(modeltype.GMRF=="CONSTANT"){
    formula.inla = as.formula(paste0("z ~ -1 + ", 
            paste(covar, collapse="+"),
            "+ f(i, model=spde)"))} 
  if(modeltype.GMRF=="EXCHANGEABLE"){
    formula.inla = as.formula(paste0("z ~ -1 + ", 
            paste(covar, collapse="+"), 
            "+ f(i, model=spde, group=i.group, control.group=list(model='exchangeable'))"))}
  
  # Call INLA
  # quick run to find posterior mode, using gaussian approximation and 
  #   empirical Bayes integration strategy over the hyperparameters
  start.inla <- inla(formula.inla, num.threads=12, family = family.inla,
                     data = inla.stack.data(sdat.full), 
                     control.predictor = list(link=1, compute=FALSE, A=inla.stack.A(sdat.full)),
                     verbose = TRUE, debug=TRUE, keep=FALSE,
                     control.inla = list(strategy="gaussian", int.strategy="eb"),
                     control.compute = list(dic=TRUE, cpo=TRUE), 
                     control.fixed = list(expand.factor.strategy='inla', correlation.matrix=TRUE), 
                     control.results=list(return.marginals.random=FALSE,return.marginals.predictor=FALSE))
  # longer run using more accurate approximation, uses posterior mode found in previous step
  out.inla <- inla(formula.inla, num.threads=12, family = family.inla, 
                   data=inla.stack.data(sdat.full), 
                   control.predictor=list(link=1, compute=TRUE, A=inla.stack.A(sdat.full)), 
                   verbose = TRUE, debug=TRUE, keep=FALSE, 
                   control.compute = list(dic=TRUE,cpo=TRUE), 
                   control.fixed = list(expand.factor.strategy='inla',correlation.matrix=TRUE), 
                   control.mode = list(theta=start.inla$mode$theta, restart=FALSE), 
                   control.results=list(return.marginals.random=FALSE,return.marginals.predictor=FALSE))
  
  etime <- Sys.time()
	rtime <- etime - btime
  
  # Calculate and return performance metrics on test data (binomial)
  if(modeltype=="binomial"){
  	# Get predicted and observed
  	ind.pred <- inla.stack.index(sdat.full,'sdat.test')$data
  	pred <- out.inla$summary.fitted.values[ind.pred,"mean"]
  	obs <- dat[test.id, sp.ind]
  	AUC <- calc_AUC(pred, obs)
  	# Return AUC and model objects
  	fit <- list("AUC"=AUC,"out.inla"=out.inla,"pred"=pred,"obs"=obs,"rtime"=rtime,"mesh1"=mesh1,
  	            "iset"=iset,"sdat.full"=sdat.full,"test.id"=test.id,"fit.id"=fit.id,
  	            "n.test"=n.test,"n.fit"=n.fit)
  }
  
  # Calculate and return performance metrics (positive)
  if(modeltype=="positive"){
  	# Get predicted and observed
  	ind.pred <- inla.stack.index(sdat.full,'sdat.test')$data
  	pred <- out.inla$summary.fitted.values[ind.pred,"mean"]
  	obs <- dat[test.id, sp.ind]
  	RMSE <- calc_RMSE(pred, obs)
  	# Return RMSE and model objects
  	fit <- list("RMSE"=RMSE,"out.inla"=out.inla,"pred"=pred,"obs"=obs,"rtime"=rtime,
  	            "mesh1"=mesh1,"iset"=iset,"sdat.full"=sdat.full,"test.id"=test.id,
  	            "fit.id"=fit.id,"n.test"=n.test,"n.fit"=n.fit)
  }
  return(fit)
}

# RF BASE is the default randomForest function
fit_RF_BASE <- function(dat, sp.ind, covar, modeltype, fit.id, test.id){
  btime <- Sys.time()
  # keep forest for prediction at test locations
  fit <- randomForest(x=dat[fit.id, covar], y=dat[fit.id, sp.ind], 
                      xtest=dat[test.id, covar], ytest=dat[test.id, sp.ind],
                      mtry=3, ntree=1000, importance=FALSE, do.trace=250, keep.forest=TRUE) 
	
	# calculate and return performance
  etime <- Sys.time()
  rtime <- etime - btime
	obs <- dat[test.id, sp.ind] # observations at test locations
	if(modeltype=="binomial"){
	  pred <- predict(fit, newdata=dat[test.id,], type='prob')[,2]
	  AUC <- calc_AUC(pred, obs)
	  return(list(AUC, fit, pred, obs, rtime))
	}
	if(modeltype=="positive"){
	  pred <- predict(fit, newdata=dat[test.id,], type='response')
	  RMSE <- calc_RMSE(pred, obs)
	  return(list(RMSE, fit, pred, obs, rtime))
	}
}

# RF DOWN downsamples the majority class (binomial component only)
# downsample = if classes are imbalanced, train RF using equal #s of 0s and 1s
fit_RF_DOWN <- function(dat, sp.ind, covar, fit.id, test.id){
  btime <- Sys.time()
	# nmin <- sum(dat[fit.id, sp.ind])-1 # number of minority class (assume 1s)
	nmin <- table(dat[fit.id,sp.ind])[2]; names(nmin) <- NULL;
	prop1 <- nmin/length(fit.id)
	prop0 <- 1-prop1
	if(prop0 < prop1) nmin <- length(fit.id) - nmin # if 0s are minority, use # of 0s
  fit <- randomForest(sampsize=rep(round(nmin/6),2), x=dat[fit.id, covar], 
                      y=dat[fit.id, sp.ind], mtry=3, ntree=1000, 
                      importance=FALSE, do.trace=250, keep.forest=TRUE) 
  
	# calculate and return performance
  etime <- Sys.time()
  rtime <- etime - btime
	obs <- dat[test.id, sp.ind] # observations at test locations
  pred <- predict(fit, newdata=dat[test.id,], type='prob')[,2]
  AUC <- calc_AUC(pred, obs)
  return(list(AUC, fit, pred, obs, rtime))
}

# RF SMOTE is also designed to improve RF for imbalanced data (binomial component only)
# SMOTE = Synthetic Minority Over-sampling Technique
#   combines downsampling of majority class with oversampling of minority class
#   creates synthetic minority class instances by generating random linear combinations
fit_RF_SMOTE <- function(dat, sp.ind, covar, fit.id, test.id){
  btime <- Sys.time()
  prop <- table(dat[fit.id, sp.ind])[2] / length(fit.id) # get percent minority class (1s)
  p.over <- round(50/prop) # percent to oversample to get to 50%
  p.under <- round(100/(1-prop)) # percent to undersample to get to 50%
  X <- cbind(dat[fit.id, covar], dat[fit.id, sp.ind])
  names(X) <- c(covar, "z")
  formula.rf <- formula(paste0("z ~ ", paste0(covar, collapse=" + ")))
  X.SMOTE <- SMOTE(formula.rf, data=X, k=5, perc.over=p.over, perc.under=p.under)
  # table(X.SMOTE$z) # check now we roughly have class balance
  fit <- randomForest(x=X.SMOTE[,covar], y=X.SMOTE[,"z"], mtry=3, ntree=1000, importance=FALSE, do.trace=250, keep.forest=TRUE)

  # calculate and return performance
  etime <- Sys.time()
  rtime <- etime - btime
	obs <- dat[test.id, sp.ind] # observations at test locations
  # pred <- predict(fit, newdata=dat[test.id, covar], type='prob', predict.all=TRUE)
  pred <- predict(fit, newdata=dat[test.id, covar], type='prob')[,2]
  AUC <- calc_AUC(pred, obs)
  return(list(AUC, fit, pred, obs, rtime))
}

# Set up binomial model storage
species.bin <- c("DBRK","PHLB","YEYE")
n.species.bin <- length(species.bin)
models.bin <- c("GLM","GAM CONSTANT","GAM IID","GMRF CONSTANT","GMRF EXCHANGEABLE","RF BASE","RF DOWN","RF SMOTE")
n.models.bin <- length(models.bin)
AUC <- array(NA,dim=c(n.species.bin, n.models.bin, k))
fits.bin <- vector("list", n.species.bin)

# Set up positive model storage
species.pos <- c("DBRK","PHLB")
n.species.pos <- length(species.pos)
models.pos <- c("GLM","GAM CONSTANT","GAM IID","GMRF CONSTANT","GMRF EXCHANGEABLE","RF BASE")
n.models.pos <- length(models.pos)
RMSE <- array(NA,dim=c(n.species.pos, n.models.pos, k))
fits.pos <- vector("list", n.species.pos)

# Use same environmental covariates for all models
covar <- c("logDEPTH", "logDEPTH2", "sst", "sst2", "inRCA", "DAY")
# Random forest includes lat/lon as covariates instead of spatial structures
rf.covar <- c(covar, "YEAR", "LAT", "LON")

# Run binomial models
for(sp in 1:n.species.bin){ # for each species
  modeltype <- "binomial"
  sp.lab <- species.bin[sp]
  sp.col <- paste0(sp.lab,"_01")
  sp.ind <- match(sp.col, names(dat))
  fits.bin[[sp]] <- vector("list", n.models.bin)
  
  for(f in 1:k){ # for each fold
    test.id <- ind[[sp]][[f]] # get test rows for this species and fold (10% of data)
    fit.id <- dat[-test.id,"id"] # get rows to fit models (90% of data)
    
    # Fit GLM
    fits.bin[[sp]][[1]] <- fit_GLM(dat, sp.ind, covar, modeltype, fit.id, test.id)
    # Fit GAM CONSTANT
    fits.bin[[sp]][[2]] <- fit_GAM_CONSTANT(dat, sp.ind, covar, modeltype, fit.id, test.id)
    # Fit GAM IID
    fits.bin[[sp]][[3]] <- fit_GAM_IID(dat, sp.ind, covar, modeltype, fit.id, test.id)
    # Fit GMRF CONSTANT
    fits.bin[[sp]][[4]] <- fit_GMRF(dat, sp.ind, covar, modeltype, 
                                    modeltype.GMRF="CONSTANT", fit.id, test.id)
    # Fit GMRF EXCHANGEABLE
    fits.bin[[sp]][[5]] <- fit_GMRF(dat, sp.ind, covar, modeltype, 
                                    modeltype.GMRF="EXCHANGEABLE", fit.id, test.id)
    # Fit RF BASE
    fits.bin[[sp]][[6]] <- fit_RF_BASE(dat, sp.ind, rf.covar, modeltype, fit.id, test.id)
    # Fit RF DOWN
    fits.bin[[sp]][[7]] <- fit_RF_DOWN(dat, sp.ind, rf.covar, fit.id, test.id)    
    # Fit RF SMOTE
    fits.bin[[sp]][[8]] <- fit_RF_SMOTE(dat, sp.ind, rf.covar, fit.id, test.id)
  }
}

# Run positive models
for(sp in 1:n.species.pos){ # for each species
  modeltype <- "positive"  
  sp.lab <- species.pos[sp]
  sp.col <- sp.lab
  sp.ind <- match(sp.col, names(dat))
  fits.pos[[sp]] <- vector("list", n.models.pos)
  
  for(f in 1:k){ # for each fold
    test.id <- ind[[sp]][[f]] # get test rows for this species and fold (10% of data)
    fit.id <- dat[-test.id,"id"] # get rows to fit models (90% of data)
    
    # Only want to fit non-zero points in the positive model
    pos.fit <- which(dat[fit.id, sp.ind] > 0)
  	fit.id <- fit.id[pos.fit]
  	pos.test <- which(dat[test.id, sp.ind] > 0)
  	test.id <- test.id[pos.test]    
    
    # Fit GLM
    fits.pos[[sp]][[1]] <- fit_GLM(dat, sp.ind, covar, modeltype, fit.id, test.id)
    # Fit GAM CONSTANT
    fits.pos[[sp]][[2]] <- fit_GAM_CONSTANT(dat, sp.ind, covar, modeltype, fit.id, test.id)
    # Fit GAM IID
    fits.pos[[sp]][[3]] <- fit_GAM_IID(dat, sp.ind, covar, modeltype, fit.id, test.id)
    # Fit GMRF CONSTANT
    fits.pos[[sp]][[4]] <- fit_GMRF(dat, sp.ind, covar, modeltype, 
                                    modeltype.GMRF="CONSTANT", fit.id, test.id)
    # Fit GMRF EXCHANGEABLE
    fits.pos[[sp]][[5]] <- fit_GMRF(dat, sp.ind, covar, modeltype, 
                                    modeltype.GMRF="EXCHANGEABLE", fit.id, test.id)
    # Fit RF BASE
    fits.pos[[sp]][[6]] <- fit_RF_BASE(dat, sp.ind, rf.covar, modeltype, fit.id, test.id)
  }
}

save.image("/home/brian/Dropbox/bycatch/manuscript/spatial-bycatch/wcann_models_finished.RData")
