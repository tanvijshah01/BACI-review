########################################################
# JAGS code and R script applied to the wren dataset
# Article: New measures for evaluation of environmental perturbations using BACI analyses
# Authors: Mathieu Chevalier, James C Russell and Jonas Knape
# Journal: Ecological Applications
########################################################

rm(list=ls());gc() # clean R environment

### loading packages
require(runjags)
require(boa)
require(mcmcplots)
require(R2jags)



# navigate two folders upgetwd
wren_data_filename <- 'data_wren_EcolApp.Rdata'
data_path <- file.path(getwd(), "data", "raw_data", "Chevalier", wren_data_filename)
load(data_path)
data<-data.wren

# top_dir <- normalizePath(file.path(getwd(), "..", ".."))
#
# load(file.path(top_dir, "data", "raw_data", "Chevalier", wren_data_filename))
#------------------ Start of modeling section

### JAGS model
cat('model{

	### Likelihood
	for(i in 1:Nts){
		for(j in 1:Nyears){
			Y[i,j] ~ dnegbin(p[i,j],r) # Negative Binomial distribution
			p[i,j] <- r/(r+lambda[i,j])
			log(lambda[i,j]) <- alpha[i] + alpha.T[stormTime[j],stormSite[i]] + beta1*ObsAge[i,j] + beta2*FirstSurv[i,j] + eps[j] # Process eaquation
		}
		alpha[i] ~ dnorm(0,tau.alpha) # site-specific effects
	}

	### Define PRIORS
	for(j in 1:Nyears){
		eps[j] ~ dnorm(0,tau.eps) # time-specific effects
	}
	tau.eps <- pow(sd.eps,-2)
	sd.eps ~ dt(0,1,1)T(0,10) # standard deviation of time-specific effects (Half-Cauchy prior)

	tau.alpha <- pow(sd.alpha,-2)
	sd.alpha ~ dt(0,1,1)T(0,10)  # standard deviation of site-specific effects (Half-Cauchy prior)

	# Priors for observer effects
	beta1 ~ dnorm(0,0.1)
	beta2 ~ dnorm(0,0.1)

	r ~ dunif(0,100) # overdispersion parameter

	# Priors for intercepts in each period and treatment
	alpha.T[1,1] ~ dnorm(0,0.1) # Control-Before
	alpha.T[2,1] ~ dnorm(0,0.1) # Impact-Before
	alpha.T[1,2] ~ dnorm(0,0.1) # Control-After
	alpha.T[2,2] ~ dnorm(0,0.1) # Impact-After

	### Assess model fit using a sum-of-squared-type discrepancy (posterior predictive checks)
	for(i in 1:Nts){
		for(j in 1:Nyears){
			predicted[i,j] <- ((r*(1-p[i,j]))/p[i,j]+0.01) # Predicted values
			res.1[i,j] <- (Y[i,j]-predicted[i,j])/sqrt((r*(1-p[i,j]))/pow(p[i,j],2)+0.01) # standardized residuals for the current dataset

			Y.new[i,j] ~ dnegbin(p[i,j],r) # Generate a new dataset at each MCMC iteration
			res.2[i,j] <- (Y.new[i,j]-predicted[i,j])/sqrt((r*(1-p[i,j]))/pow(p[i,j],2)+0.01) # standardized residuals for the new dataset
		}
	}
	fit <- sum(pow(res.1[,],2)) # Sum of squared standardized residuals for the current dataset
	fit.new <- sum(pow(res.2[,],2)) # Sum of squared standardized residuals for the new dataset
	test <- step(fit.new-fit)
	bpvalue <- mean(test) # Bayesian P-values

}',file="jags_code_EcolApp.bugs")

### Define JAGS settings
store=1000
nadap=1000
chains=3
nburn=2000
thin=5

### Variables to monitor
my_var=c("alpha.T","alpha","sd.alpha","beta1","beta2","r","eps","sd.eps","bpvalue","fit","fit.new")

### data to feed JAGS
Nts=data$Nts # Number of time series (i.e. squares)

Nyears=data$Nyears # Number of years (i.e. 17 years from 1996 to 2012)

Y=data$Y # Count matrix (nrow=Nts; ncol=Nyears)
ObsAge=data$ObsAge # ObsAge matrix (nrow=Nts; ncol=Nyears); Standardized to zero mean and unit variance; Missing values have been replaced by zeroes after the standardization
FirstSurv=data$FirstSurv # Binary variable indicating whether the square was surveyed for the first time by an observer; Missing values have been replaced by zeroes
stormTime=data$stormTime # Indicator variable indicating whether the year considered is in the Before (coded 1) or the After (coded 2) period
stormSite=data$stormSite # Indicator variable indicating whether the square considered is in the control (coded 1) or the impact treatment (coded 2)
data.jags=list("Nts","Nyears","Y","ObsAge","stormTime","FirstSurv","stormSite")

### Parallel computation of the model (takes around 15min on a Intel(R) core(TM) i7-6700CPU @? 3.4GHz computer)
start=Sys.time()
out.mod=jags.parallel(data=data.jags,parameters.to.save=my_var,
                      model.file="jags_code_EcolApp.bugs",
                      export_obj_names=c("chains","store","thin","nburn"),
                      n.chains=chains,n.burnin=nburn,n.iter=(store*thin)+nburn,n.thin=thin)
save(out.mod,file=paste("out.mod.",species,".Rdata",sep=""))
end=Sys.time()
end-start # time elapsed for computation

#------------------ End of modeling section

#------------------ Start of post-processing section

##########
### Model diagnostics
##########
species = 'wren'

# Check model convergence for all parameters -> Gelman and Rubin diagnostic + check of MCMC chains
# load(paste("out.mod.",species,".Rdata",sep=""))
load(paste("out.mod.",".Rdata",sep=""))
out.mod=as.mcmc(out.mod)
gelman=gelman.diag(out.mod,multivariate=FALSE,autoburnin=FALSE)
max(gelman$psrf[,1],na.rm=T)# Rhat must be less than 1.1 for all parameters

traplot(out.mod,"alpha.T")
caterplot(out.mod,"alpha",reorder=F)
traplot(out.mod,"sd.alpha")
traplot(out.mod,"beta1")
traplot(out.mod,"beta2")
traplot(out.mod,"r")
traplot(out.mod,"eps")
traplot(out.mod,"sd.eps")

# Check model fit
comb=combine.mcmc(out.mod) # combine MCMC chains
plot(as.vector(comb[,"fit"]),as.vector(comb[,"fit.new"]))
abline(0,1)
bpvalue=round(mean(comb[,"bpvalue"]),2)
legend("topright",paste("b-pvalue=",bpvalue,sep=""),bty="n")

##########
### Extract posterior distribution of intercepts and draw figures
##########

# Posterior distribution of the four intercepts
# Since the estimates are on the log-scale, we need to add random site and random year effects to match with the exepected mean of a log-normal distribution; mean=exp(mu + sigma^2 / 2))
bef.ctrl=exp(comb[,"alpha.T[1,1]"] + (comb[,"sd.alpha"]^2)/2 + (comb[,"sd.eps"]^2)/2)
bef.imp=exp(comb[,"alpha.T[1,2]"] + (comb[,"sd.alpha"]^2)/2 + (comb[,"sd.eps"]^2)/2)
aft.ctrl=exp(comb[,"alpha.T[2,1]"] + (comb[,"sd.alpha"]^2)/2 + (comb[,"sd.eps"]^2)/2)
aft.imp=exp(comb[,"alpha.T[2,2]"] + (comb[,"sd.alpha"]^2)/2 + (comb[,"sd.eps"]^2)/2)

# Median and 95% HPD intervals of the four intercepts
all.inter=cbind(bef.ctrl,bef.imp,aft.ctrl,aft.imp)
est.inter=apply(all.inter,2,median)
HPD.inter=apply(all.inter,2,function(x)boa.hpd(x,alpha=0.05))
low.bound=HPD.inter[1,]
high.bound=HPD.inter[2,]

# Graphical representation of model estimates with respect to data (Fig. 4a)
yr=1996:2012
ctrl.Y=apply(Y[which(stormSite==1),],2,function(x)mean(x,na.rm=T)) # average counts in control sites for each year
imp.Y=apply(Y[which(stormSite==2),],2,function(x)mean(x,na.rm=T)) # average counts in impact sites for each year

plot(NULL,ylim=c(0,4),xlim=c(yr[1],yr[length(yr)]),xlab="Time",ylab="Average population counts",cex.lab=1,xaxt="n",
     yaxt="n",main=species)
axis(side=1,cex.axis=1,at=yr)
axis(side=2,cex.axis=1,las=2)
abline(v=2004.5,lty=2)
legend("topright",col=c("#1b9e77","#d95f02"),lty=1,legend=c("Control","Impact"),pch=16,bty="n",cex=1,lwd=1)

len=c(9,8) # vector containing the number of years before (1996-2004 = 9 years) and after (2005-2012 = 8 years) the perturbation
pos.A=seq(yr[len[1]]+0.5,yr[length(yr)]+0.5,length.out=len[2])
pos.B=seq(yr[1]-0.5,yr[len[1]]+0.5,length.out=len[1])

lines(pos.B,rep(est.inter[1],len[1]),col="#1b9e77",lwd=3)
polygon(c(rev(pos.B),pos.B),c(rev(rep(low.bound[1],len[1])),rep(high.bound[1],len[1])),col=adjustcolor("#1b9e77",0.3),border=NA)
lines(pos.A,rep(est.inter[3],len[2]),col="#1b9e77",lwd=3)
polygon(c(rev(pos.A),pos.A),c(rev(rep(low.bound[3],len[2])),rep(high.bound[3],len[2])),col=adjustcolor("#1b9e77",0.3),border=NA)
points(ctrl.Y~yr,pch=16,col="#1b9e77",cex=1.5)

lines(pos.B,rep(est.inter[2],len[1]),col="#d95f02",lwd=3)
polygon(c(rev(pos.B),pos.B),c(rev(rep(low.bound[2],len[1])), rep(high.bound[2],len[1])),col=adjustcolor("#d95f02",0.3),border=NA)
lines(pos.A,rep(est.inter[4],len[2]),col="#d95f02",lwd=3)
polygon(c(rev(pos.A),pos.A),c(rev(rep(low.bound[4],len[2])), rep(high.bound[4],len[2])),col=adjustcolor("#d95f02",0.3),border=NA)
points(imp.Y~yr,pch=16,col="#d95f02",cex=1.5)

##########
### Extract posterior distribution of the three measures of impact and draw figures
##########

# Posterior distribution of the three measures of impact
CI.div=abs(aft.imp-aft.ctrl)-abs(bef.imp-bef.ctrl)
CI.cont=abs(aft.imp-bef.imp)-abs(aft.ctrl-bef.ctrl)
BACI=(aft.imp-bef.imp)-(aft.ctrl-bef.ctrl)

# Median and 95% HPD intervals of the three measures of impact
all.mes=cbind(CI.div,CI.cont,BACI)
est.mes=apply(all.mes,2,median)
HPD.mes=apply(all.mes,2,function(x)boa.hpd(x,alpha=0.05))

# Graphical representation of the three measures of impact (Fig. 4b)
low.bound=HPD.mes[1,]
high.bound=HPD.mes[2,]
plot(NULL,xlab="",xlim=c(0.5,3.5),ylim=c(-0.1,1.3),xaxt="n",yaxt="n",cex.lab=1,ylab="Change in population counts",main=species)
axis(side=1,cex.axis=1,at=1:3,label=c("CI-divergence","CI-contribution","BACI contrast"))
axis(side=2,cex.axis=1,las=2)
points(x=1:3,est.mes,pch=19,cex=1.5)
segments(1:3,low.bound,1:3,high.bound,lwd=1.5)
abline(h=0,lty=2)

#------------------ End of post-processing section
