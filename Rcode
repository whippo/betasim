## R code- PART 1
## Simulation to assemble communities at random from
## a species pool of size gamma and the calculate beta diversity
## among communities
## first block of code plots the algebraic solution
## to the beta-partition- gamma diversity relationship under random sampling (i.e. Fig 2A),
## second block of code is a simulation approach to the same issue,
## while the third block plots the results of the simulation (i.e.- Fig S2A)
###############################################################
## Algebraic solution for the beta partition: (i.e. fig 2A) ###
################################################################
# vector of regional gamma diversities
gamma_vect <- seq(1,400, by = 5)
# ncom is the number of transects per site
ncom <- 10
#stems per transect:
num_stems  <- c( 5, 15, 30, 50,  75, 100 )
gamma<-NULL
gamma_obs<-NULL
beta_obs_matrix<-matrix(0,length(num_stems), length(gamma_vect))
probs <- NULL
for (j in 1:length(num_stems)) {
  beta_obs <- NULL
  for(i in 1:length(gamma_vect)){
    pool <- gamma_vect[i]
    quantile_vect <- seq(0,1,length.out = (pool +2))
    quantile_vect <- quantile_vect[2:(pool+1)]; #this will make evenly spaced draws from the lognormal species abundance distribution
    abundances <- qlnorm(quantile_vect)
    probs <- abundances/sum(abundances); #converts the abundance to probabilities
    prob_in_one_comm <- (1 - exp(num_stems[j]*log(1 - probs)))
    gamma_observed <- sum(1 - exp(num_stems[j]*ncom *log(1 - probs)
    ))
    in_how_many = NULL
    for(k in 1:ncom) {
      in_how_many = rbind((dbinom(k,size = ncom, prob =
                                    prob_in_one_comm)), in_how_many)
    }
    apply(in_how_many, MARGIN = 1, FUN=sum)/ncom->pi_vect
    beta_partition <- (1 - sum((seq(ncom,1,by=-
                                      1))*pi_vect/gamma_observed))
    beta_obs <- c(beta_obs, beta_partition)
    gamma_obs<-c(gamma_obs, gamma_observed)
    gamma<-c(gamma, gamma_vect[i])
  }
  beta_obs_matrix[j,] <- beta_obs
}
plot(gamma_vect, beta_obs_matrix[1,], type = "l", ylim = c(0,1), xlim=c(0,300),xlab="γ diversity", ylab="β diversity", cex.lab=1.5, lwd=2, main="Algebraic solution")
for (m in 2:length(num_stems)) {
  points(gamma_vect, beta_obs_matrix[m,], type = "l", lwd=2)
}
#########################
## start of simulation ##
#########################
## vector of gamma diversites:
gamma_vect=c(5,10,30,40,50,60,70,80,90, 100, 150, 200, 300, 400)
## number of plots per location/ simulation run (we use 10 plots here because it matches the Gentry observed datasets)
ncom<-10
##use uniform abundance distribution? If not, a lognormal abundance distribution is used
unif_abund=FALSE
## simulation replicates per gamma:
## set low here to run quickly- 1000 reps used in paper
reps=100
##number of stems per plot ('n' in the figures):
local_stem_density = c( 5, 15, 30, 50,  75, 100 )
gamma_obs<-NULL
means<-NULL
vars<-NULL
gamma<-NULL
stems<-NULL
mean_alpha<-NULL
beta_partition<-NULL
for(l in 1:length(local_stem_density)){
  for(i in 1:length(gamma_vect)){
    for(j in 1:reps){
      ## generate the species pool of size gamma, with abundances distribued as a lognormal: set to true: data:
        occur<-rlnorm(n=gamma_vect[i])
      #use uniform abundance distribution if logical is
      if(unif_abund){
        occur<-rep(1, length(occur))
      }
      ##build an empty species by site matrix to hold the
      matrix(data=rep(0, ncom*gamma_vect[i]), nrow=ncom,
             ncol=gamma_vect[i])->results
      dimnames(results)<-list(seq(1:ncom),seq(1:gamma_vect[i]))
      ## populate each community (row of matrix) by sampling individuals with replacement
      ## from the species pool:
      for(k in 1:ncom){
        local<-sample(1:gamma_vect[i],
                      size=local_stem_density[l], replace=TRUE, prob=occur)
        table(local)]->tab
        names(tab)<-tab
results[k, unlist(dimnames(results)[2]) %in%
      }
## remove any species from the species by site matrix that do not occur in any community:
  if(sum(apply(results, MARGIN=2, FUN=sum)==0)>0    ){
    results[,-(which(apply(results, MARGIN=2,
                           FUN=sum)==0))]->results
  }
## convert the abundance-based species by site matrix into a presence/ absence matrix:
  results->occur_mat
occur_mat[occur_mat>0]<-1
## calculate the alpha diversity of each community:
apply(occur_mat, MARGIN=1, FUN=sum)->alphas
## calculate the gamma diversity of the set of transects:
  gamma_observed<-ncol(results)
gamma_obs<-c(gamma_obs, gamma_observed)
gamma<-c(gamma, gamma_vect[i])
stems<-c(stems, local_stem_density[l])
mean_alpha<-c(mean_alpha, mean(alphas))
beta_partition<-c(beta_partition, 1-
                    mean(alphas)/ncol(results))
    }
  } }
##save results as a data frame:
data.frame(gamma_obs, gamma, stems, mean_alpha, beta_partition)->df
################
## end of simulation ##
################
#################################################
##plot results of simulation (e.g. Fig. S2A) ####
#################################################
plot(c(0,0), xlab="γ diversity", ylab="β diversity", cex.lab=1.5, xlim=c(0, 300), ylim=c(0, 1.0), type="n", main="Simulations")
for(l in 1:length(local_stem_density)){
  subset(df, df$stems==local_stem_density[l])->sub
  tapply(sub$gamma_obs, INDEX=list(sub$gamma), FUN=mean)->gamma_means
  data.frame(gamma_means)->gm
  gm$gamma<-row.names(gm)
  merge(sub, gm)->sub
  tapply(sub$beta_partition, INDEX=list(sub$gamma_means),
         FUN=mean)->means
  data.frame(as.numeric(names(means)), means)->summary
  names(summary)[1]<-"gamma_regional"
  lines(c(0,means)~c(0,gamma_regional), data=summary, lwd=2)
  xpos<-max(gamma_means)-30
  #if(xpos<80){xpos<-80}
  #if(xpos>700){xpos<-700}
  text(xpos, max(means)+.03, paste("n =",local_stem_density[l]),
       pos=4, cex=1.2)
}
# R Code Part 2
# Function to calculate beta deviations (a standard effect size)
# Calculates observed gamma, alpha and beta diversities, mean expected (null) beta, and
# the standard effect size (SES), which measures the deviation of observed beta from expected beta
# where, SES = (observed.beta-mean.expected.beta)/sd.expected.beta
# The function takes the data for a location that contains multiple plots (e.g., a single location from the Gentry dataset)
#  the data file should have one row for each individual, with two columns:
  #  'transect' column contains the sample/transect name where the individual is located
#  'spp' column contains the species identification of the individual
# Note that in this format, there will be multiple rows with identical information,
# since species can have multiple individuals present in a transect.
# This function is for use with abundance data  (not presence/absence data)!
  #see example data file:  "sample.Gentry.plot.data.for.Rcode.txt"
  #################################################################
####################################################
## START FUNCTION
#################################################################
####################################################
ses.beta.function=function(gdata, Nsim=999)  # 'gdata' is the
  data file, 'Nsim' is the number of randomizations for
calculating expected/null beta
{
  plot.gamma=length(unique(gdata$spp))   #calculate the total number of species at the site
  transect.spp=tapply(gdata$spp,gdata$transect,unique)  #generate species list for each transect
  obs.mean.alpha=mean(sapply(transect.spp,length))  #calculate average number of species per transect
  obs.beta=1-obs.mean.alpha/plot.gamma  #calculate observed beta partition
  rand.mean.alpha=vector(length=Nsim)  #create empty vector to be filled in with randomly generated alpha values
  for(j in 1:Nsim)   #start loop for simulations
  {
    samp=sample(gdata$transect,length(gdata$transect),replace=F)
    #swaps order of plotnames
    swap.data=data.frame("transect"=samp,"spp"=gdata$spp)
    #assigns random plotnames to individuals
    rand.transect.spp=tapply(swap.data$spp,swap.data$transect,unique
    )  #generate species list for each transect
    rand.mean.alpha[j]=mean(sapply(rand.transect.spp,length))
    #calculate average number of species per transect
  }  #end loop for simulations
  null.plot.beta=1-rand.mean.alpha/plot.gamma  #calculates the 1000 random beta values
  mean.null.beta=mean(null.plot.beta) #calculates the mean of the random beta values
  sd.null.beta=sd(null.plot.beta)  #calculates the sd of the random beta values
  ses.beta=(obs.beta-mean.null.beta)/sd.null.beta  #calculates the deviation of the observed from expected (random) beta
  return(c("gamma"=plot.gamma,"obs.mean.alpha"=obs.mean.alpha,"obs.
           beta"=obs.beta,"mean.null.beta"=mean.null.beta,"sd.null.beta"=sd
           .null.beta,"ses.beta"=ses.beta))
}
#################################################################
####################################################
## END FUNCTION
#################################################################
####################################################
