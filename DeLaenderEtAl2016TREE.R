
print("The quartz function only works on mac")

#Community model function
#changes = Vector containing level(s) 'l' of environmental change. Any length.
#slope = Vector containing responses to l. Lenght equal to nr of species.
#pool = Nr of species initially present.
#mu = Vector with growth rates. Lenght equal to nr of species.
#alpha = Matrix containing species interactions. 
#effect = Vector with max effects of the species on the EF (fmax in eq IV). 
#.........Lenght equal to nr of species.

Model      <- function(changes, slope, pool, mu, alpha, effect, Xi=rep(1,pool))
{
  Results  <- NULL
  for (change in changes)
  {
    resp   <- 1 + slope*change                   #Bracketed part of eq III
    respEF <- 1 + slope*2*change                 #Bracketed part of eq IV
    respEF <- respEF - (respEF<0)*respEF         #Convertion of negative functioning to zero
    X      <- Xi                                 #Initialize vector with species densities
    timestep <- 0.01                             #Define timesteps
    for (t in c(1:10000))                        #Solve model with loop
    {
      dXdt <- X*(mu*resp+alpha%*%X)              #Matrix notation of eq I
      X    <- X + dXdt*timestep                  
    }
    EF <- sum(effect*respEF*X)                   #EqII
    Richness <- sum(X>0.01)                      #Only consider species with density above threshold
    Density <- sum(X)                            #Total density
    Results <- rbind(Results,                    #Keep results
                     c(change, Richness, EF, 
                       Density))
  }
  return(Results)
}

#Running the model for the 4 scenarios
#1.1. defining the 4 cases 
pool        <- 10                                #max nr of species
alpha       <- diag(pool)                        #construction of alpha matrix begins here
inter       <- -0.2                              #interspecific interaction strength
alpha       <- alpha+inter-(2+inter)*alpha       
alpha       <- data.matrix(alpha)                #construction of alpha matrix done
#define amounts of env change for 4 cases (always the same)
changesAll  <- list(caseA=seq(0,1,length.out=10), 
               caseB=seq(0,1,length.out=10),
               caseC=seq(0,1,length.out=10),
               caseD=seq(0,1,length.out=10))
effects     <- rep(10, pool)                     #always fmax of 10
#growth rates (always the same for 4 cases)
mus         <- seq(10,20, length.out = pool)
#responses for 4 cases (different)
slopes      <- list(caseA=(seq(1,0,length.out=pool)-1)/1, #negative responses; low mu, small response
                    caseB=(seq(0,1,length.out=pool)-1)/1, #negative responses; low mu, large response
                    caseC=(seq(1,2,length.out=pool)-1)/1, #positive responses; low mu, small response
                    caseD=(seq(2,1,length.out=pool)-1)/1) #positive responses; low mu, large response

#1.2. random direct experiment
#10 sp to 5 sp
TenSp        <- c(1:pool) #Make vector of 1:10
#create all combinations of n out of 10 with n=10 to 5
Combinations <- list(TenSp=c(1:pool), NineSp=combn(TenSp, 9),
                     EightSp=combn(TenSp, 8), SevenSp=combn(TenSp, 7),
                     SixSp=combn(TenSp, 6), FiveSp=combn(TenSp, 5))
AllResults   <- NULL #initiate object
#Loop over them and run model in absence of env change
for (Comb in c(1:length(Combinations)))
{
  Combination <- Combinations[[Comb]]
  if (is.null(dim(Combination))) #Calculate separately for initial richness = 10
  {
    init              <- rep(0, pool) 
    init[Combination] <- 1 #only the species that should be present are sown
    Results           <- Model(changes=0, 
                               slope=1, #random choice, doens't matter since changes=0
                               pool=pool, 
                               mu=mus, 
                               alpha=alpha,
                               effect=effects,
                               Xi=init)
    AllResults        <- rbind(AllResults,Results)
    
  } else #Now calculate for the other richnesses
  {
    for (column in c(1:ncol(Combination)))
    {
      init                       <- rep(0, pool) 
      init[Combination[,column]] <- 1
      Results                    <- Model(changes=0, 
                                          slope=1, #only the species that should be present are sown
                                          pool=pool, 
                                          mu=mus, 
                                          alpha=alpha,
                                          effect=effects,
                                          Xi=init)
      AllResults                 <- rbind(AllResults,Results)
    }
  }
}

AllResults <- as.data.frame(AllResults)
MinResults <- aggregate(V3~V2,min,data=AllResults) #minimal EF per richness level 
MaxResults <- aggregate(V3~V2,max,data=AllResults) #maximal EF per richness level

#1.3. Indirect manipulations
labels     <- c("tolerant=rare", "tolerant=abundant", 
                "tolerant=rare", "tolerant=abundant")

quartz("",8,4.2,type="pdf",file="Fig.pdf")
par(mfrow=c(2,4), mar=c(3, 3, 2, 1) + 0.1, mgp=c(2, 0.5, 0))
AllResults <- NULL
for (case in c(1,3,2,4))
{
  Results  <- Model(changes=changesAll[[case]], 
                    slope=slopes[[case]], 
                    pool=pool, 
                    mu=mus, 
                    alpha=alpha,
                    effect=effects)
  #Calculate per-capita functioning f
  Results   <- cbind(Results, Results[,3]/Results[,4])
  
  #Plot level of change vs. effect on richness
  plot(Results[,1],
       (Results[,2]-Results[1,2])/Results[1,2]*100,t="l",
       lwd=3,ylim=c(-50,50), main="Level-dependent effects",
       xlab=expression(paste("Level of change ", italic("l"))), yaxt="n",
       ylab="Effect (%)")
  axis(2, at=c(-40,0,40), labels=c(-40,0,40))
  
  #Add legend for the first case only
  if (case==1) {legend("topleft", 
                       c("Per-capita", "Density", "Richness"), 
                       pch=c(NA, NA, NA), lwd=c(1, 2, 3))}
  
  #Add level of change vs. density 
  points(Results[,1],
         (Results[,4]-Results[1,4])/Results[1,4]*100,t="l", lwd=2)
  
  #Add level of change vs. per-capita functioning
  points(Results[,1],
         (Results[,ncol(Results)]-Results[1,ncol(Results)])/Results[1,ncol(Results)]*100,
         t="l",lwd=1)
  
  #The B-EF plot
  plot(Results[,2], Results[,3],t="p",
       xlim=c(5,10),
       xlab="Richness", ylab="Function",
       pch=21, cex=1.5, col="black", 
       main="Resulting B-EF", xaxt="n",
       bg=rev(heat.colors(nrow(Results), alpha=1)))
  axis(1, at=c(1:10), labels=c(1:10))
  
  #Add the random direct B-EF plot
  polygon(c(MinResults[,1],rev(MinResults[,1])), c(MinResults[,2], rev(MaxResults[,2])),
          col="grey", border = NA)
  #Add a point with the initial B and EF
  points(Results[1,2], Results[1,3],t="p", pch=4,cex=1.5)
  AllResults <- rbind(AllResults, Results)
}
dev.off()




