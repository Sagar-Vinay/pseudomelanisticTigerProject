##Additional supplemental material for -
##High frequence of an otherwise rare phenotype in a small and isolated tiger population#
##Goes with Methods section - simulations

##Future projection simulations to estimate the time to fixation of mutant allele or wildtype allele under three population growth models and with or without migration ##
## Drift with constant size, logistic growth and 1 migrant per generation - time to fixation simulations ##
## Vinay Sagar; 09/11/2020 ##

##First, defining the variables##
Ni <- 10 ##starting population size##
p0 <- 0.5 ##starting allele frequency ##
ngen <- 10000 ##arbitrary number assigned for number of generations, the simulation will end when one allele gets fixed over the other ##
nrep <- 2000 ##Number of repeats ##

##I want to store the vector of time to fixation with each population model into a single matrix##
mat.tFix <- matrix(NA, nrow = nrep, ncol = 6) ##Creating the matrix##
colnames(mat.tFix) <- c("constN", "constN-M1", "LG-K35", "LG-K35-M1", "LG-K104", "LG-K104-M1")
##Also need to store the allele type fixed into a matrix##
mat.alleleFix <- matrix(NA, nrow = nrep, ncol = 6)
colnames(mat.alleleFix) <- c("constN", "constN-M1", "LG-K35", "LG-K35-M1", "LG-K104", "LG-K104-M1")

####*********************####
####** case 1 : ConstN **####
####*********************####

##create a vector to store the time to fixation (t.fix) in every repeatation##
t.fix <- numeric(length = nrep)
##create a vector to store the allele fixed in every repeatation##
allele.fix <- numeric(length = nrep)

for (j in 1:nrep) {
  
  ##building the population##
  geno.iter0 <- matrix(NA, nrow = Ni, ncol = 2) ##create a two column matrix to store the genotypes
  geno.iter0[,]=t(replicate(Ni, ifelse(runif(2,0,1)<=p0,"B","A"))) ##fill up the matrix
  
  start_af = length(which(geno.iter0[,]=="B"))/(Ni*2) ##defining allele frequency##
  save_af <- rep(-1,ngen) ##creating a vector to store allele frequency each generation
  save_af[1]=start_af ##Filling the first entry in the vector
  
  for (k in 2:ngen) {
    
    geno.iter <- matrix(NA, nrow = Ni, ncol = 2) ##create the genotype matrix##
    geno.iter[,]=t(replicate(Ni,ifelse(runif(2,0,1)<=start_af,"B","A"))) ##filling up the matrix
    
    save_af[k] <- length(which(geno.iter[,]=="B"))/(Ni*2) ##storing the allele frequency each generation
    start_af <- save_af[k] ##updating start_af for the next run of the loop
    
    if(save_af[k]==0||save_af[k]==1){
      break
    }
  }
  save_af <- unique(save_af)
  t.fix[j]=length(save_af)-1 ##saving the nth generation when the allele gets fixed
  allele.fix[j]=ifelse(save_af[length(save_af)-1]==0,0,1)
}

mat.tFix[,1] = t.fix
mat.alleleFix[,1] = allele.fix

####***************************************************####
####** Case 2 - constN with 1 migrant per generation **####
####***************************************************####

##create a vector to store the time to fixation (t.fix) in every repeatation##
t.fix <- numeric(length = nrep)
##create a vectore to store the allele fixed in every repeatation##
allele.fix <- numeric(length = nrep)

for (j in 1:nrep) {
  
  ##building the population##
  geno.iter0 <- matrix(NA, nrow = Ni, ncol = 2) ##create a two column matrix to store the genotypes
  geno.iter0[,]=t(replicate(Ni, ifelse(runif(2,0,1)<=p0,"B","A"))) ##fill up the matrix
  
  start_af = length(which(geno.iter0[,]=="B"))/(Ni*2) ##defining allele frequency##
  save_af <- rep(-1,ngen) ##creating a vector to store allele frequency each generation
  save_af[1]=start_af ##Filling the first entry in the vector
  
  for (k in 2:ngen) {
    
    geno.iter <- matrix(NA, nrow = Ni, ncol = 2) ##create the genotype matrix##
    geno.iter[,]=t(replicate(Ni,ifelse(runif(2,0,1)<=start_af,"B","A"))) ##filling up the matrix
    
    ##adding one additional individual in the denominator to have migration effect##
    save_af[k] <- length(which(geno.iter[,]=="B"))/((Ni+k-1)*2) ##storing the allele frequency each generation
    start_af <- save_af[k] ##updating start_af for the next run of the loop
    
    if(save_af[k]==0||save_af[k]==1){
      break
    }
  }
  save_af <- unique(save_af)
  t.fix[j]=length(save_af)-1
  allele.fix[j]=ifelse(save_af[length(save_af)-1]==0,0,1)
}

mat.tFix[,2] = t.fix
mat.alleleFix[,2] = allele.fix

####*******************************************************************####
####** Case 3 - Logistic growth with carrying capacity 35, isolation **####
####*******************************************************************####

K <- 35 ##carrying capacity, 35 or 104, i.e., 1.3 or 3.8 tigers per 100sqkm. Area of similipal = 2750sqkm##
r <- 0.15

##create a vector to store the time to fixation (t.fix) in every repeatation##
t.fix <- numeric(length = nrep)
##create a vectore to store the allele fixed in every repeatation##
allele.fix <- numeric(length = nrep)

for (j in 1:nrep) {
  
  ##building the population##
  geno.iter0 <- matrix(NA, nrow = Ni, ncol = 2) ##create a two column matrix to store the genotypes
  geno.iter0[,]=t(replicate(Ni, ifelse(runif(2,0,1)<=p0,"B","A"))) ##fill up the matrix
  
  start_af = length(which(geno.iter0[,]=="B"))/(Ni*2) ##defining allele frequency##
  save_af <- rep(-1,ngen) ##creating a vector to store allele frequency each generation
  save_af[1]=start_af ##Filling the first entry in the vector
  
  for (k in 2:ngen) {
    
    ##Here, pop size should change every generation according to logistic growth##
    
    C.int <- (K-Ni)/Ni ##integration constant in logistics growth##
    pop.size <- numeric(length = ngen-1) ##Defining a pop size change vector##
    pop.size[k-1]=round(K/(1+(C.int*(2.718^((-1*r)*(k-1))))),0) ##equation for population size change##
    
    geno.iter <- matrix(NA, nrow = pop.size[k-1], ncol = 2) ##create the genotype matrix##
    geno.iter[,]=t(replicate(pop.size[k-1],ifelse(runif(2,0,1)<=start_af,"B","A"))) ##filling up the matrix
    
    save_af[k] <- length(which(geno.iter[,]=="B"))/(pop.size[k-1]*2) ##storing the allele frequency each generation
    start_af <- save_af[k] ##updating start_af for the next run of the loop
    
    if(save_af[k]==0||save_af[k]==1){
      break
    }
  }
  save_af <- unique(save_af)
  t.fix[j]=length(save_af)-1
  allele.fix[j]=ifelse(save_af[length(save_af)-1]==0,0,1)
}

mat.tFix[,3] = t.fix
mat.alleleFix[,3] = allele.fix

####**********************************************************************************####
####** case 4 - Logistic growth with carrying capacity 35, 1 migrant per generation **####
####**********************************************************************************####

K <- 35 ##carrying capacity, 35 or 104, i.e., 1.3 or 3.8 tigers per 100sqkm. Area of similipal = 2750sqkm##
r <- 0.15

##create a vector to store the time to fixation (t.fix) in every repeatation##
t.fix <- numeric(length = nrep)
##create a vectore to store the allele fixed in every repeatation##
allele.fix <- numeric(length = nrep)

for (j in 1:nrep) {
  
  ##building the population##
  geno.iter0 <- matrix(NA, nrow = Ni, ncol = 2) ##create a two column matrix to store the genotypes
  geno.iter0[,]=t(replicate(Ni, ifelse(runif(2,0,1)<=p0,"B","A"))) ##fill up the matrix
  
  start_af = length(which(geno.iter0[,]=="B"))/(Ni*2) ##defining allele frequency##
  save_af <- rep(-1,ngen) ##creating a vector to store allele frequency each generation
  save_af[1]=start_af ##Filling the first entry in the vector
  
  for (k in 2:ngen) {
    
    ##Here, pop size should change every generation according to logistic growth##
    
    C.int <- (K-Ni)/Ni ##integration constant in logistics growth##
    pop.size <- numeric(length = ngen-1) ##Defining a pop size change vector##
    pop.size[k-1]=round(K/(1+(C.int*(2.718^((-1*r)*(k-1))))),0) ##equation for population size change##
    
    geno.iter <- matrix(NA, nrow = pop.size[k-1], ncol = 2) ##create the genotype matrix##
    geno.iter[,]=t(replicate(pop.size[k-1],ifelse(runif(2,0,1)<=start_af,"B","A"))) ##filling up the matrix
    
    ##Adding one individual here in the denominator##
    save_af[k] <- length(which(geno.iter[,]=="B"))/((pop.size[k-1]+1)*2) ##storing the allele frequency each generation
    start_af <- save_af[k] ##updating start_af for the next run of the loop
    
    if(save_af[k]==0||save_af[k]==1){
      break
    }
  }
  save_af <- unique(save_af)
  t.fix[j]=length(save_af)-1
  allele.fix[j]=ifelse(save_af[length(save_af)-1]==0,0,1)
}

mat.tFix[,4] = t.fix
mat.alleleFix[,4] = allele.fix

####********************************************************************####
####** case 5 - Logistic growth with carrying capacity 104, isolation **####
####********************************************************************####

K <- 104 ##carrying capacity, 35 or 104, i.e., 1.3 or 3.8 tigers per 100sqkm. Area of similipal = 2750sqkm##
r <- 0.15

##create a vector to store the time to fixation (t.fix) in every repeatation##
t.fix <- numeric(length = nrep)
##create a vectore to store the allele fixed in every repeatation##
allele.fix <- numeric(length = nrep)

for (j in 1:nrep) {
  
  ##building the population##
  geno.iter0 <- matrix(NA, nrow = Ni, ncol = 2) ##create a two column matrix to store the genotypes
  geno.iter0[,]=t(replicate(Ni, ifelse(runif(2,0,1)<=p0,"B","A"))) ##fill up the matrix
  
  start_af = length(which(geno.iter0[,]=="B"))/(Ni*2) ##defining allele frequency##
  save_af <- rep(-1,ngen) ##creating a vector to store allele frequency each generation
  save_af[1]=start_af ##Filling the first entry in the vector
  
  for (k in 2:ngen) {
    
    ##Here, pop size should change every generation according to logistic growth##
    
    C.int <- (K-Ni)/Ni ##integration constant in logistics growth##
    pop.size <- numeric(length = ngen-1) ##Defining a pop size change vector##
    pop.size[k-1]=round(K/(1+(C.int*(2.718^((-1*r)*(k-1))))),0) ##equation for population size change##
    
    geno.iter <- matrix(NA, nrow = pop.size[k-1], ncol = 2) ##create the genotype matrix##
    geno.iter[,]=t(replicate(pop.size[k-1],ifelse(runif(2,0,1)<=start_af,"B","A"))) ##filling up the matrix
    
    save_af[k] <- length(which(geno.iter[,]=="B"))/(pop.size[k-1]*2) ##storing the allele frequency each generation
    start_af <- save_af[k] ##updating start_af for the next run of the loop
    
    if(save_af[k]==0||save_af[k]==1){
      break
    }
  }
  save_af <- unique(save_af)
  t.fix[j]=length(save_af)-1
  allele.fix[j]=ifelse(save_af[length(save_af)-1]==0,0,1)
}

mat.tFix[,5] = t.fix
mat.alleleFix[,5] = allele.fix

####***********************************************************************************####
####** Case 6 - Logistic growth with carrying capacity 104, 1 migrant per generation **####
####***********************************************************************************####

K <- 104 ##carrying capacity, 35 or 104, i.e., 1.3 or 3.8 tigers per 100sqkm. Area of similipal = 2750sqkm##
r <- 0.15

##create a vector to store the time to fixation (t.fix) in every repeatation##
t.fix <- numeric(length = nrep)
##create a vectore to store the allele fixed in every repeatation##
allele.fix <- numeric(length = nrep)

for (j in 1:nrep) {
  
  ##building the population##
  geno.iter0 <- matrix(NA, nrow = Ni, ncol = 2) ##create a two column matrix to store the genotypes
  geno.iter0[,]=t(replicate(Ni, ifelse(runif(2,0,1)<=p0,"B","A"))) ##fill up the matrix
  
  start_af = length(which(geno.iter0[,]=="B"))/(Ni*2) ##defining allele frequency##
  save_af <- rep(-1,ngen) ##creating a vector to store allele frequency each generation
  save_af[1]=start_af ##Filling the first entry in the vector
  
  for (k in 2:ngen) {
    
    ##Here, pop size should change every generation according to logistic growth##
    
    C.int <- (K-Ni)/Ni ##integration constant in logistics growth##
    pop.size <- numeric(length = ngen-1) ##Defining a pop size change vector##
    pop.size[k-1]=round(K/(1+(C.int*(2.718^((-1*r)*(k-1))))),0) ##equation for population size change##
    
    geno.iter <- matrix(NA, nrow = pop.size[k-1], ncol = 2) ##create the genotype matrix##
    geno.iter[,]=t(replicate(pop.size[k-1],ifelse(runif(2,0,1)<=start_af,"B","A"))) ##filling up the matrix
    
    ##Adding one individual here in the denominator##
    save_af[k] <- length(which(geno.iter[,]=="B"))/((pop.size[k-1]+1)*2) ##storing the allele frequency each generation
    start_af <- save_af[k] ##updating start_af for the next run of the loop
    
    if(save_af[k]==0||save_af[k]==1){
      break
    }
  }
  save_af <- unique(save_af)
  t.fix[j]=length(save_af)-1
  allele.fix[j]=ifelse(save_af[length(save_af)-1]==0,0,1)
}

mat.tFix[,6] = t.fix
mat.alleleFix[,6] = allele.fix

write.csv(mat.tFix, file = "tFix-6popmodels-Ni10-pi0.5-2.csv")
write.csv(mat.alleleFix, file = "nBFix-6popmodels-Ni10-pi0.5-2.csv")
