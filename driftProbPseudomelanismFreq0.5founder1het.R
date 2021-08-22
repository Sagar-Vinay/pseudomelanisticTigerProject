##Additional supplemental material for -
##High frequence of an otherwise rare phenotype in a small and isolated tiger population#
##Goes with Methods section - Simulations ##
##Random genetic drift after bottleneck - discrete generation time forward simulations ##
##** Vinay Sagar, Aug 20, 2020 **##

##Defining variables##
##1. Starting population size, or bottleneck population size, Ni. Ni goes from 2 to 20 ##
##2. Time allowed for drift, ngen. ngen goes from 10 to 50 generations##
##3. The growth rate, r = 0.15
##4. The allele frequency at the time of bottleneck, p0. p0 has very low values; for one heterozygote in the founder poupulation situation, p0 is obsolete##
##5. The number of times each run has to be repeated, nrep. nrep = 1000
##6. The threshold probability to reach, af_thresh = 0.5

Ni <- c(2:20)
ngen <- c(10:50)
p0 <- 0.05
nrep <- 1000
af_thresh <- 0.5
K <- 104 ##carrying capacity, 35 or 104, i.e., 1.3 or 3.8 tigers per 100sqkm. Area of similipal = 2750sqkm##
r <- 0.15

##I need to repeat the simulation nrep times for each value of Ni and ngen##
##After each run, I need to save the allele frequency into a matrix and at the end of the loop calculate the probability of allele frequency attaining af_thresh value##

prob.mel.mat <- matrix(-1,nrow=length(Ni),ncol = length(ngen))
finalPopSize.mat <- matrix(-1,nrow=length(Ni),ncol = length(ngen))

for (i in seq_along(Ni)) {
  print(Ni[i])
  for (j in seq_along(ngen)) {
    
    C.int <- (K-Ni[i])/Ni[i] ##integration constant in logistics growth##
    pop.size <- numeric(length = ngen[j]) ##Defining a pop size change vector##
    pop.size=round(K/(1+(C.int*(2.718^((-1*r)*c(1:ngen[j]))))),0) ##equation for population size change##
    
    ##I want to save the population size in final generation for each value of Ni and ngen##
    finalPopSize.mat[i,j]= pop.size[length(pop.size)]
    
    ##Now that the population size per generation vector is ready, I need to create a population of allels each generation##
    ##two columns for two alleles, pop.size rows for that many individuals in the population##
    ##after each population is created, allele frequency is calculated and stored in a vector, and next generation is built using that value##
    
    ##The whole thing needs to be repeated nrep times and the allele frequency in the last generation needs to be stored in a matrix##
    
    af_ngen <- numeric(length = nrep) ##vector to store the allele frequency in the final generation for nrep repeatations
    
    for (m in 1:nrep) {
      
      ##first, I need to build the starting bottleneck population, with size Ni[i]
      geno.iter0 <- matrix(NA, nrow = Ni[i], ncol = 2) ##create a two column matrix to store the genotypes
      geno.iter0[,]=t(replicate(Ni[i], ifelse(runif(2,0,1)<=p0,"B","A"))) ##fill up the matrix
      
      #for one heterozygote in the founder population, above two lines of code should be replaced with the following two lines
      #geno.iter0 <- matrix("A", nrow = Ni[i], ncol = 2) ##create a two column matrix to store the genotypes and fill it with wildtype allele
      #geno.iter0[1,1]="B" ##Introduce a mutant allele in the popultion##
      
      start_af = length(which(geno.iter0[,]=="B"))/(Ni[i]*2) ##defining start_af for the begining of the loop
      save_af <- numeric(length = length(pop.size)) ##creating a vector to store allele frequency each generation
      
      for (k in seq_along(pop.size)) {
        geno.iter <- matrix(NA, nrow = pop.size[k], ncol = 2) ##create the genotype matrix##
        geno.iter[,]=t(replicate(pop.size[k],ifelse(runif(2,0,1)<=start_af,"B","A"))) ##filling up the matrix
        
        save_af[k] <- length(which(geno.iter[,]=="B"))/(pop.size[k]*2) ##storing the allele frequency each generation
        start_af <- save_af[k] ##updating start_af for the next run of the loop
      }
      af_ngen[m]=save_af[length(save_af)]
    }
    prob.mel<-length(which(af_ngen[]>=af_thresh))/nrep
    prob.mel.mat[i,j]=prob.mel
    
  }
}
##Need to assign names (as per values) to the rows and columns##
rnames <- paste0("Ni=",Ni) ##Row name vector
cnames <- paste0("ngen=",ngen) ##Column names vector
rownames(prob.mel.mat)=rnames ##assigning row names
colnames(prob.mel.mat)=cnames ##assigning column names

write.csv(prob.mel.mat,file = "mel-prob-matrix-startaf-0.05.csv")
