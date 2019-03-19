#1while

N = 50 # number of diploid individuals
N.chrom = 2*N # number of chromosomes
#p = .5; q = 1-p
#N.gen = 100 # number of generations
N.sim = 100 # number of simulations
w<-c(1,1,1) #a vector of three fitness values for genotype AA, Aa, and aa

w<-(w/max(w))

probabilities <- c()

for (freq in 1:99) {
  p <- freq/100
  q <- 1-p
  X = array(0, dim=c(1,N.sim)) #each column is a simulation, each row is a generation
  X[1,] = rep(N.chrom*p,N.sim) # initialize number of A1 alleles in first generation
  numfix <- 0
  for(j in 1:N.sim){
    i <- 2
    fixed <- FALSE
    lost <- FALSE
    while (!(fixed || lost)) {
      X <- rbind(X, i)
      genotypes<-matrix(rbinom(N.chrom,1,prob=X[i-1,j]/N.chrom),N,2) #generate an N x 2 matrix of genotypes, where the probability of drawing the A1 allele is given by its frequency in the previous generation
      fitnessv<-w[3:1][rowSums(genotypes)+1] #produce a probability vector of weights that each individual will contribute to the next generation, based on its genotype. Note that w needs to be inverted to keep proper track because the first value in w is for the AA genotype, which is given by allele combinations 1 and 1. 
      parents<-c(sample(N,prob=fitnessv,replace=T),sample(N,prob=fitnessv,replace=T)) #generate a vector of all parents who will contribute to the next generation (two for each offspring). Some parents contribute more than once (especially if they have selectively advantageous genotypes), and some are not represented at all. 
      nextgen<-sapply(parents,function(x) sample(genotypes[x,],1)) #sample one of the alleles from each parent, each time it contributed
      X[i,j] = sum(nextgen)
      if (all(nextgen == 1)) {
        fixed <- TRUE
      }
      if (all(nextgen == 0)) {
        lost <- TRUE
      }
      i <- i+1
    }
    if (fixed) {
      numfix <- numfix + 1
    }
  }
  probabilities = c(probabilities, numfix/N.sim)
}

probabilities

plot(probabilities, xlab = 'Initial allele frequency', ylab = 'Probability of fixation')



#1for

N = 50 # number of diploid individuals
N.chrom = 2*N # number of chromosomes
#p = .5; q = 1-p
N.gen = 100 # number of generations
N.sim = 100 # number of simulations
w<-c(1,1,1) #a vector of three fitness values for genotype AA, Aa, and aa
w<-(w/max(w))

probabilities = c()

for (freq in 1:99) {
  p <- freq/100
  q <- 1-p
  X = array(0, dim=c(N.gen,N.sim)) #each column is a simulation, each row is a generation
  X[1,] = rep(N.chrom*p,N.sim) # initialize number of A1 alleles in first generation
  numfix <- 0
  for(j in 1:N.sim){
    for(i in 2:N.gen){
      genotypes<-matrix(rbinom(N.chrom,1,prob=X[i-1,j]/N.chrom),N,2) #generate an N x 2 matrix of genotypes, where the probability of drawing the A1 allele is given by its frequency in the previous generation
      fitnessv<-w[3:1][rowSums(genotypes)+1] #produce a probability vector of weights that each individual will contribute to the next generation, based on its genotype. Note that w needs to be inverted to keep proper track because the first value in w is for the AA genotype, which is given by allele combinations 1 and 1. 
      parents<-c(sample(N,prob=fitnessv,replace=T),sample(N,prob=fitnessv,replace=T)) #generate a vector of all parents who will contribute to the next generation (two for each offspring). Some parents contribute more than once (especially if they have selectively advantageous genotypes), and some are not represented at all. 
      nextgen<-sapply(parents,function(x) sample(genotypes[x,],1)) #sample one of the alleles from each parent, each time it contributed
      X[i,j] = sum(nextgen)
    } 
    if (all(nextgen == 1)) {
      numfix <- numfix+1
    }
  }
  probabilities = c(probabilities, numfix/N.sim)
}

probabilities

plot(probabilities, xlab = 'Initial allele frequency', ylab = 'Probability of fixation')



#2

#N = 50 # number of diploid individuals
#N.chrom = 2*N # number of chromosomes
p = .5; q = 1-p
N.gen = 50 # number of generations
N.sim = 100 # number of simulations
w<-c(1,1,1) #a vector of three fitness values for genotype AA, Aa, and aa
w<-(w/max(w))

probabilities = c()
inds = c()
for (i in 1:50) {
  inds = c(inds,i*5)
}
for (individuals in inds) {
  N = individuals # number of diploid individuals
  N.chrom = 2*N # number of chromosomes
  X = array(0, dim=c(N.gen,N.sim)) #each column is a simulation, each row is a generation
  X[1,] = rep(N.chrom*p,N.sim) # initialize number of A1 alleles in first generation
  numunfix <- 0
  for(j in 1:N.sim){
    for(i in 2:N.gen){
      genotypes<-matrix(rbinom(N.chrom,1,prob=X[i-1,j]/N.chrom),N,2) #generate an N x 2 matrix of genotypes, where the probability of drawing the A1 allele is given by its frequency in the previous generation
      fitnessv<-w[3:1][rowSums(genotypes)+1] #produce a probability vector of weights that each individual will contribute to the next generation, based on its genotype. Note that w needs to be inverted to keep proper track because the first value in w is for the AA genotype, which is given by allele combinations 1 and 1. 
      parents<-c(sample(N,prob=fitnessv,replace=T),sample(N,prob=fitnessv,replace=T)) #generate a vector of all parents who will contribute to the next generation (two for each offspring). Some parents contribute more than once (especially if they have selectively advantageous genotypes), and some are not represented at all. 
      nextgen<-sapply(parents,function(x) sample(genotypes[x,],1)) #sample one of the alleles from each parent, each time it contributed
      X[i,j] = sum(nextgen)
    } 
    if (!(all(nextgen == 1) || (all(nextgen == 0)))) {
      numunfix <- numunfix+1
    }
  }
  probabilities = c(probabilities, numunfix/N.sim)
}

probabilities

plot(inds, probabilities, xlab = 'Number of individuals', ylab = 'Proportion of unfixed loci')



#3

N = 5000 # number of diploid individuals  50, 500, 5000
N.chrom = 2*N # number of chromosomes
p = 1/N; q = 1-p
N.gen = 100 # number of generations
N.sim = 100 # number of simulations
#w<-c(1,1,1) #a vector of three fitness values for genotype AA, Aa, and aa
#w<-(w/max(w))

probabilities = c()
s <- 1
h <- 0.5 #0, 0.5, 1
deltas <- 0.005
threshold <- FALSE
prevnumunlost <- 0
while (!threshold) {
  w<-c(1,s+(h*(1-s)),s) #a vector of three fitness values for genotype AA, Aa, and aa
  X = array(0, dim=c(N.gen,N.sim)) #each column is a simulation, each row is a generation
  X[1,] = rep(N.chrom*p,N.sim) # initialize number of A1 alleles in first generation
  numunlost <- 0
  for(j in 1:N.sim){
    for(i in 2:N.gen){
      genotypes<-matrix(rbinom(N.chrom,1,prob=X[i-1,j]/N.chrom),N,2) #generate an N x 2 matrix of genotypes, where the probability of drawing the A1 allele is given by its frequency in the previous generation
      fitnessv<-w[3:1][rowSums(genotypes)+1] #produce a probability vector of weights that each individual will contribute to the next generation, based on its genotype. Note that w needs to be inverted to keep proper track because the first value in w is for the AA genotype, which is given by allele combinations 1 and 1. 
      parents<-c(sample(N,prob=fitnessv,replace=T),sample(N,prob=fitnessv,replace=T)) #generate a vector of all parents who will contribute to the next generation (two for each offspring). Some parents contribute more than once (especially if they have selectively advantageous genotypes), and some are not represented at all. 
      nextgen<-sapply(parents,function(x) sample(genotypes[x,],1)) #sample one of the alleles from each parent, each time it contributed
      X[i,j] = sum(nextgen)
    } 
    if (!(all(nextgen == 0))) {
      numunlost <- numunlost+1
    }
  }
  probabilities = c(probabilities, numunlost/N.sim)
  if ((numunlost/N.sim > 0.25) && (prevnumunlost/N.sim > 0.25)) {
    threshold = TRUE
    s <- s+deltas
  }
  s <- s-deltas
  prevnumunlost <- numunlost
}

probabilities
slabels = c()
for (i in 1:length(probabilities)) {
  slabels = c(slabels, -0.005+0.005*i)
}

plot(slabels, probabilities, main = 'Proportion of loci retaining advantageous allele vs. selection coefficient', xlab = 'Selection coefficient', ylab = 'Proportion of loci with advantageous allele unlost', sub = '5000 individuals, completely additive')