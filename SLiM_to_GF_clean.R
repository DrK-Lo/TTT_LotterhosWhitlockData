#Processing of SLiM output from FullFathomFive_SP_100_###.slim codes from 
# output visualization to Gradient Forest analysis 

require(adegenet)
require(gdm)
require(gradientForest)
require(foreach)
require(doParallel)
require(OutFLANK)
require(pbapply)
require(gdata)
require(data.table)
require(PresenceAbsence)
require(ROCR)
require(modEvA)
require(ggplot2)
require(grid)
require(gridExtra)
require(gtools)
require(cowplot)

##############################################
#Read in the SLiM simulation output files
#and prepare data for GF
##############################################

#setwd("/Users/akijarl/.../SLiM_output/")

#Specify plot title to reflect specific of simulation to be visualized
plotTitle <- "FFF_output_SP_100_1653512648948\n recomb. 1e-6, mig. 0.1, sel. 0.1, Env. rate 0.01"

#Read in the table of M2 allele frequencies over time
fit<-read.table(paste("Fit_output_SP_100_",substr(plotTitle, start=19, stop=31),".txt",sep=""),col.names = paste0("V", 1:301),fill=T) 

#Read in the tables of 1000 neutral alleles sampled at different time points
freq1000<-read.table(paste("Freq_output_SP_100_",substr(plotTitle, start=19, stop=31),"_Gen1000.txt",sep=""))
freq1100<-read.table(paste("Freq_output_SP_100_",substr(plotTitle, start=19, stop=31),"_Gen1100.txt",sep=""))
freq1200<-read.table(paste("Freq_output_SP_100_",substr(plotTitle, start=19, stop=31),"_Gen1200.txt",sep="")) 
freq1300<-read.table(paste("Freq_output_SP_100_",substr(plotTitle, start=19, stop=31),"_Gen1300.txt",sep=""))

#Create Pop. fit column names
fit_nam <- NULL
for(i in 1:100){
  fit_nam <- c(fit_nam,print(paste("P",i,"_fit",sep="")))
}

#Create Pop. freq column names
freq_nam <- NULL
for(i in 1:100){
  freq_nam <- c(freq_nam,print(paste("P",i,"_freq",sep="")))
}

#Create Pop. environment column names
env_nam <- NULL
for(i in 1:100){
  env_nam <- c(env_nam,print(paste("P",i,"_env",sep="")))
}

colnames(fit) <- c("Gen",fit_nam,freq_nam,env_nam)

#Transpose data and add column names
fitt<-data.frame(t(fit))
fitt<-fitt[-1,]

#Create Generation column names
gen_nam <- NULL
for(i in seq(800,1300,10)){
  gen_nam <- c(gen_nam,print(paste("Gen",i,sep="")))
}
colnames(fitt)<-gen_nam

#Add column to track location
fitt$Location <- as.factor(rep(seq(0,9,1),30))

#Define which part of dataset are Fitness values, which are Frequencies, and which are Envrionmental values
fitt$Type <- as.factor(c(rep("Fit",100),rep("Freq",100),rep("Hab",100)))

#Subset M2 allele frequencies, change which Generation you want (i.e. "Gen1000", "Gen1200", etc.)
M2 <- data.frame(fitt[101:200,c("Gen1000","Location")]) #Subset M2 allele frequencies and location by Generation
colnames(M2)[1]<-"M2"

#Generate column names for the null mutants
mut_nam <- NULL
for(i in 1:1000){
  mut_nam <- c(mut_nam,paste("M1",i,sep="_"))
}

#Ascribe column names generated above to the null datasets 
colnames(freq1000)<-mut_nam
colnames(freq1100)<-mut_nam
colnames(freq1200)<-mut_nam
colnames(freq1300)<-mut_nam

#Merge the M2 allele frequency dataset with the appropriate null allele dataframe. Make sure M2 and null allele datasets are from the same time period (if you subset "Gen1200" M2 allele frequency above, specify freq1200 dataset here)
alFreq<-cbind(freq1000,M2) #Keep track of which M2 generation (1000, 1100, 1200, 1300) is being compared
alFreq<-alFreq[,-1002] #Strip off location column

#Subset the environmental variables for the generation you're considering (make sure M2, null, and environmental data are not being compare across generations)
envPop<-data.frame(fitt[201:300,"Gen1000"])
names(envPop) <- "envSelect" #Specify column heading for the environmental data


##############################################
#Set uf GF functions and run data
##############################################
# Build output GF data frames
##############################################
#Set the number of cores to use for parallel processing
cores <- 3

gfR2tab <- function(gfMods.list, alFreqs){
  i=1
  while(is.null(gfMods.list[[i]])){i=i+1}
  tab <- do.call(rbind, gfMods.list)
  vrNm <- rep(row.names(tab)[1:nrow(gfMods.list[[i]])], 
              nrow(tab)/nrow(gfMods.list[[i]]))
  tab <- data.frame(variable=vrNm, tab)
  tab <- dcast(tab, SNPnames~variable, value.var="imps")
  envR2 <- rowSums(data.frame(tab[,-1]))
  R2Tab <- data.frame(tab, envR2=envR2)
  
  # get name of SNP if it has a positive R2
  posR2 <- unlist(lapply(gfMods.list, function(x){
    return(as.character(unique(x[,2])))}))
  
  # Find which loci have R2 < 0 (no GF model for those) & assign R2=0
  negR2 <- !(colnames(alFreqs) %in% posR2)
  negR2 <- colnames(alFreqs)[negR2] 
  
  noGF <- data.frame(matrix(0, nrow=length(negR2), ncol=ncol(R2Tab)))
  colnames(noGF) <- colnames(R2Tab)           
  noGF$SNPnames <- negR2
  
  R2Tab <- rbind(R2Tab, noGF)
  return(R2Tab[order(R2Tab$SNPnames),])}

##############################################
# Chunk to fit GF models to minor allele frequencies at the level of
# populations
# GF is fit to each SNP individually to 
# ease computational / memory burden

##### added by MCF, running all loci in on model ##########
# I add a random variables in hopes of getting some of the 
# gf plotting functions to work (which seem to fail when 
# the model uses only one variable
envRand <- runif(nrow(envPop), -1, 1)
env <- data.frame(envPop, envRand)
gfMod <- gradientForest(data=data.frame(env, alFreq),
                        predictor.vars=colnames(env),
                        response.vars=colnames(alFreq),
                        corr.threshold=0.5, 
                        ntree=500, 
                        trace=T)


plot(gfMod, plot.type="C", common.scale=T)
plot(gfMod, plot.type="S", common.scale=T)
ci.o <- cumimp(gfMod, predictor="envSelect", type="Overall")
#plot(ci.o$x, ci.o$y, type="l")
##### added by MCF, running all loci in on model ##########

cl <- makeCluster(cores)
registerDoParallel(cl)

#Use for each to run in parallel - each SNP modeled by GF on different processor
gfAllele.freq <- foreach(k=1:ncol(alFreq), .verbose=F, .packages=c("gradientForest", "data.table")) %dopar% {
  # get locus k & name with SNP ID
  locus <- data.frame(alFreq[,k])
  names(locus) <- colnames(alFreq)[k]
  
  # here's the GF modeling function
  gfLocus <- gradientForest(data=data.frame(envPop, locus),
                            #gfLocus <- gradientForest(data=data.frame(locPop, locus),
                            predictor.vars=colnames(envPop), 
                            #predictor.vars=colnames(locPop),
                            response.vars=colnames(alFreq)[k], 
                            corr.threshold=0.5, 
                            ntree=500, 
                            trace=F)
  
  # get result is the model is not = NULL (no fit)
  if(!is.null(gfLocus)){
    # cumulative importance values (to plot turnover functions)
    cImp <- cumimp(gfLocus, "envSelect", type="Species")
    cImp <- data.frame(rbindlist(cImp, idcol="allele"))
    # return data frame with allele ID, cumulative importance results,
    # and R2 value for GF model
    data.frame(cImp, r2=gfLocus$result)
  }
}
stopCluster(cl)

# Prep table of importance values for each SNP (rows) 
# and each var (columns)
gfAF <- lapply(gfAllele.freq, function(x){
  if(!is.null(x)){
    ttt <- x[1,c("r2", "allele")]
    return(data.frame(imps=ttt$r2, SNPnames=ttt$allele, row.names = "envSelect",
                      stringsAsFactors=F))}})

# run function to convert to table # then write to file
gfAllele.R2.af <- gfR2tab(gfAF, alFreq)
gfAllele.R2.af <- gfAllele.R2.af[match(mixedsort(gfAllele.R2.af$SNPnames), 
                                       gfAllele.R2.af$SNPnames),]
gfAllele.R2.af$SNPnames <- gsub("S", "", gfAllele.R2.af$SNPnames)

finalOut <- data.frame(SNPnames=gfAllele.R2.af$SNPnames, 
                       #cImp.envSelect.alleleFreq=gfAllele.R2.af$envSelect,
                       r_squared.alleleFreq=gfAllele.R2.af$envR2)

#Save output, make sure to change Gen designation ("Gen1000_","Gen1300_",etc.) before saving 
write.table(finalOut, paste("alFreq_",substr(plotTitle, start=19, stop=31),"Gen1000_gradientforests.Cpval",sep=""), row.names=F)
save(gfAllele.freq, file=paste("alFreq_",substr(plotTitle, start=19, stop=31),"Gen1000_gradientforests.Rdata",sep=""))

#Strip all unused alleles from dataset
impDatList <- gfAllele.freq[unlist(lapply(gfAllele.freq, function(x){!is.null(x)}))]
impDat <- do.call(rbind, impDatList)

#Pull out the cumulative importance for the M2 allele
cImpMAF.sel<-impDat[which(impDat$allele=="M2"),]

#Pot Cumulative Importance plot, make sure to change Generation header in ggtitle
p.imp <- ggplot() + 
  geom_line(data=impDat, aes(x=x, y=y, group=allele), colour=rgb(0,0,0,0.4), lwd=0.5) +
  #facet_grid(. ~ strSel) +
  labs(y="Cumulative Importance", x="Environment") +
  #geom_line(data=cImpMAF.neut, aes(x=x, y=y),
  #          colour=rgb(0,0,1, 0.75), lwd=1) +
  geom_line(data=cImpMAF.sel, aes(x=x, y=y),
            colour=rgb(1,0,0, 0.75), lwd=1) +
  #geom_line(data=cImpMAF, aes(x=x, y=y),
  #          colour=rgb(0,0,0, 0.75), lwd=1) +
  theme(plot.margin = unit(c(1.25,1.25,1.25,1.25), "in")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, colour = "grey60"),
        axis.title.x = element_text(size=24)) +
  theme(axis.text.y = element_text(size = 16, colour = "grey60"),
        axis.title.y = element_text(size=24, vjust=1)) +
  theme(strip.text = element_text(size=16)) +
  ggtitle(paste(plotTitle,"Gen1000")) +
  #scale_x_continuous(limits=c(-4,7))+
  theme(plot.title = element_text(size=14, face="bold.italic"))

p.imp


##############################################
# Chunk to fit GF models to minor allele frequencies 
# All loci
#cl <- makeCluster(cores)
#registerDoParallel(cl)

gfMAF <- gradientForest(data=data.frame(envPop, alFreq),
                        predictor.vars=colnames(envPop), 
                        response.vars=colnames(alFreq), 
                        corr.threshold=0.5, 
                        ntree=500, 
                        trace=T)

cImpMAF <- cumimp(gfMAF, "envSelect", type="Overall")
cImpMAF <- data.frame(x=cImpMAF$x, y=cImpMAF$y)

#stopCluster(cl)

cImpMAF <- data.frame(allele="All", x=cImpMAF$x, y=cImpMAF$y, r2=1, strSel="A")

ggCand <- cImpMAF

p.imp <- ggplot() + geom_line(data=ggCand, aes(x=x, y=y, group=allele),
                              colour=rgb(0,0,0,0.4), lwd=0.5) +
  #facet_grid(. ~ strSel) +
  labs(y="Cumulative Importance", x="Environment") +
  geom_line(data=cImpMAF, aes(x=x, y=y),
            colour=rgb(0,0,1, 0.75), lwd=1) +
  #geom_line(data=cImpMAF.sel, aes(x=x, y=y),
  #         colour=rgb(1,0,0, 0.75), lwd=1) +
  #geom_line(data=cImpMAF, aes(x=x, y=y),
  #          colour=rgb(0,0,0, 0.75), lwd=1) +
  theme(plot.margin = unit(c(1.25,1.25,1.25,1.25), "in")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, colour = "grey60"),
        axis.title.x = element_text(size=24)) +
  theme(axis.text.y = element_text(size = 16, colour = "grey60"),
        axis.title.y = element_text(size=24, vjust=1)) +
  theme(strip.text = element_text(size=16)) +
  ggtitle(paste(plotTitle,", Gen1000")) +
  theme(plot.title = element_text(size=14, face="bold.italic"))

p.imp
