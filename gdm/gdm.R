# load libraries and data tables -----------------------------------------------
library(gdm)

# Fst data pre climate change
preFst <- read.csv("/Volumes/Samsung_T5/Projects/testTheTests/1674984747052_Pre_FST_allele.csv")
#preFst <- read.csv("/Volumes/Samsung_T5/Projects/testTheTests/Pre_FST.csv")
preFst[preFst<0] <- 0
#names(preFst)[1] <- "popID"
# rescale 0-1 to assist with model convergence
#preFst[,-1] <- preFst[,-1]/(max(preFst[,-1]))


# Fst data post climate change
#postFst <- read.csv("/Volumes/Samsung_T5/Projects/testTheTests/Post_FST.csv")
#names(postFst)[1] <- "popID"
# rescale 0-1 to assist with model convergence
#postFst[,-1] <- postFst[,-1]/(max(postFst[,-1]))


# env data pre / post climate change & offset
#offset <- read.csv("/Volumes/Samsung_T5/Projects/testTheTests/FST_bT1T2-1.csv")
dat <- read.csv("/Volumes/Samsung_T5/Projects/testTheTests/1674984747052_sum_Pop_subset.csv")

# offset
offset <- dat$F_ST_M2_bef.aft.

# env at time=1
#envT1 <- read.table("/Volumes/Samsung_T5/Projects/testTheTests/T1_Env.txt")
envT1 <- dat$Env_before
envT1rast <- raster(nrows=10, ncols=10, res=1, xmn=0, xmx=10, ymn=0, ymx=10)
envT1rast[] <- envT1
plot(envT1rast)

# env at time=2 (after climate change)
#envT2 <- read.table("/Volumes/Samsung_T5/Projects/testTheTests/T2_Env.txt")
envT2 <- dat$Env_after
envT2rast <- raster(nrows=10, ncols=10, res=1, xmn=0, xmx=10, ymn=0, ymx=10)
envT2rast[] <- envT2
plot(envT2rast)

# generate x, y positions for each population
xy <- xyFromCell(envT1rast, 1:100)#expand.grid(x=1:10, y=10:1)
#xy <- data.frame(popID=1:100, xy)#expand.grid(x=1:10, y=10:1))
#xy <- xy[match(preFst$popID, xy$popID),]
#preFst <- cbind(xy, preFst[,-1]) 

# format data frames for model prep
envT1 <- data.frame(popID = 1:100, xy, env=envT1)
envT2 <- data.frame(popID = 1:100, xy, env=envT2)


# gdm stuff --------------------------------------------------------------------
# format Fst and env for GDM
gdmData <- formatsitepair(bioData=preFst,
                          bioFormat = 3,
                          siteColumn = "popID",
                          XColumn = "x",
                          YColumn = "y",
                          predData = envT1)

# fit & plot model
gdmData <- subset(gdmData, distance>0.05) # only way to get model to fit is remove small values...
mod <- gdm(gdmData, geo=T)
plot(mod)
dev.off()

# predict genetic offset between time=1 and time=2
predGDM <- predict(mod, data=envT1rast, predRasts=envT2rast, time=T)
plot(predGDM) # a raster map of offset
predOffset <- extract(predGDM, envT1[,c("x", "y")])
plot(predOffset, offset)


