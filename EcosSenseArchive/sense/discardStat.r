setwd("/home/sgaichas/Data/Projects/MSE/HerringMSE/GOM/sense/")

filename <- "vuldisc"
statname <- "vulstat"

require("batch")
parseCommandArgs()

discarded <- read.csv(paste(filename, ".csv", sep=""), header=F)
stat <- read.csv(paste(statname, ".csv", sep=""), header=F)

nkept <- stat[1,2]

groupname <- levels(discarded[,5])

dat <- discarded[,4:5]
colnames(dat) <- c("year","groupdied")
dat <- dat[!is.na(dat$year),]
nfail <- dim(dat)[1]

#summary(dat$groupdied)

pdf(paste(filename, "_hist.pdf", sep=""), width=9, height=10.5)
par(mfrow=c(5,4))
par(mar=c(2,2,5,2)+0.1)
par(oma=c(2,2,2,0))
hist(dat$year, breaks=c(0:50), xlim=c(0,50), main=paste((nkept + nfail), "systems drawn", "\n", nfail, "total failing"), xlab="")
for (i in 1:length(groupname)){
    nspfail <- length(dat$groupdied[dat$groupdied==groupname[i]])
    hist(dat$year[dat$groupdied==groupname[i]], breaks=c(0:50), xlim=c(0,50), main=paste(groupname[i],"\n", nspfail, "failures"), xlab="")
    }
mtext("Year failing", side=1, outer=T)
mtext("Frequency", side=2, outer=T)
dev.off()

#cleanup workspace, if desired
rm(list=objects()) 
