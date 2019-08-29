setwd("/home/sgaichas/Data/Projects/GLOBEC_PRS/GOAagg/diffs")
groupnames <- read.csv("GOA_names.csv")

# CHANGE TWO LINES BELOW to read in file, save filename for plotting
filename <- "vulpinnX5_diff"
outname <- "vul, toothed whales x5 (no pinnipeds)"

require("batch")
parseCommandArgs()

results <- read.csv(paste(filename, ".csv", sep=""))

# get rid of last year with 0s, find new max year, get only last 10 years for equilibrium boxplots
endYr <- max(results$Year)
results <- results[results$Year<endYr,]

maxYr <- max(results$Year)
results_end <- results[results$Year>(maxYr-10),]

# get means of last 10 years, save in new object mean_diff2 
sum_diff <- rowsum(results_end[,6:26], results_end$Species)
mean_diff2 <- sum_diff/length(levels(as.factor(results_end$Year)))

#output all quantiles to csv
diffout <- mean_diff2
rownames(diffout) <- groupnames$Group[1:dim(mean_diff2)[1]]
write.csv(diffout, file=paste(filename, "_meanquants.csv"))

# make a single boxplot of means of last 10 years for all species
pdf(paste(filename, "_box.pdf", sep=""), width=10.4, height=9.5)
par(mfrow=c(1,1))
par(mar=c(14,4,4,2)+0.1)
par(oma=c(2,0,0,0))
plot((1:(dim(mean_diff2)[1])), seq(min(mean_diff2$CI_5), max(mean_diff2$CI_95), length.out=dim(mean_diff2)[1]), 
      type="n", xlab="", ylab="Proportional Difference from Base", axes=F)
abline(h=0.0, lty=2)
mtext(paste(outname))
axis(1, at=1:(dim(mean_diff2)[1]), lab=groupnames$Group[1:dim(mean_diff2)[1]], las=2)
mtext("Species", side=1, outer=T)
axis(2)
for (j in 1:dim(mean_diff2)[1]){
    rect((j-0.2),mean_diff2$CI_25[j],(j+0.2),mean_diff2$CI_75[j], col="grey")
    arrows(j, mean_diff2$CI_25[j], j, mean_diff2$CI_5[j], angle=90, length=0.1)
    arrows(j, mean_diff2$CI_75[j], j, mean_diff2$CI_95[j], angle=90, length=0.1)
    }
dev.off()
    
# make time series plots of differences
pdf(paste(filename, "_ts.pdf", sep=""), width=8, height=10.5)
par(mfrow=c(6,4))
par(mar=c(2,2,3,2)+0.1)
par(oma=c(2,2,2,0))
for (i in 1:max(results$Species)){
    plotmax<-max(results$CI_95[results$Species==i], na.rm=T)
    plotmin<-min(results$CI_5[results$Species==i], na.rm=T)
    plot(results$Year[results$Species==i], results$CI_50[results$Species==i], , ylim=c(plotmin,plotmax), xlab="", ylab="", type="n")
    lines(results$Year[results$Species==i], results$CI_95[results$Species==i], lty=2)
    lines(results$Year[results$Species==i], results$CI_75[results$Species==i])
    lines(results$Year[results$Species==i], results$CI_50[results$Species==i], lwd=2)
    lines(results$Year[results$Species==i], results$CI_25[results$Species==i])
    lines(results$Year[results$Species==i], results$CI_5[results$Species==i], lty=2)
    mtext(groupnames$Group[i])
    }
mtext(paste(outname), side=3, outer=T)
mtext("Year", side=1, outer=T)
mtext("Biomass change (t/km2)", side=2, outer=T)
dev.off()
    
#cleanup workspace, if desired
rm(list=objects()) 

