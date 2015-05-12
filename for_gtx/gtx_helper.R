library(gtx)

merpfile = "for_gtx/merp_for_gtx.txt";

###read IV data file from MeRP pipeline
mymerp <- read.table(file=merpfile, header=T)

###initial scan
#calculate the risk score and causal effect estimate within gtx
raw_score_sum <- grs.summary(mymerp$BETA_TRAIT,mymerp$BETA_ENDPOINT,mymerp$SE_ENDPOINT,n=NA)

#make plot, check for outliers
grs.plot(mymerp$BETA_TRAIT,mymerp$BETA_ENDPOINT,mymerp$SE_ENDPOINT)


####heterogeneity filtering
#filter IV based on heterogeneity stats
keep_list <- grs.filter.Qrs(mymerp$BETA_TRAIT,mymerp$BETA_ENDPOINT,mymerp$SE_ENDPOINT, p.thresh = 0.05)

#created filtered IV file
filt_mymerp <-mymerp[keep_list,]

#redo score analysis with filtered instrument
filt_score_sum <- grs.summary(filt_mymerp$BETA_TRAIT,filt_mymerp$BETA_ENDPOINT,filt_mymerp$SE_ENDPOINT,n=NA)

#replot as sanity check
grs.plot(filt_mymerp$BETA_TRAIT,filt_mymerp$BETA_ENDPOINT,filt_mymerp$SE_ENDPOINT)