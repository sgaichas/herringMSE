#!/bin/sh

#get n kept in all
grep kept *_novul_*log.csv | wc -l > novul.csv
grep kept *_vul_*log.csv | wc -l > vul.csv
sed 's/^/kept,/' novul.csv > novulstat.csv 
sed 's/^/kept,/' vul.csv > vulstat.csv 

#get n discarded in all
grep discarded *_novul_*log.csv | wc -l > novul.csv
grep discarded *_vul_*log.csv | wc -l > vul.csv
sed 's/^/discarded,/' novul.csv >> novulstat.csv 
sed 's/^/discarded,/' vul.csv >> vulstat.csv 

#put all discarded into one file
grep discarded *_novul_*log.csv  > novuldisc.csv
grep discarded *_vul_*log.csv  > vuldisc.csv

#run R script to make pdf histograms
R --slave --args filename "novuldisc" statname "novulstat" < discardStat.r
R --slave --args filename "vuldisc" statname "vulstat" < discardStat.r







