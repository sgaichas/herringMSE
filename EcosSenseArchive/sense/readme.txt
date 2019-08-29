December 14 2016

running script SenseStats.sh produces a vuldisc.csv that does not read into R correctly
grep on each Seed's log file revealed that for Seed6, three groups failed in year 1
all other seeds read in correctly so their csvs can be combined into one file
but for seed 6 read in fails:
expecting only 6 columns, this 7th is misaligned into col 1 thus ruining everything downstream
the fix is to remove Megabenthos- other from line 1230 of vuldisc_6.csv
this is system 1496
now has only 2 groups dying
the reporting is only for the first group anyway
saved as vuldisc_6_1.csv; use this one to combine with others and continue with script in discardStat.r

also, system 2915 (line 2385) had 8 groups fail so obviously this will break it too
removed Megabenthos- other,Shrimp et al.,Small Pelagics- commercial,Small Pelagics- squid,Demersals- benthivores,Demersals- omnivores,


