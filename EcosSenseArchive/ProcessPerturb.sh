# >>> Sense Runs processing for 10% decrease in M (increase in P) perturbation runs, June 2012
# modified from June 2008 "file handling for SenseJune11.doc"
#
cd sense/perturb

head -n1 GOM_novul_perturb_Seed0_actor1.csv > GOM_actorheader.csv

#for broken files, strip off last incomplete line; comment out until cat statement for complete files
#for f in *_novul_perturb_Seed*_actor*.csv
#do
#	sed '$d' < $f > temp
#	mv $f $f.old
#	mv temp $f
#done
cat *novul_perturb_Seed*_actor*.csv | grep 'Diff' > Diffs_novul.csv
cat *_vul_perturb_Seed*_actor*.csv | grep 'Diff' > Diffs_vul.csv

cat GOM_actorheader.csv Diffs_novul.csv > GOM_Diffs_novul.csv
cat GOM_actorheader.csv Diffs_vul.csv > GOM_Diffs_vul.csv
cp GOM_Diffs_novul.csv ../../results/GOM_Diffs_novul.csv
cp GOM_Diffs_vul.csv ../../results/GOM_Diffs_vul.csv

cd ../../results

#run r script ActorQuants.r to get species specific quantiles of diff from base--automate this
R --slave --args filename "GOM_Diffs_novul" groupnames "GOMactors" < ActorQuants.r
R --slave --args filename "GOM_Diffs_vul" groupnames "GOMactors" < ActorQuants.r

#all below now done within the r script ActorQuants.R

#head -n1 GOM_Diffs_vul_HMS_quants.csv > diffheader.csv
#for file in GOM_Diffs_novul_*.csv; do sed '1d' "$f"; done > allnovul_diffquants.csv
#cat *_novul_*_quants.csv > allnovul_diffquants.csv
#cat *_vul_*_quants.csv > allvul_diffquants.csv

#cat diffheader.csv, allnovul_diffquants.csv > GOM_allnovul_diffquants.csv
#cat diffheader.csv, allvul_diffquants.csv > GOM_allvul_diffquants.csv


