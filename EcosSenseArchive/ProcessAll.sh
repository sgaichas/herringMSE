# >>> Sense Runs processing for GOM scenario runs, December 2016
#
cd sense/base
sh preprocessSense.sh
cp vulbase1.csv ../diffs/vulbase1.csv
cp novulbase1.csv ../diffs/novulbase1.csv

cd ../herrup
sh preprocessSense.sh
cp vulherrup1.csv ../diffs/vulherrup1.csv
cp novulherrup1.csv ../diffs/novulherrup1.csv

cd ../herrdown
sh preprocessSense.sh
cp vulherrdown1.csv ../diffs/vulherrdown1.csv
cp novulherrdown1.csv ../diffs/novulherrdown.csv

cd ../diffs 
sh SenseProcess.sh
#sh SenseProcess_guilds.sh
#sh Rplot_guilds.sh

#cd ../sense
#sh SenseStats.sh

 
