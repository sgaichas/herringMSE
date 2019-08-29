#!/bin/sh
#Usage: php SenseDiffs.php [difference scenario] [base scenario] > [Output Diff File]

php SenseDiffs.php vulherrdown1.csv vulbase1.csv > vulherrdown_diff.csv
php SenseDiffs.php novulherrdown1.csv novulbase1.csv > novulherrdown_diff.csv

php SenseDiffs.php vulherrup1.csv vulbase1.csv > vulherrup_diff.csv
php SenseDiffs.php novulherrup1.csv novulbase1.csv > novulherrup_diff.csv


#php SenseDiffsP.php vulherrdown1.csv vulbase1.csv > vulherrdown_diffP.csv
#php SenseDiffsP.php novulherrdown1.csv novulbase1.csv > novulherrdown_diffP.csv

#php SenseDiffsP.php vulherrup1.csv vulbase1.csv > vulherrup_diffP.csv
#php SenseDiffsP.php novulherrup1.csv novulbase1.csv > novulherrup_diffP.csv



