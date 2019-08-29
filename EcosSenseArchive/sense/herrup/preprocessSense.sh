# >>> 
#
head -n1 GOM_vul_herrup_Seed0_end.csv > endheader.csv
sed '1,1d' GOM_vul_herrup_Seed0_end.csv > vul_s0.csv
sed '1,1d' GOM_vul_herrup_Seed1_end.csv > vul_s1.csv
sed '1,1d' GOM_vul_herrup_Seed2_end.csv > vul_s2.csv
sed '1,1d' GOM_vul_herrup_Seed3_end.csv > vul_s3.csv
sed '1,1d' GOM_vul_herrup_Seed4_end.csv > vul_s4.csv
sed '1,1d' GOM_vul_herrup_Seed5_end.csv > vul_s5.csv
sed '1,1d' GOM_vul_herrup_Seed6_end.csv > vul_s6.csv
sed '1,1d' GOM_vul_herrup_Seed7_end.csv > vul_s7.csv
cat vul_s*.csv  > vulherrup.csv
cat endheader.csv vulherrup.csv > vulherrup1.csv
sed '1,1d' GOM_novul_herrup_Seed0_end.csv > novul_s0.csv
sed '1,1d' GOM_novul_herrup_Seed1_end.csv > novul_s1.csv
sed '1,1d' GOM_novul_herrup_Seed2_end.csv > novul_s2.csv
sed '1,1d' GOM_novul_herrup_Seed3_end.csv > novul_s3.csv
sed '1,1d' GOM_novul_herrup_Seed4_end.csv > novul_s4.csv
sed '1,1d' GOM_novul_herrup_Seed5_end.csv > novul_s5.csv
sed '1,1d' GOM_novul_herrup_Seed6_end.csv > novul_s6.csv
sed '1,1d' GOM_novul_herrup_Seed7_end.csv > novul_s7.csv
cat novul_s*.csv > novulherrup.csv
cat endheader.csv novulherrup.csv > novulherrup1.csv
rm vul_s*.csv
rm novul_s*.csv
#php Cons_process_scenarios.php vulherrup1.csv vulherrup > vulherrup_summ.csv
#php Cons_process_scenarios.php novulherrup1.csv novulherrup > novulherrup_summ.csv
#
#cat CGOM_vul_herrup_Seed1_end.csv CGOM_vul_herrup_Seed2_end.csv CGOM_vul_herrup_Seed4_end.csv| grep '^[^Seed]' > vulherrup_sub.csv
#cat endheader.csv vulherrup_sub.csv > vulherrup1_sub.csv
#cat CGOM_novul_herrup_Seed1_end.csv CGOM_novul_herrup_Seed2_end.csv CGOM_novul_herrup_Seed4_end.csv| grep '^[^Seed]' > novulherrup_sub.csv
#cat endheader.csv novulherrup_sub.csv > novulherrup1_sub.csv
#php GOM_process_scenarios.php vulherrup1_sub.csv vulherrup_sub > vulherrup_sub_summ.csv
#php GOMherr_process_scenarios.php novulherrup1_sub.csv novulherrup_sub > novulherrup_sub_summ.csv
                                                                       
