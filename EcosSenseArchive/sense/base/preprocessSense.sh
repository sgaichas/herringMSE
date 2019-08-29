# >>> 
#
head -n1 GOM_vul_baserun_Seed0_end.csv > endheader.csv
sed '1,1d' GOM_vul_baserun_Seed0_end.csv > vul_s0.csv
sed '1,1d' GOM_vul_baserun_Seed1_end.csv > vul_s1.csv
sed '1,1d' GOM_vul_baserun_Seed2_end.csv > vul_s2.csv
sed '1,1d' GOM_vul_baserun_Seed3_end.csv > vul_s3.csv
sed '1,1d' GOM_vul_baserun_Seed4_end.csv > vul_s4.csv
sed '1,1d' GOM_vul_baserun_Seed5_end.csv > vul_s5.csv
sed '1,1d' GOM_vul_baserun_Seed6_end.csv > vul_s6.csv
sed '1,1d' GOM_vul_baserun_Seed7_end.csv > vul_s7.csv
cat vul_s*.csv  > vulbase.csv
cat endheader.csv vulbase.csv > vulbase1.csv
sed '1,1d' GOM_novul_baserun_Seed0_end.csv > novul_s0.csv
sed '1,1d' GOM_novul_baserun_Seed1_end.csv > novul_s1.csv
sed '1,1d' GOM_novul_baserun_Seed2_end.csv > novul_s2.csv
sed '1,1d' GOM_novul_baserun_Seed3_end.csv > novul_s3.csv
sed '1,1d' GOM_novul_baserun_Seed4_end.csv > novul_s4.csv
sed '1,1d' GOM_novul_baserun_Seed5_end.csv > novul_s5.csv
sed '1,1d' GOM_novul_baserun_Seed6_end.csv > novul_s6.csv
sed '1,1d' GOM_novul_baserun_Seed7_end.csv > novul_s7.csv
cat novul_s*.csv > novulbase.csv
cat endheader.csv novulbase.csv > novulbase1.csv
rm vul_s*.csv
rm novul_s*.csv
#php Cons_process_scenarios.php vulbase1.csv vulbase > vulbase_summ.csv
#php Cons_process_scenarios.php novulbase1.csv novulbase > novulbase_summ.csv
#
#cat CGOM_vul_baserun_Seed1_end.csv CGOM_vul_baserun_Seed2_end.csv CGOM_vul_baserun_Seed4_end.csv| grep '^[^Seed]' > vulbase_sub.csv
#cat endheader.csv vulbase_sub.csv > vulbase1_sub.csv
#cat CGOM_novul_baserun_Seed1_end.csv CGOM_novul_baserun_Seed2_end.csv CGOM_novul_baserun_Seed4_end.csv| grep '^[^Seed]' > novulbase_sub.csv
#cat endheader.csv novulbase_sub.csv > novulbase1_sub.csv
#php GOM_process_scenarios.php vulbase1_sub.csv vulbase_sub > vulbase_sub_summ.csv
#php GOM_process_scenarios.php novulbase1_sub.csv novulbase_sub > novulbase_sub_summ.csv
                                                                       
