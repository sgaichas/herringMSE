# >>> 
#
head -n1 GOM_vul_herrdown_Seed0_end.csv > endheader.csv
sed '1,1d' GOM_vul_herrdown_Seed0_end.csv > vul_s0.csv
sed '1,1d' GOM_vul_herrdown_Seed1_end.csv > vul_s1.csv
sed '1,1d' GOM_vul_herrdown_Seed2_end.csv > vul_s2.csv
sed '1,1d' GOM_vul_herrdown_Seed3_end.csv > vul_s3.csv
sed '1,1d' GOM_vul_herrdown_Seed4_end.csv > vul_s4.csv
sed '1,1d' GOM_vul_herrdown_Seed5_end.csv > vul_s5.csv
sed '1,1d' GOM_vul_herrdown_Seed6_end.csv > vul_s6.csv
sed '1,1d' GOM_vul_herrdown_Seed7_end.csv > vul_s7.csv
cat vul_s*.csv  > vulherrdown.csv
cat endheader.csv vulherrdown.csv > vulherrdown1.csv
sed '1,1d' GOM_novul_herrdown_Seed0_end.csv > novul_s0.csv
sed '1,1d' GOM_novul_herrdown_Seed1_end.csv > novul_s1.csv
sed '1,1d' GOM_novul_herrdown_Seed2_end.csv > novul_s2.csv
sed '1,1d' GOM_novul_herrdown_Seed3_end.csv > novul_s3.csv
sed '1,1d' GOM_novul_herrdown_Seed4_end.csv > novul_s4.csv
sed '1,1d' GOM_novul_herrdown_Seed5_end.csv > novul_s5.csv
sed '1,1d' GOM_novul_herrdown_Seed6_end.csv > novul_s6.csv
sed '1,1d' GOM_novul_herrdown_Seed7_end.csv > novul_s7.csv
cat novul_s*.csv > novulherrdown.csv
cat endheader.csv novulherrdown.csv > novulherrdown1.csv
rm vul_s*.csv
rm novul_s*.csv
#php Cons_process_scenarios.php vulherrdown1.csv vulherrdown > vulherrdown_summ.csv
#php Cons_process_scenarios.php novulherrdown1.csv novulherrdown > novulherrdown_summ.csv
#
#cat CGOM_vul_herrdown_Seed1_end.csv CGOM_vul_herrdown_Seed2_end.csv CGOM_vul_herrdown_Seed4_end.csv| grep '^[^Seed]' > vulherrdown_sub.csv
#cat endheader.csv vulherrdown_sub.csv > vulherrdown1_sub.csv
#cat CGOM_novul_herrdown_Seed1_end.csv CGOM_novul_herrdown_Seed2_end.csv CGOM_novul_herrdown_Seed4_end.csv| grep '^[^Seed]' > novulherrdown_sub.csv
#cat endheader.csv novulherrdown_sub.csv > novulherrdown1_sub.csv
#php GOM_process_scenarios.php vulherrdown1_sub.csv vulherrdown_sub > vulherrdown_sub_summ.csv
#php GOMherr_process_scenarios.php novulherrdown1_sub.csv novulherrdown_sub > novulherrdown_sub_summ.csv
                                                                       
