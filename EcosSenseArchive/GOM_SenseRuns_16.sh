# >>> Sense Runs for Herring MSE, Northeast US shelf
#    Old GOA script for reference below

# Test run--newly compiled code for GOM under Linux, no stanzas, no fitting file
#../src_fitted/ecosim -Mo -Otest/basetest -Fconstant_effort.csv -JGOM_stanzas_in.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0

# First generate base set of ecosystems with and without varying functional response (vul) parameters

../src_fitted/ecosim -Me -Vv0.0 -Osense/GOM_basedat_novul -Fconstant_effort.csv -JGOM_stanzas_in.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S1000 -B50 -R0 -Vs1.0 -Vh1.0
../src_fitted/ecosim -Me -Vv6.0 -Osense/GOM_basedat_vul   -Fconstant_effort.csv -JGOM_stanzas_in.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S4000 -B50 -R0 -Vs1.0 -Vh1.0

# Run R script to evaluate how many food webs kept, who died in not kept ones
# this script has a hard coded working directory in it! point it to the correct place!

cd sense
sh SenseStats.sh
cd ..

# Then perturb the systems both with and without varying vuls, outputting only average relative change in P (last 10 years)

../src_fitted/ecosim -Mp -Vv0.0 -Isense/GOM_basedat_novul  -Osense/perturb/GOM_novul_perturb -Fconstant_effort.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S1000 -B50 -R0 -Vs1.0 -Vh1.0
../src_fitted/ecosim -Mp -Vv6.0 -Isense/GOM_basedat_vul   -Osense/perturb/GOM_vul_perturb   -Fconstant_effort.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S4000 -B50 -R0 -Vs1.0 -Vh1.0


# Test with herring forced down and up drivers (baseline 2000-2003, 5 year slope to either 50% of base or 150% of base until 2050)

#../src_fitted/ecosim -Mo -Otest/basefish_herrup -Fconstant_effort_herringup.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mo -Otest/basefish_herrdown -Fconstant_effort_herringdown.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0

# Then run a base scenario for comparison with herring up and down manipulation runs

../src_fitted/ecosim -Mi -Vv0.0 -Isense/GOM_basedat_novul -Osense/base/GOM_novul_baserun -Fconstant_effort.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S1000 -B50 -R0 -Vs1.0 -Vh1.0
../src_fitted/ecosim -Mi -Vv6.0 -Isense/GOM_basedat_vul   -Osense/base/GOM_vul_baserun   -Fconstant_effort.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S4000 -B50 -R0 -Vs1.0 -Vh1.0

#run with herring up

../src_fitted/ecosim -Mi -Vv0.0 -Isense/GOM_basedat_novul -Osense/herrup/GOM_novul_herrup -Fconstant_effort_herringup.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S1000 -B50 -R0 -Vs1.0 -Vh1.0
../src_fitted/ecosim -Mi -Vv6.0 -Isense/GOM_basedat_vul   -Osense/herrup/GOM_vul_herrup  -Fconstant_effort_herringup.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S4000 -B50 -R0 -Vs1.0 -Vh1.0

#run with herring down

../src_fitted/ecosim -Mi -Vv0.0 -Isense/GOM_basedat_novul -Osense/herrdown/GOM_novul_herrdown -Fconstant_effort_herringdown.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S1000 -B50 -R0 -Vs1.0 -Vh1.0
../src_fitted/ecosim -Mi -Vv6.0 -Isense/GOM_basedat_vul   -Osense/herrdown/GOM_vul_herrdown  -Fconstant_effort_herringdown.csv -JGOM_stanzas_in_SAagg.csv -PGOM_fitVectorsStart.csv -KGOM_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S4000 -B50 -R0 -Vs1.0 -Vh1.0

#Process all runs using scripts in each directory

sh ProcessAll.sh
sh ProcessPerturb.sh







# Run with code compiled for full GOA model, with old GOA stanza file ("corrected" one breaks)
# see changes made to code for sense runs outputting P, no consumption outputs, 10% decrease in M (increase in P)
# in ReadMe.txt under main Beer Belly directory

# Test run--newly compiled code under Linux, old stanzas, no fitting file
#../src_fitted/ecosim -Mo -Otest/basetest -FnoGOA_fitData_09PT.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0

# Test with effort drivers ("real" 1990-2009, then hold 2009 steady till 2059)
#../src_fitted/ecosim -Mo -Osense/basefish -FGOA_baserun_19902059.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0

# Test poll recruitment up in 2009
#../src_fitted/ecosim -Mo -Osense/pollclim/GOA_novul_pollrecup -FGOA_baserun_19902059_pollrecup2009.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0

# Scenarios: same as EBS, poll rec down 30% by 2050 due to climate, fishing up to 40% above base, and both
#../src_fitted/ecosim -Mo -Osense/pollclim/GOA_novul_pollreclow -FGOA_baserun_19902059_pollrecdown.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mo -Osense/pollfish/GOA_novul_pollfishup -FGOA_baserun_19902059_pollfishup.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mo -Osense/pollrecfish/GOA_novul_both -FGOA_baserun_19902059_pollfishuprecdown.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0

# First generate base set of ecosystems with and without varying functional response (vul) parameters

#../src_fitted/ecosim -Me -Vv0.0 -Osense/GOA_basedat_novul -FGOA_baserun_19902059.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Me -Vv6.0 -Osense/GOA_basedat_vul   -FGOA_baserun_19902059.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S40000 -B50 -R0 -Vs1.0 -Vh1.0

# Then perturb the systems both with and without varying vuls, outputting only average relative change in P (last 10 years)

#../src_fitted/ecosim -Mp -Vv0.0 -Isense/GOA_basedat_novul  -Osense/perturb/GOA_novul_perturb -FGOA_baserun_19902059.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mp -Vv6.0 -Isense/GOA_basedat_vul   -Osense/perturb/GOA_vul_perturb   -FGOA_baserun_19902059.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S40000 -B50 -R0 -Vs1.0 -Vh1.0

# Then run a base scenario for comparison with climate and fishing manipulation runs

#../src_fitted/ecosim -Mi -Vv0.0 -Isense/GOA_basedat_novul -Osense/base/GOA_novul_baserun -FGOA_baserun_19902059.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mi -Vv6.0 -Isense/GOA_basedat_vul   -Osense/base/GOA_vul_baserun   -FGOA_baserun_19902059.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S40000 -B50 -R0 -Vs1.0 -Vh1.0

#run with "climate" pollock effect--30% decrease in recruitment by 2050

#../src_fitted/ecosim -Mi -Vv0.0 -Isense/GOA_basedat_novul -Osense/pollclim/GOA_novul_pollreclow -FGOA_baserun_19902059_pollrecdown.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mi -Vv6.0 -Isense/GOA_basedat_vul   -Osense/pollclim/GOA_vul_pollreclow -FGOA_baserun_19902059_pollrecdown.csv-JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S40000 -B50 -R0 -Vs1.0 -Vh1.0

#run with increased fishing, no climate--bring effort up to same as EBS, 40% above base

#../src_fitted/ecosim -Mi -Vv0.0 -Isense/GOA_basedat_novul -Osense/pollfish/GOA_novul_pollfishup -FGOA_baserun_19902059_pollfishup.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mi -Vv6.0 -Isense/GOA_basedat_vul   -Osense/pollfish/GOA_vul_pollfishup -FGOA_baserun_19902059_pollfishup.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S40000 -B50 -R0 -Vs1.0 -Vh1.0

#run with both

#../src_fitted/ecosim -Mi -Vv0.0 -Isense/GOA_basedat_novul -Osense/pollrecfish/GOA_novul_both -FGOA_baserun_19902059_pollfishuprecdown.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S10000 -B50 -R0 -Vs1.0 -Vh1.0
#../src_fitted/ecosim -Mi -Vv6.0 -Isense/GOA_basedat_vul   -Osense/pollrecfish/GOA_vul_both -FGOA_baserun_19902059_pollfishuprecdown.csv -JGOA_stanzas_in_SAagg.csv -PGOA_fitVectorsStart.csv -KGOA_nofit_fitControl.csv -Cno_climate_in.csv -DnoDietFits.csv -S40000 -B50 -R0 -Vs1.0 -Vh1.0


#Process all runs using scripts in each directory
#
#sh ProcessAll.sh
