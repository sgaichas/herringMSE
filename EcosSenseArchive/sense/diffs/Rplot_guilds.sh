#!/bin/sh
#run R script to make pdf boxplots and time series plots for sense outputs

R --slave --args filename "vulherrdown_diff" outname "vul, herring down, B" < diffPlots.r
R --slave --args filename "novulherrdown_diff" outname "novul, herring down, B" < diffPlots.r

R --slave --args filename "vulherrup_diff" outname "vul, herring up, B" < diffPlots.r
R --slave --args filename "novulherrup_diff" outname "novul, herring up, B" < diffPlots.r


R --slave --args filename "vulherrdown_diffP" outname "vul, herring down, P" < diffPlots.r
R --slave --args filename "novulherrdown_diffP" outname "novul, herring down, P" < diffPlots.r

R --slave --args filename "vulherrup_diffP" outname "vul, herring up, P" < diffPlots.r
R --slave --args filename "novulherrshup_diffP" outname "novul, herring up, P" < diffPlots.r




