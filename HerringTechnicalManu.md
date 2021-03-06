---
title: 'The dream and the reality: meeting decision-making time frames while incorporating
  ecosystem and economic models into management strategy evaluation'
author: Jonathan Deroba, Sarah Gaichas, Min-Yang Lee, Rachel G. Feeney, Deirdre Boelke,
  Brian J. Irwin
#date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  html_document: null
  pdf_document:
    includes:
      in_header: preamble-latex.tex
    keep_tex: yes
    pandoc_args: --latex-engine=xelatex
  word_document: null
bibliography: HerringMSE.bib
csl: canadian-journal-of-fisheries-and-aquatic-sciences.csl
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, fig.align='center')
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)
```

#Introduction

Management strategy evaluation (MSE) uses simulation to evaluate the trade-offs resulting from alternative management options in the face of uncertainty (Punt et al. 2014).  MSEs require time, however, for stakeholder input, data collection, and model development (Butterworth 2007; Punt et al., 2014).  As such, the process can take much longer than "traditional" management time frames (Butterworth 2007).  The development time is also likely to lengthen when explicit ecosystem, multi-species, or socioeconomic considerations are desired because the data and modeling needs, and subsequent uncertainties, are all greater than in a single species approach.  This manuscript chronicles the development of an MSE done on a truncated timetable (~12 months) required to meet management time frames.

In January 2016, the New England Fishery Management Council  (NEFMC), the political body responsible for federally managed species in the northeast US, approved the conduct of an MSE to evaluate harvest control rules (HCRs) for Atlantic herring  (hereafter herring) *Clupea harengus*.  The primary goal of this MSE was to evaluate the effects of alternative HCRs on the herring fishery, predators of herring, and the human environment.   The MSE was also to include a public stakeholder process.   The NEFMC desired to have results from the MSE ready to inform fisheries management decisions within one year, which left little time to develop the technical aspects of the MSE, especially in the context of multiple public stakeholder workshops.  A particularly challenging aspect of the time frame was deciding what technical aspects of the MSE (e.g., operating models) could be compromised (i.e., perhaps not ideal from a scientific or best-practices standpoint) while ensuring the results would still accurately portray the trade-offs among objectives and remain relevant for decision-making. The objectives of this manuscript were to:

1.	Evaluate the relative performance of HCRs at meeting herring fishery objectives, including those related to predators of herring, as informed by stakeholder input, and,

2.	Discuss our approach to developing an MSE on a relatively truncated timetable in order to meet management time frames, and identify the lessons learned throughout the process, especially as they relate to using MSE as a tool to advance an ecosystem based approach to management (Plag?nyi et al., 2014).


#Methods
**Herring** 
 
*Basics.--* An MSE was developed specific to Gulf of Maine - Georges Bank Atlantic herring.  The MSE was a modified version of that used in Deroba (2014), and symbols were largely consistent with Deroba (2014; Table 1). The MSE was based on an age-structured simulation that considered fish from age-1 through age-8+ (age-8 and older), which is consistent with the age ranges used in the 2012 and 2015 Atlantic herring stock assessments (NEFSC 2012; Deroba 2015).  The abundances at age in year one of all simulations equaled the equilibrium abundances produced by the fishing mortality rate that would reduce the population to 40% of $SSB_{F=0}$.  Abundance in each subsequent age and year was calculated assuming that fish died exponentially according to an age and year specific total instantaneous mortality rate (Table XX T1-T2).

Recruitment followed Beverton-Holt dynamics (Francis 1992; Table XX T3-T5).  The variance of recruitment process errors ($\sigma_R^2$) equaled 0.36 and the degree of autocorrelation ($\omega$) equaled 0.1, which are values consistent with recruitment estimates from a recent Atlantic herring stock assessment (Deroba 2015).

*Assessment Error.--* A stock assessment was approximated (i.e., assessment errors) similar to Punt et al. (2008) and Deroba (2014).  Assessment error was modeled as a year-specific lognormal random deviation common to all ages, with first-order autocorrelation and a term that created the option to include bias (Table XX T6-T7).  The variance of assessment errors ($\sigma_\phi^2$ ) equaled 0.05 and autocorrelation ($\vartheta$) equaled 0.7.  Rho ($\rho$) allowed for the inclusion of bias in the assessed value of abundance (see below; Deroba 2014).  Assessed spawning stock biomass ($\widehat{SSB}_y$) was calculated similarly to $SSB_y$ except with $N_{a,y}$ replaced with $\widehat{N}_{a,y}$ (Table XX T5), and assessed total biomass ($\widehat{B}_y$) was calculated as the sum across ages of the product of $\widehat{N}_{a,y}$ and $W_a$.

*Operating Models.--* The stakeholder workshops identified uncertainties about herring life history traits and stock assessment, and the effect of some of these uncertainties on harvest control rule performance was evaluated by simulating the control rules for each of eight operating models (Table 2; Figures 1-2).  The uncertainties addressed by the eight operating models included: Atlantic herring natural mortality and recruitment , Atlantic herring weight-at-age, and possible bias in the stock assessment beyond the unbiased measurement error ($\epsilon_{\phi y}$). 

The specific values used in the operating models for each of the uncertainties were premised on data used in recent stock assessments or estimates from fits of stock assessment models (Deroba 2015).  Natural mortality in recent stock assessments has varied among ages and years, with $M$ being higher during 1996-2014 than in previous years (NEFSC 2012; Deroba 2015).  Natural mortality, however, has also been identified as an uncertainty in the stock assessments and sensitivity runs have been conducted without higher $M$ during 1996-2014, such that $M$ was constant among years (NEFSC 2012; Deroba 2015).  To capture uncertainty in $M$ in the MSE, operating models were run with either relatively high or low $M_a$ (Table 2; Figure 1).  Relatively high $M_a$ values equaled the age-specific natural mortality rates used for the years 1996-2014 in the stock assessment.  Relatively low $M_a$ values in the MSE equaled the age-specific natural mortality rates used for the years 1965-1995 in the stock assessment.  In the MSE, $M_a$ was always time invariant.

Uncertainty in estimates of stock-recruit parameters were represented in the MSE by using the parameters estimated by stock assessments fit with and without the higher $M$ during 1996-2014.  Stock assessment fits with higher $M$ during 1996-2014 produced estimates of steepness and $SSB_{F=0}$ that were lower than in stock assessment fits without higher $M$ during 1996-2014 (Table 3; Figure 1).  Thus, operating models with relatively high $M_a$ always had relatively low steepness and $SSB_{F=0}$, and the opposite held with relatively low $M_a$ (Table 2).

Uncertainty in Atlantic herring size-at-age was accounted for by having operating models with either fast or slow growth (i.e., weights-at-age; Table 2; Figure 3).  Atlantic herring weight-at-age generally declined from the mid-1980s through the mid-1990s, and has been relatively stable since.  Reasons for the decline are speculative and no causal relationships have been established.  Thus, fast growth operating models had weights-at-age that equaled the January 1 weights-at-age from the most recent stock assessment averaged over the years 1976-1985, while the slow growth operating models averaged over the years 2005-2014 (Deroba 2015).  In the MSE, weight-at-age was always time invariant.
Differences in $M$, stock-recruit parameters, and weights-at-age led to differences in unfished and $MSY$ reference points among operating models (Table 3).  The effect of $M$ and stock-recruit parameters was larger than the effect of differences in weight-at-age (Table 3).

To address concerns about possible stock assessment bias, operating models with and without a positive bias were included.  In operating models without bias, $\rho$=0 and the only assessment error was that caused by the unbiased measurement errors ($\epsilon_{\phi y}$).  In operating models with bias, $\rho$=0.6, which was based on the degree of retrospective pattern in $SSB$ from the most recent stock assessment (Deroba 2015).

*Harvest Control Rules.--* Several  basic control rules were evaluated, including a biomass based control rule (Katsukawa 2004), a constant catch rule, and a conditional constant catch rule (Figure 3; Clark and Hare 2004; Deroba and Bence 2012).  The biomass based control rule was defined by three parameters: the proportion ($\psi$) of $F_{MSY}$ that dictates the maximum desired fishing mortality rate ($\tilde{F}$), an upper SSB threshold ($SSB_{up}$), and a lower SSB threshold ($SSB_{low}$).  The $\tilde{F}$ equaled the maximum when $\widehat{SSB}$ was above the upper threshold, declined linearly between the upper and lower thresholds, and equaled zero below the lower threshold:


$\tilde{F}_y= \bigg\{ {\substack{F_{MSY} \psi \; \; if \; \widehat{SSB}_y >= SSB_{up}\\(F_{MSY} \psi) \frac{\widehat{SSB}_y - SSB_{low}}{SSB_{up}-SSB_{low}} \; \; if \; SSB_{low}<\widehat{SSB}_y<SSB_{up} \\0 \; \; if \; \widehat{SSB}_y<=SSB_{low}}}$


The $\tilde{F}_y$ was then used to set a quota in year y + 1 (Table XX T8).  $\tilde{F}_{ay}$ equaled $\tilde{F}_y$ times $S_a$, and $S_a$ was time and simulation invariant selectivity at age equal to the values for the mobile gear fishery reported in Deroba (2015; Table 1).  $\tilde{F}_y$ was used to set a quota in the following year to approximate the practice of using projections based on an assessment using data through year y - 1 to set quotas in the following year(s).  Furthermore, although $\tilde{F}_y$ was set using $\widehat{SSB}_y$, the quota was based on $\widehat{B}_y$ because the fishery selects some immature ages.  The fully selected fishing mortality rate that would remove the quota from the true population ($\bar{F}_y$) was found using Newton-Raphson iterations.
Several variations of the biomass based rule were also evaluated.  These variations included applying the control rule annually, using the same quota for three year blocks such that the control rule is applied every fourth year (i.e., $Q_{y+1}=Q_{y+2}=Q_{y+3}$), using the same quota for 5 year blocks, and using the same quota for three year blocks but restricting the change in the quota to 15% in either direction when the control rule was reapplied in the fourth year.  Thus, four variants of the biomass based control rule were evaluated: 1) annual application, 2) three year blocks, 3) five year blocks, and 4) 3 year blocks with a 15% restriction.

For each biomass based control rule variant, a range of values for the three parameters defining the control rule were evaluated.  The proportion ($\psi$) of $F_{MSY}$ that dictates the maximum desired fishing mortality rate was varied from 0.1$F_{MSY}$ to 1.0$F_{MSY}$ in increments of 0.1, while the upper and lower SSB threshold parameters ($SSB_{up} \; SSB_{low}$) were varied from 0.0$SSB_{MSY}$ to 4$SSB_{MSY}$ but with inconsistent increments (i.e., 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3, 3.5, 4).  The full factorial of combinations for the three biomass based control rule parameters produced 1,360 shapes (note $SSB_{low}$ must be <= $SSB_{up}$) and each of these shapes was evaluated for each of the four biomass based control rule variants described above.

The constant catch control rule is defined by one parameter, a desired constant catch (i.e., quota) amount (Figure 3).  The constant catch amounts were varied from 0.1$MSY$ to 1.0$MSY$ in increments of 0.1.  

The conditional constant catch rule used a constant desired catch amount unless removing that desired catch from the assessed biomass caused the fully selected fishing mortality rate to exceed a pre-determined maximum, in which case the desired catch was set to the value produced by applying the maximum fully selected fishing mortality rate to the assessed biomass (Figure 3).  Thus, the conditional constant catch rule has two policy parameters: a desired constant catch amount, and a maximum fishing mortality rate.  The constant catch amounts were varied from 0.1$MSY$ to 1.0$MSY$ in increments of 0.1, while the maximum fishing mortality rate equaled 0.5$F_{MSY}$.  When the maximum fishing mortality rate portion the conditional constant catch rule was invoked, a quota was set in the same manner as when $\widehat{SSB}_y >= SSB_{up}$ in the biomass based control rule described above.

*Implementation Error.--* Implementation errors were also included in a similar way as in Punt et al. (2008) and Deroba and Bence (2012), as year-specific lognormal random deviations (Table XX T9).  The variance of implementation errors ($\sigma_\theta^2$) equaled 0.001.

**Predators**

There are two components of predator modeling for the herring MSE: a predator population model, and a herring-predator relationship model to link herring with predator populations. Here, we give an overview of the modeling process, and we describe the decisions made in parameterizing individual predator models and herring-predator relationships in the following sections. The overall population in numbers for each predator each year $N_{y}$ is modeled with a delay-difference function: 
\begin{equation}
N_{y+1} = N_{y}S_{y} +  R_{y+1} \label{delaydiffN_equation},
\end{equation}
where annual predator survival $S_{y}$ is based on annual natural mortality $v$ and exploitation $u$ 
\begin{equation}
S_{y} =  (1-v_{y})(1-u) \label{survival_equation},
\end{equation}
and annual recruitment $R_{y}$ (delayed until recruitment age a) is a Beverton-Holt function:
\begin{equation}
R_{y+a} = \frac{\alpha B_{y}}{\beta + B_{y}} \label{delaydiffrec_equation}.
\end{equation}

Predator recruitment parameters are defined with steepness = $h$, unfished recruitment $R_{F=0}$, and unfished spawning biomass $B_{F=0}$ as
\begin{equation}
\alpha = \frac{4hR_{F=0}}{5h-1} \label{predalpha_equation}, and
\end{equation}

\begin{equation}
\beta = \frac{(B_{F=0}/R_{F=0})((1-h)/(4h))}{(5h-1)/(4hR_{F=0})} \label{predbeta_equation}
\end{equation}. 

Predator population biomass is defined with Ford-Walford plot intercept ($FWint$) and slope ($FWslope$) growth parameters
\begin{equation}
B_{y+1} = S_{y} (FWint N_{y} + FWslope B_{y}) + FWint R_{y+1} \label{delaydiffB_equation}
\end{equation}

Parameterizing this model requires specification of the stock-recruitment relationship (steepness h and unfished spawning stock size in numbers or biomass), the natural mortality rate, the fishing mortality (exploitation) rate, the initial population size, and the weight at age of the predator (Ford-Walford plot intercept and slope parameters). For each predator, population parameters were derived from different sources (Tab. 1).

Table 1. Predator population model specification and parameter sources
```{r popsources, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'}
tabl <- "
|                      | Highly migratory     | Seabird              | Groundfish           | Marine mammal        |
|:---------------------|:---------------------|:---------------------|:---------------------|:---------------------|
| Stakeholder preferred species | Bluefin tuna | Common tern | not specified | not specified |
| Species modeled     | Bluefin tuna (western Atlantic stock) | Common tern (Gulf of Maine colonies as defined by the GOM Seabird Working Group) | Spiny dogfish (GOM and GB cod stocks also examined) | none, data limited (Minke & humpback whales, harbor porpoise, harbor seal examined) |
| Stock-recruitment   | Current assessment and literature | Derived from observations | Current assessment and literature | No time series data for our region | 
| Natural mortality | Current assessment | Literature | Current assessment | Derivable from assessment? |
| Fishing mortality | Current assessment | n/a | Current assessment | Derivable from assessment? |
| Initial population | Current assessment | Derived from observations | Current assessment | Derivable from assessment? |
| Weight at age | Literature | Literature | Literature | Literature | 
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion
```

Predator population models were based on either the most recent stock assessment for the predator or from observational data from the Northeast US shelf. Herring-predator relationships were based on either peer-reviewed literature or observational data specific to the Northeast US shelf. We did not include process or observation error in any of these modeled relationships. This is obviously unrealistic, but the primary obective of the herring MSE is to evaluate the effect of herring management on predators. Leaving out variability driven by anything other than herring is intended to clarify the effect of herring managment. 

To develop the herring-predator relationship model, specific herring population characteristics (e.g. total abundance or biomass, or abundance/biomass of certain ages or sizes) were related to either predator growth, predator reproduction, or predator survival. Our aim was to use information specific to the Northeast US shelf ecosystem, either from peer-reviewed literature, from observations, or a combination. 

In general, if support for a relationship between herring and predator recruitment was evident, it was modeled as a predator recruitment multiplier based on the herring population ($Hpop_{y}$) relative to a specified threshold $Hthresh$:
\begin{equation}
R_{y+a} = \frac{\alpha B_{y}}{\beta + B_{y}} * \frac{\gamma(Hpop_{y}/Hthresh)}{(\gamma-1)+(Hpop_{y}/Hthresh)} \label{recwithherring_equation}, 
\end{equation}

where $\gamma$ > 1 links herring population size relative to the threshold level to predator recruitment.

If a relationship between predator growth and herring population size was evident, annual changes in growth were modeled by modifying either the Ford-Walford intercept ($AnnAlpha$) or slope ($AnnRho$):

\begin{equation}
B_{y+1} = S_{y} (AnnAlpha_{y} N_{y} + FWslope B_{y}) + AnnAlpha_{y}R_{y+1} \label{delaydiffB_equation}, or
\end{equation}
\begin{equation}
B_{y+1} = S_{y} (FWint N_{y} + AnnRho_{y} B_{y}) + FWint R_{y+1} \label{delaydiffB_equation}.
\end{equation}

Finally, herring population size $Hpop_{y}$ could be related to predator survival using a multiplier on constant predator annual natural mortality $v$: 

\begin{equation}
v_{y} =  v e ^ {-(\frac{Hpop_{y}}{Hpop_{F=0}})\delta} \label{varmort_equation},
\end{equation}
where 0 < $\delta$ <1 links herring population size to predator survival.

After specifying the population model parameters and herring-predator relationship, we applied the [@hilborn_quantitative_2003] equilibruim calculation for the delay difference model with F=0 to get the unfished spawners per recruit ratio. This ratio was then used in a new equilibruim calculation with the current predator exploitation rate to estimate Beverton-Holt stock recruitment parameters, equilibrium recruitment and equilibrium individual weight under exploitation. Then, each model was run forward for 150 years with output from the herring operating model specifying the herring population characteristics. 

## Highly migratory species
Bluefin tuna were identified at the stakeholder workshop as the recommended highly migratory herring predator. 

### Tuna population model
Western Atlantic bluefin tuna population parameters were drawn from the 2014 stock assessment [@iccat_report_2015], the growth curve from [@restrepo_updated_2010], and recruitment parameters from a detailed examination of alternative stock recruit relationships [@porch_making_2016]. Ultimately, the “low recruitment” scenario was selected to represent bluefin tuna productivity in the Gulf of Maine, which defines Bmsy as 13,226 t and therefore affects measures of status relative to Bmsy. Continuation of the current tuna fishing strategy (F<0.5Fmsy under the low recruitment scenario) is assumed. All predator population model parameters are listed in Table 2. 

### Herring-tuna relationship model
Tuna diets are variable depending on location and timing of foraging [@chase_differences_2002; @golet_changes_2013; @logan_diet_2015; @golet_paradox_2015], but for the purposes of this analysis, we assumed that herring is an important enough prey of tuna to impact tuna growth in the Northeast US shelf ecosystem. A relationship between bluefin tuna growth and herring average weight was implemented based on information and methods in @golet_paradox_2015. The relationship between tuna condition anomaly (defined as proportional departures from the weight-at-length relationship used in the assessment) and average weight of tuna-prey-sized herring ($Havgwt_{y}$, herring >180 mm collected from commercial herring fisheries) was modeled as a generalized logistic function with lower and upper bounds on tuna growth parameters:
\begin{equation}
AnnAlpha_{y} = (0.9 FWint) + \frac{(1.1 FWint) - (0.9 FWint)}{1+e^{(1-\lambda)*(100(Havgwt_{y}-T)/T)}} \label{tunagrow_equation},
\end{equation} 
where $\lambda$ > 1 links herring average weight anomalies to tuna growth.
   
The inflection point of $T$ = 0.15 kg average weight matches the 0 tuna weight anomaly in @golet_paradox_2015 (p. 186, Fig 2C), and upper and lower bounds were determined by estimating the growth intercept with weight at age 10% higher or lower, respectively from the average weight at age obtained by applying the length to weight conversion reported in the 2014 stock assessment [@iccat_report_2015] to the length at age estimated from the @restrepo_updated_2010 growth curve (Fig \ref{herringtuna}). When included in the model with $\lambda$ = 1.1 in equation \ref{tunagrow_equation}, the simulated variation in tuna weight at age covered the observed range reported in @golet_paradox_2015. 

```{r, fig.cap="Modeled herring average weight-tuna growth relationship \\label{herringtuna}", out.width='300pt', message = FALSE, warning = FALSE}
# from MSE_TunaLink.R
FWalpha       <-  0.020605  #ford-walford plot intercept parameter 
FWrho         <-  0.9675  #ford-walford plot slope parameter
BFTGrowThresh <-  0.15
preypredgrow  <-  1.1 #strengh of effect of prey N on predator growth, >=1, 1=no effect

preyAvgWt     <-  seq(0,0.3,by=0.01)
AnnualAlpha   <-  c()

AnnualAlpha<-(0.9*FWalpha) + ((1.1*FWalpha) - (0.9*FWalpha))/(1+exp((1-preypredgrow)*(100*(preyAvgWt-BFTGrowThresh)/BFTGrowThresh))) #alpha changes with herring avg wt, not slope

qplot(preyAvgWt, AnnualAlpha, geom="line") 

```

## Seabirds
Common terns were identified at the stakeholder workshop as the recommended seabird herring predator. 

### Tern population model
There is no published stock assessment or population model for most seabirds in the Northeast US. Therefore, Gulf of Maine Common and Arctic tern population parameters were drawn from accounts in the Birds of North America [@hatch_arctic_2002; @nisbet_common_2002] and estimated from counts of breeding pairs and estimates of fledgling success summarized by the Gulf of Maine Seabird Working Group (GOMSWG; data at http://gomswg.org/minutes.html), as corrected and updated by seabird experts from throughout Maine. While we analyzed both Arctic and Common tern information, the stakeholder workshop identified Common terns as the example species for modeling, and this species has more extensive data and a generally higher proportion of herring in its diet based on that data. Therefore, the model is based on common terns in the Gulf of Maine. 

Adult breeding pairs by colony were combined with estimated productivity of fledglings per nest to estimate the annual number of fledglings for each year. A survival rate of 10% was applied to fledglings from each year to represent “recruits” to the breeding adult population age 4 and up [@nisbet_common_2002]. This “stock-recruit” information was used to estimate steepness for the delay difference model based on common tern information only. Fitting parameters with R nls [@Rcite] had variable success, with the full dataset unable to estimate a significant beta parameter (cyan line, Fig. \ref{ternSR}) for common terns, and a truncated dataset resulting in low population production rates inconsistent with currently observed common tern trends (bright green line overlaid with black, Fig. \ref{ternSR}). Therefore, steepness was estimated to give a relationship (black line, Fig. \ref{ternSR}) falling between these two lines. The resulting stock recruit relationship set steepness at 0.26, a theoretical maximum breeding adult population of 45,000 pairs (@nisbet_common_2002, 1930’s New England population), and a theoretical maximum recruitment of 4,500 individuals annually (reflecting approximately a productivity of 1.0 at “carrying capacity” resulting in a stable population). Average common tern productivity is 1.02 (all GOM colony data combined). Adult mortality was assumed to be 0.1 for the delay difference model (survival of 90% [@nisbet_common_2002] for adults).

```{r, fig.cap="Stock-recruitment function for Gulf of Maine common terns \\label{ternSR}", out.width='400pt', message = FALSE, warning = FALSE}
ternpopproddiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMTernPopsProdDiet_Corrected.csv")

#reconcile names
names(ternpopproddiet)[5]<-"productivity"

#get rid of trailing white spaces
ternpopproddiet$majority.diet <- sub(" +$", "", ternpopproddiet$majority.diet)

terntotalpop <- ternpopproddiet %>%
  group_by(yr, species) %>%
  filter(yr<2016) %>%
  select(breedpairs, nfledged) %>%
  filter(!is.na(nfledged)) %>%
  summarise(totpairs = sum(breedpairs, na.rm=T), totfledge = sum(nfledged, na.rm=T))

ternrec <- ggplot(terntotalpop, aes(x=totpairs, y=totfledge*.1, colour=species)) +
  geom_point() +
  scale_y_continuous(limits = c(0,5000)) +
  scale_x_continuous(limits = c(0,45000)) +
  ggtitle("Tern stock and recruitment assuming 10% fledgling->adult survival")

BHfun <- function(x,a,b){
  (x/a)/(1+(1/a*b)*x)
}

BHsteep <- function(Sdat, steep, Rmax, Smax){
  4*steep*Rmax*Sdat/
    (Smax*(1-steep)+Sdat*(5*steep-1))
}

ternrec+
  stat_function(fun=BHfun, args=list(a= 5.445551,b=0.0002266453), col="green", lwd=2)+
  stat_function(fun=BHfun, args=list(a= 9.786020,b=-0.0001570145), col="cyan3")+
  stat_function(fun=BHsteep, args=list(steep=0.41, Rmax=2900, Smax=45000), col="black") +
  stat_function(fun=BHsteep, args=list(steep=0.26, Rmax=4500, Smax=45000), col="black")

```

The resulting model based on common tern population dynamics in the Gulf of Maine (with no link to herring) predicts that the population will increase to its carrying capacity under steady conditions over a 150 year simulation. The actual population has increased at ~2% per year between 1998 and 2015 (Fig. \ref{terntrend}). Given the lack of detailed demographic information in the delay-difference model, this was considered a good representation of the average observed trend in current common tern population dynamics. 

```{r, fig.cap="Population trends for Gulf of Maine terns, no herring link \\label{terntrend}", out.width='400pt', message = FALSE, warning = FALSE}
ternpopproddiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMTernPopsProdDiet_Corrected.csv")

#reconcile names
names(ternpopproddiet)[5]<-"productivity"

#get rid of trailing white spaces
ternpopproddiet$majority.diet <- sub(" +$", "", ternpopproddiet$majority.diet)

terntotalpop <- ternpopproddiet %>%
  group_by(yr, species) %>%
  filter(yr<2016) %>%
  select(breedpairs, nfledged) %>%
  filter(!is.na(nfledged)) %>%
  summarise(totpairs = sum(breedpairs, na.rm=T), totfledge = sum(nfledged, na.rm=T))

terntotalpop19982015 <- terntotalpop %>%
  filter(yr>1997)

ternpoptrend <- ggplot(terntotalpop19982015, aes(x=yr, y=totpairs, colour=species))

# plot with actual data and simulated data from pop model
# WARNING only works with PredN from MSE_TernLink.R reproduced here
#######Some Functions#########
###Function to calculate unfished "SSBR", and abundance or biomass, or equilibrium characteristics given exploitation rate for delay-difference dynamics with Beverton-Holt SR.
#function(unfishedspr,desired unfished N or B, steepness,units of numbs or bio?, annual natural mortality, annual fishery exploitation rate, ford walford parameter intercept, ford walford slope param, do unfished or exploited equilib?)
unfishspr<-function(predunfishspr=NULL,targzero=NULL,steep=NULL,numbs_weightb=NULL,annualnatmort=NULL,annualexploit=NULL,FWalpha=NULL,FWrho=NULL,unfishflag=NULL) {
  #see delay-diff chapter of hilborn and walters and common pubs on steepness (Mangel et al 2010)
  predrzero<-1/predunfishspr*targzero  #unfished recruitment=f(unfishedSSBR,unfished N or B specified by user
  predalpha<-(4*steep*predrzero)/(5*steep-1) #alpha of BH SR #1/beta? alpha in Mangel et al 2010 is (targzero/predrzero)*((1-steep)/(4*steep)) #
  predbeta<-((targzero/predrzero)*((1-steep)/(4*steep)))/((5*steep-1)/(4*steep*predrzero)) #beta of BH SR #alpha/beta? beta ibid (5*steep-1)/(4*steep*predrzero)
  #predalpha<- (targzero/predrzero)*((1-steep)/(4*steep)) #direct from Mangel et al 2010
  #predbeta<-(5*steep-1)/(4*steep*predrzero) #direct from Mangel et al 2010
  if(numbs_weightb==1){ #solving in units of numbers or biomass
    kappa<-1-(1-annualnatmort)*(1-annualexploit) #growth-survival constant (pg 339 Hilborn and Walters)
    eqvalue<-(predalpha/kappa)-(predbeta/(1-annualexploit)) #equilibrium N
    Ceq<-eqvalue*annualexploit #equilibrium catch
    Req<-kappa*eqvalue #equilibrium recruitment
    SSeq<-eqvalue-Ceq  #equilibrium spawning stock
  } else {
    kappa<-((1-(1+FWrho)*(1-annualnatmort)*(1-annualexploit))+(FWrho*(1-annualnatmort)^2*(1-annualexploit)^2))/FWalpha  #growth-survival constant (pg 339 Hilborn and Walters)
    #kappa<-((1-(1+FWrho)*(1-annualnatmort)*(1-annualexploit))+(FWrho*(1-annualnatmort)^2*(1-annualexploit)^2))/(Recwt-((Prerecwt*FWrho)*(1-annualnatmort)*(1-annualexploit)))
    eqvalue<-(predalpha/kappa)-(predbeta/(1-annualexploit)) #equilibrium Biomass
    Yeq<-eqvalue*annualexploit #equilibrium yield
    Req<-kappa*eqvalue #equilibrium recruitment
    SSeq<-eqvalue-Yeq #equilibrium spawning stock
    Neq<-Req/(1-(1-annualnatmort)*(1-annualexploit)) #equilibrium N
    Weq<-eqvalue/Neq #equilibrium average weight
  }  
  if(unfishflag==0){ #if doing unfished then just return unfished spr to optimize function; required for equlibrium calcs
  return((targzero-eqvalue)^2)
  } else { #if not doing unfished (exploit>0) then return various equilibrium values
   if(numbs_weightb==1) { #units of numbers or biomass?
     return(data.frame("Neq"=eqvalue,"Ceq"=Ceq,"Req"=Req,"SSeq"=SSeq,"predrzero"=predrzero,"predalpha"=predalpha,"predbeta"=predbeta))
   }else {
     return(data.frame("Beq"=eqvalue,"Yeq"=Yeq,"Req"=Req,"SSeq"=SSeq,"Neq"=Neq,"Weq"=Weq,"predrzero"=predrzero,"predalpha"=predalpha,"predbeta"=predbeta))
   }
  }
}
#####end function unfishspr#####

#### Common Tern Gulf of Maine pars (see HerringMSEpredspars.xlsx for derivation

numbs_weight  <-  1  #2 #1=numbers, else weight.  Will predator be tracked just in numbers or numbers and weight?  Weight allows for growth effects.
Predzero      <-  45000 #5.90E+05 # #5.60E+09  #50000 #predator unfished total N or biomass
Predsteep     <-  0.26 #0.25  #0.8 #predator steepness
PredannualA   <-  0.1  #0.6 #predator annual natural mortality rate (0-1); serves as max if time varying
Predexploit   <-  0.00  #0.3 #predator annual exploitation rate (0-1)
FWalpha       <-  0.00015#1.5 #  #0.243 #ford-walford plot intercept parameter if doing biomass units tons
FWrho         <-  0.0  #1.45 #ford-walford plot slope parameter if doing biomass units; serves as max if time varying

PredN         <-  c() #for predator abundance vector
PredN[1]      <-  3000  #7000 #initial abundance in year 1
PredB         <-  c() #for predator biomass vector if doing in those units
PredB[1]      <-  1.5  #6000 #initial biomass in year 1 
PredRecdelay  <-  4  #1 #delay in years for when recruits in year y are added to population

preypredrec<-1 #strength of effect of prey N on predator recruitment (like Plaganyi and Butterworth) >=1; 1=no effect
preyprednatm<-0 #strength of effect of prey N on predator annual natural mort (0-1) 0=no effect
preypredgrow<-1 #strengh of effect of prey N on predator growth, >=1, 1=no effect
  
COTEprodThresh<-400000

#####End user inputs########

#Do predator unfished equilibrium calculation and get unfished spr
predunfishsprfxn<-optimize(f=unfishspr,interval=c(0,50),targzero=Predzero,steep=Predsteep,numbs_weightb=numbs_weight,annualnatmort=PredannualA,annualexploit=0,FWalpha=FWalpha,FWrho=FWrho,unfishflag=0,tol=0.0000001)
predunfishspr<-predunfishsprfxn$minimum
#send unfished spr through function again with an exploitation rate to get equlibrium conditions
PredEquilib<-unfishspr(predunfishspr=predunfishspr,targzero=Predzero,steep=Predsteep,numbs_weightb=numbs_weight,annualnatmort=PredannualA,annualexploit=Predexploit,FWalpha=FWalpha,FWrho=FWrho,unfishflag=1)
predalpha<-PredEquilib$predalpha #predator BH SR parm
predbeta<-PredEquilib$predbeta #predator BH SR parm
Req<-PredEquilib$Req
Weq<-PredEquilib$Weq
##############################################
#### Read in Herring "data"
preyB<-read.table("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Unadj1111TotBioSimYear.txt", header=T)  
preysim<-read.table("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Unadj1111NAASimYear.txt", header=T)
preyNsim<-transmute(preysim, preyN=Age1+Age2+Age3+Age4+Age5+Age6+Age7+Age8, Sim=Sim)
     
by_simN<-group_by(preyNsim, Sim)
Nyears<-summarise(by_simN, n=n())

preychar<-read.table("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Unadj1111SimCharAA.txt", header=T)
preychar<-cbind(preychar, Age=rep(1:8,max(preychar$Sim)))
#BASE ON TOTAL B
Ternforage<-cbind(preyB, Yr=rep(1:Nyears$n, max(preychar$Sim)))
h<-1
  preyN<-Ternforage$TotalBio[Ternforage$Sim==h]
  nyears<-Nyears$n[h]
  #preyNzero for terns is the threshold total B where prod drops, FIXED based on GOM data
  preyNzero<-COTEprodThresh

  ###predator dynamics###
  #stuff I need
  Recmult<-(preypredrec*(preyN/preyNzero))/((preypredrec-1)+(preyN/preyNzero)) #fraction of expected predator BH recruitment dependent on prey N
  PredAnnualNatM<-PredannualA*exp(-(preyN/preyNzero)*preyprednatm) #needed if annual nat mort is time varying
  TotalSurv<-(1-PredAnnualNatM)*(1-Predexploit) #total annual predator survival
  catch<-c()
  Spawn<-c() #spawners in numbers or biomass depending
  yield<-c()
  Predrec<-c()
  Predrec[1:PredRecdelay]<-Req #set recruitment in initial years at equlibrium to account for delay
  ###year loop for predator dynamics
  for(y in 1:(nyears-1)){
    if(numbs_weight==1){
        catch[y]<-PredN[y]*Predexploit
        Spawn[y]<-PredN[y]-catch[y]
      } else {
        yield[y]<-PredB[y]*Predexploit
        Spawn[y]<-PredB[y]-yield[y]
      }
      Predrec[y+PredRecdelay]<-Recmult[y]*((predalpha*Spawn[y])/(predbeta+Spawn[y])) #SR
      PredN[y+1]<-PredN[y]*TotalSurv[y]+Predrec[y+1]
      }  

simTern <- data.frame(yr=seq(1998,2015,1), PredN=PredN[6:23])

ternpoptrend +   geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  geom_point(data=simTern, aes(x=yr, y=PredN, colour="Simulated"))


```

### Herring-tern relationship model
The relationship between herring abundance and tern reproductive success was built based on information from individual colonies on annual productivity, proportion of herring in the diet, and amount of herring in the population as estimated by the current stock assessment. Since little of this information has appeared in the peer-reviewed literature, we present it in detail here.  First, productivity information was evaluated by major diet item recorded for chicks over all colonies and years. In general, common tern productivity was higher when a streamlined fish species was the major diet item relative to invertebrates, but having herring as the major diet item resulted in about the same distribution of productivities as having hake or sandlance as the major diet item for these colonies (Fig. \ref{chickdiet}). 

```{r, fig.cap="Major diet items for Gulf of Maine tern fledgelings \\label{chickdiet}",  out.width='400pt', message = FALSE, warning = FALSE}
#knitr::include_graphics("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/TernMajorDietItems.png")
ternpopproddiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMTernPopsProdDiet_Corrected.csv")

#reconcile names
names(ternpopproddiet)[5]<-"productivity"

#get rid of trailing white spaces
ternpopproddiet$majority.diet <- sub(" +$", "", ternpopproddiet$majority.diet)

majorterndiet <- ternpopproddiet %>% 
  filter(!majority.diet %in% c("", "hake and sandlance",
                             "hake and butterfish",
                             "herring and hake",
                             "herring present"))

dietitem <- ggplot(majorterndiet, aes(x=majority.diet, y=productivity, colour=species)) +
  geom_hline(yintercept = 1) + 
  geom_boxplot()

dietitem + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))                  

```

Individual colonies showed different trends in number of nesting pairs, productivity, and proportion of herring in the diet (plots available upon request). When both Arctic and Common terns shared a colony, interannual changes in productivity were generally similar between species, suggesting that conditions at and around the colony (weather, predation pressure, and prey fields) strongly influenced productivity rather than species-specific traits. Only two colonies (Machias Seal Island near the Canadian Border and Stratton Island in southern Maine) showed a significant positive correlation between the proportion of herring in the chick diet and productivity. Other islands showed either non-significant (no) relationships, or in one case (Metinic Island) a significant negative relationship (Fig. \ref{proddiet}). 

```{r, fig.cap="Herring proporiton in diet and tern productivity by colony \\label{proddiet}", out.width='400pt', message = FALSE, warning = FALSE}
#knitr::include_graphics("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/ternproddiet.png")
ternpopproddiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMTernPopsProdDiet_Corrected.csv")

#reconcile names
names(ternpopproddiet)[5]<-"productivity"

#get rid of trailing white spaces
ternpopproddiet$majority.diet <- sub(" +$", "", ternpopproddiet$majority.diet)

ternproddiet <- ggplot(ternpopproddiet, aes(x=prop.herring, y=productivity, colour=species)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Herring proportion in diet and tern productivity GOM")

ternproddiet + facet_wrap(~colony) #+ theme(text=element_text(size=16))

```

The estimated population size of herring on the Northeast US shelf had some relationship to the amount of herring in tern diet at several colonies (4 of 13 common tern colony diets related to herring Age 1 recruitment, 6 of 13 common tern colony diets related to herring total B, and 4 of 13 common tern colony diets related to herring SSB; detailed statistics and plots available upon request). However, statistically significant direct relationships between herring population size and tern productivity were rare, with only Ship Island productivity increasing with herring total B, and Eastern Egg Rock, Matinicus Rock, Ship, and Monomoy Islands productivity increasing with herring SSB. Given that Monomoy Island tern chicks consistently received the lowest proportion of herring in their diets of any colony (0-11%), we don’t consider this relationship further to build the model.  

Based on tern feeding observations, we would expect the number of age 1 herring in the population to be most related to tern productivity since that is the size class terns target, but this relationship was not found in analyzing the data. Herring total biomass was positively related to tern diets at nearly half of the colonies, and reflects all size classes including the smaller sizes most useful as tern forage, but was only directly related to tern productivity at one colony.  Herring SSB was not considered further as an index of tern prey because it represents sizes larger than tern forage. 

<!--Examining all colonies together shows some relationship between herring population size and herring in tern diets, but no consistent relationship between the proportion of herring in tern diets and tern productivity, or between herring population size and tern productivity. Therefore, a reasonable lower bound on the potential for herring to influence tern productivity in the Gulf of Maine is “no relationship,” and any herring control rule (short of one that severely overfishes or eliminates the stock) would not expected to affect tern productivity under this assumption.--> 

To represent the potential for herring to influence tern productivity, we parameterized a tern “recruitment multiplier” based on herring assessed total biomass and common tern productivity across all colonies (except Monomoy where terns eat sandlance). This relationship includes a threshold herring biomass where common tern productivity would drop below 1.0, and above that threshold productivity exceeds 1.0 (Fig. \ref{herrternmod}). The threshold of ~400,000 tons is set where a linear relationship between herring total biomass and common tern productivity crosses productivity=1 (black dashed line in Fig. \ref{herrternmod}). However, the selected threshold is uncertain because there are few observations of common tern productivity at low herring total biomass (1975-1985). The linear relationship does not have a statistically significant slope; a curve was fit to represent a level contribution of herring total biomass to common tern productivity above the threshold. The curve descends below the threshold, dropping below 0.5 productivity at around 50,000 tons and representing the extreme assumption that herring extinction would result in tern productivity of 0. Although the relationship of tern productivity to herring biomass at extremely low herring populations has not been quantified, control rules that allow herring extinction do not meet stated management objectives for herring, so this extreme assumption for terns will not change any decisions to include or exclude control rules. 

```{r, fig.cap="Modeled influence of herring total biomass on tern reproductive success \\label{herrternmod}", out.width='400pt', message = FALSE, warning = FALSE}
ternpopproddiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMTernPopsProdDiet_Corrected.csv")

#reconcile names
names(ternpopproddiet)[5]<-"productivity"

#get rid of trailing white spaces
ternpopproddiet$majority.diet <- sub(" +$", "", ternpopproddiet$majority.diet)

herring <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/herringassessout2015.csv")
herring <- herring[,1:6]

ternpopproddiet_herring <- left_join(ternpopproddiet, herring)

terntotalpop <- ternpopproddiet %>%
  group_by(yr, species) %>%
  filter(yr<2016) %>%
  select(breedpairs, nfledged) %>%
  filter(!is.na(nfledged)) %>%
  summarise(totpairs = sum(breedpairs, na.rm=T), totfledge = sum(nfledged, na.rm=T))

# Monomoy Island is a third of the population but eats mainly sandlance and is south of GOM proper--leave out?

terntotpop_noMon <- ternpopproddiet %>%
  group_by(yr, species) %>%
  filter(colony != "Monomoy") %>%
  filter(yr<2016) %>%
  select(breedpairs, nfledged) %>%
  filter(!is.na(nfledged)) %>%
  summarise(totpairs = sum(breedpairs, na.rm=T), totfledge = sum(nfledged, na.rm=T))

ternpop_herringpop <- left_join(terntotalpop, herring)
ternpopnomon_herringpop <- left_join(terntotpop_noMon, herring)

totprodnomonherr <- ggplot(ternpopnomon_herringpop, aes(x=HerringTotalB, y=totfledge/totpairs, colour=species)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  ggtitle("Tern production (no Monomoy) and assessed herring tot B")

#totprodnomonherr + stat_smooth(method = "lm", formula = y ~ x, size = 1) 

HerrAlpha <- 1.09
HerringThresh <- 400000

recmult<-(HerrAlpha*(ternpopnomon_herringpop$HerringTotalB/HerringThresh))/((HerrAlpha-1)+(ternpopnomon_herringpop$HerringTotalB/HerringThresh))

#totprodnomonherr + geom_line(aes(x=HerringTotalB, y=recmult), lwd=1.5)

lmtotnomon <- lm(totfledge/totpairs ~ HerringTotalB, species=="Common", data=ternpopnomon_herringpop)

#totprodnomonherr + geom_line(aes(x=HerringTotalB, y=recmult), lwd=1.5) +
  #geom_line(aes(x=HerringTotalB, y=9.439e-01 + 1.210e-07*HerringTotalB), color="black", lty=3)

totprodnomonherr + 
  scale_y_continuous(limits = c(0,1.8)) +
  scale_x_continuous(limits = c(0,2500000)) +
  geom_line(aes(x=HerringTotalB, y=recmult), lwd=1.5) +
  geom_line(aes(x=HerringTotalB, y=9.439e-01 + 1.210e-07*HerringTotalB), color="black", lty=3) +
  geom_line(aes(x=seq(0,2500000,length.out=51), y=(HerrAlpha*(seq(0,2500000,length.out=51)/HerringThresh))/((HerrAlpha-1)+(seq(0,2500000,length.out=51)/HerringThresh))))


```

When included in the model using $\gamma$ = 1.09 in equation \ref{recwithherring_equation}, this relationship adjusts the modeled common tern population increase to match the current average increase in common tern nesting pairs observed in the data (Fig. \ref{terntrendwherring}). There is still considerable uncertainty around this mean population trajectory which cannot be reflected in our simple model. 

```{r, fig.cap="Population trends for Gulf of Maine terns with simulated herring-common tern productivity relationship \\label{terntrendwherring}", out.width='400pt', message = FALSE, warning = FALSE}
ternpopproddiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMTernPopsProdDiet_Corrected.csv")

#reconcile names
names(ternpopproddiet)[5]<-"productivity"

#get rid of trailing white spaces
ternpopproddiet$majority.diet <- sub(" +$", "", ternpopproddiet$majority.diet)

terntotalpop <- ternpopproddiet %>%
  group_by(yr, species) %>%
  filter(yr<2016) %>%
  select(breedpairs, nfledged) %>%
  filter(!is.na(nfledged)) %>%
  summarise(totpairs = sum(breedpairs, na.rm=T), totfledge = sum(nfledged, na.rm=T))

terntotalpop19982015 <- terntotalpop %>%
  filter(yr>1997)

ternpoptrend <- ggplot(terntotalpop19982015, aes(x=yr, y=totpairs, colour=species))

# plot with actual data and simulated data from pop model
# WARNING only works with PredN from MSE_TernLink.R reproduced here
#######Some Functions#########
###Function to calculate unfished "SSBR", and abundance or biomass, or equilibrium characteristics given exploitation rate for delay-difference dynamics with Beverton-Holt SR.
#function(unfishedspr,desired unfished N or B, steepness,units of numbs or bio?, annual natural mortality, annual fishery exploitation rate, ford walford parameter intercept, ford walford slope param, do unfished or exploited equilib?)
unfishspr<-function(predunfishspr=NULL,targzero=NULL,steep=NULL,numbs_weightb=NULL,annualnatmort=NULL,annualexploit=NULL,FWalpha=NULL,FWrho=NULL,unfishflag=NULL) {
  #see delay-diff chapter of hilborn and walters and common pubs on steepness (Mangel et al 2010)
  predrzero<-1/predunfishspr*targzero  #unfished recruitment=f(unfishedSSBR,unfished N or B specified by user
  predalpha<-(4*steep*predrzero)/(5*steep-1) #alpha of BH SR #1/beta? alpha in Mangel et al 2010 is (targzero/predrzero)*((1-steep)/(4*steep)) #
  predbeta<-((targzero/predrzero)*((1-steep)/(4*steep)))/((5*steep-1)/(4*steep*predrzero)) #beta of BH SR #alpha/beta? beta ibid (5*steep-1)/(4*steep*predrzero)
  #predalpha<- (targzero/predrzero)*((1-steep)/(4*steep)) #direct from Mangel et al 2010
  #predbeta<-(5*steep-1)/(4*steep*predrzero) #direct from Mangel et al 2010
  if(numbs_weightb==1){ #solving in units of numbers or biomass
    kappa<-1-(1-annualnatmort)*(1-annualexploit) #growth-survival constant (pg 339 Hilborn and Walters)
    eqvalue<-(predalpha/kappa)-(predbeta/(1-annualexploit)) #equilibrium N
    Ceq<-eqvalue*annualexploit #equilibrium catch
    Req<-kappa*eqvalue #equilibrium recruitment
    SSeq<-eqvalue-Ceq  #equilibrium spawning stock
  } else {
    kappa<-((1-(1+FWrho)*(1-annualnatmort)*(1-annualexploit))+(FWrho*(1-annualnatmort)^2*(1-annualexploit)^2))/FWalpha  #growth-survival constant (pg 339 Hilborn and Walters)
    #kappa<-((1-(1+FWrho)*(1-annualnatmort)*(1-annualexploit))+(FWrho*(1-annualnatmort)^2*(1-annualexploit)^2))/(Recwt-((Prerecwt*FWrho)*(1-annualnatmort)*(1-annualexploit)))
    eqvalue<-(predalpha/kappa)-(predbeta/(1-annualexploit)) #equilibrium Biomass
    Yeq<-eqvalue*annualexploit #equilibrium yield
    Req<-kappa*eqvalue #equilibrium recruitment
    SSeq<-eqvalue-Yeq #equilibrium spawning stock
    Neq<-Req/(1-(1-annualnatmort)*(1-annualexploit)) #equilibrium N
    Weq<-eqvalue/Neq #equilibrium average weight
  }  
  if(unfishflag==0){ #if doing unfished then just return unfished spr to optimize function; required for equlibrium calcs
  return((targzero-eqvalue)^2)
  } else { #if not doing unfished (exploit>0) then return various equilibrium values
   if(numbs_weightb==1) { #units of numbers or biomass?
     return(data.frame("Neq"=eqvalue,"Ceq"=Ceq,"Req"=Req,"SSeq"=SSeq,"predrzero"=predrzero,"predalpha"=predalpha,"predbeta"=predbeta))
   }else {
     return(data.frame("Beq"=eqvalue,"Yeq"=Yeq,"Req"=Req,"SSeq"=SSeq,"Neq"=Neq,"Weq"=Weq,"predrzero"=predrzero,"predalpha"=predalpha,"predbeta"=predbeta))
   }
  }
}
#####end function unfishspr#####

#### Common Tern Gulf of Maine pars (see HerringMSEpredspars.xlsx for derivation

numbs_weight  <-  1  #2 #1=numbers, else weight.  Will predator be tracked just in numbers or numbers and weight?  Weight allows for growth effects.
Predzero      <-  45000 #5.90E+05 # #5.60E+09  #50000 #predator unfished total N or biomass
Predsteep     <-  0.26 #0.25  #0.8 #predator steepness
PredannualA   <-  0.1  #0.6 #predator annual natural mortality rate (0-1); serves as max if time varying
Predexploit   <-  0.00  #0.3 #predator annual exploitation rate (0-1)
FWalpha       <-  0.00015#1.5 #  #0.243 #ford-walford plot intercept parameter if doing biomass units tons
FWrho         <-  0.0  #1.45 #ford-walford plot slope parameter if doing biomass units; serves as max if time varying

PredN         <-  c() #for predator abundance vector
PredN[1]      <-  3000  #7000 #initial abundance in year 1
PredB         <-  c() #for predator biomass vector if doing in those units
PredB[1]      <-  1.5  #6000 #initial biomass in year 1 
PredRecdelay  <-  4  #1 #delay in years for when recruits in year y are added to population

preypredrec<-1.09 #strength of effect of prey N on predator recruitment (like Plaganyi and Butterworth) >=1; 1=no effect
preyprednatm<-0 #strength of effect of prey N on predator annual natural mort (0-1) 0=no effect
preypredgrow<-1 #strengh of effect of prey N on predator growth, >=1, 1=no effect
  
COTEprodThresh<-400000

#####End user inputs########

#Do predator unfished equilibrium calculation and get unfished spr
predunfishsprfxn<-optimize(f=unfishspr,interval=c(0,50),targzero=Predzero,steep=Predsteep,numbs_weightb=numbs_weight,annualnatmort=PredannualA,annualexploit=0,FWalpha=FWalpha,FWrho=FWrho,unfishflag=0,tol=0.0000001)
predunfishspr<-predunfishsprfxn$minimum
#send unfished spr through function again with an exploitation rate to get equlibrium conditions
PredEquilib<-unfishspr(predunfishspr=predunfishspr,targzero=Predzero,steep=Predsteep,numbs_weightb=numbs_weight,annualnatmort=PredannualA,annualexploit=Predexploit,FWalpha=FWalpha,FWrho=FWrho,unfishflag=1)
predalpha<-PredEquilib$predalpha #predator BH SR parm
predbeta<-PredEquilib$predbeta #predator BH SR parm
Req<-PredEquilib$Req
Weq<-PredEquilib$Weq
##############################################
#### Read in Herring "data"
preyB<-read.table("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Unadj1111TotBioSimYear.txt", header=T)  
preysim<-read.table("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Unadj1111NAASimYear.txt", header=T)
preyNsim<-transmute(preysim, preyN=Age1+Age2+Age3+Age4+Age5+Age6+Age7+Age8, Sim=Sim)
     
by_simN<-group_by(preyNsim, Sim)
Nyears<-summarise(by_simN, n=n())

preychar<-read.table("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Unadj1111SimCharAA.txt", header=T)
preychar<-cbind(preychar, Age=rep(1:8,max(preychar$Sim)))
#BASE ON TOTAL B
Ternforage<-cbind(preyB, Yr=rep(1:Nyears$n, max(preychar$Sim)))
h<-1
  preyN<-Ternforage$TotalBio[Ternforage$Sim==h]
  nyears<-Nyears$n[h]
  #preyNzero for terns is the threshold total B where prod drops, FIXED based on GOM data
  preyNzero<-COTEprodThresh

  ###predator dynamics###
  #stuff I need
  Recmult<-(preypredrec*(preyN/preyNzero))/((preypredrec-1)+(preyN/preyNzero)) #fraction of expected predator BH recruitment dependent on prey N
  PredAnnualNatM<-PredannualA*exp(-(preyN/preyNzero)*preyprednatm) #needed if annual nat mort is time varying
  TotalSurv<-(1-PredAnnualNatM)*(1-Predexploit) #total annual predator survival
  catch<-c()
  Spawn<-c() #spawners in numbers or biomass depending
  yield<-c()
  Predrec<-c()
  Predrec[1:PredRecdelay]<-Req #set recruitment in initial years at equlibrium to account for delay
  ###year loop for predator dynamics
  for(y in 1:(nyears-1)){
    if(numbs_weight==1){
        catch[y]<-PredN[y]*Predexploit
        Spawn[y]<-PredN[y]-catch[y]
      } else {
        yield[y]<-PredB[y]*Predexploit
        Spawn[y]<-PredB[y]-yield[y]
      }
      Predrec[y+PredRecdelay]<-Recmult[y]*((predalpha*Spawn[y])/(predbeta+Spawn[y])) #SR
      PredN[y+1]<-PredN[y]*TotalSurv[y]+Predrec[y+1]
      }  

simTern <- data.frame(yr=seq(1998,2015,1), PredN=PredN[6:23])

ternpoptrend +   geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  geom_point(data=simTern, aes(x=yr, y=PredN, colour="Simulated"))


```

## Groundfish
Because no specific groundfish was identified as a representative herring predator during the stakeholder workshop, the first decision was which groundfish to model. Annual diet estimates (based on sample sizes of ~100+ stomachs) are available for the top three groundfish predators of herring (those with herring occurring in the diets most often in the entire NEFSC food habits database): spiny dogfish, Atlantic cod, and silver hake. Cod and spiny dogfish were considered first because their overall diet proportions of herring are higher, and because silver hake has the least recently updated assessment. Diet compositions by year were estimated for spiny dogfish, Georges Bank cod, and Gulf of Maine cod to match the scale of stock assessments. Full weighted diet compositions were estimated, and suggest considerable interannual variability in the herring proportion in groundfish diets (filled blue proportions of bars in Fig. \ref{gfishdiets}).  

```{r,fig.cap="Annual diet compositions for major groundfish predators of herring estimated from NEFSC food habits database \\label{gfishdiets}", message = FALSE, warning = FALSE}

dogdiet   <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/GroundfishDat/dogfishdietcomp.csv")
GBcoddiet <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/GroundfishDat/GBcoddietcomp.csv")
GOMcoddiet <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/GroundfishDat/GOMcoddietcomp.csv")

names(dogdiet)[1]<-"Year"
names(GBcoddiet)[1]<-"Year" 
names(GOMcoddiet)[1]<-"Year"

dogdiet <- cbind(dogdiet, Pred=rep("dogfish", dim(dogdiet)[1]))
GBcoddiet <- cbind(GBcoddiet, Pred=rep("GBcod", dim(GBcoddiet)[1]))
GOMcoddiet <- cbind(GOMcoddiet, Pred=rep("GOMcod", dim(GOMcoddiet)[1]))

dogdiet <- gather(dogdiet, Prey, Percent, -Year, -Pred)
GBcoddiet <- gather(GBcoddiet, Prey, Percent, -Year, -Pred)
GOMcoddiet <- gather(GOMcoddiet, Prey, Percent, -Year, -Pred)

dogcoddiet <- bind_rows(dogdiet, GBcoddiet, GOMcoddiet)

dogcoddiet$fillwhite <- with(dogcoddiet, ifelse(Prey=="CLUHAR", 1,0))

compplot2 <- ggplot(dogcoddiet, aes(Year, Percent, colour=Prey)) + 
  geom_bar(stat = "identity", aes(fill=fillwhite)) 
  #geom_bar(stat = "identity") 

compplot2 + facet_wrap("Pred", nrow=3) + theme(legend.position="none")

```

Some interannual variation in diet may be explained by changing herring abundance. Dogfish and both cod stocks had positive relationships between the amount of herring observed in annual diets and the size of the herring population according to the most recent assessment (statistics and plots available upon request). This suggests that these groundfish predators are opportunistic, eating herring in proportion to their availability in the ecosystem. However, monotonically declining cod populations for both GOM and GB cod stocks resulted in either no herring-cod relationship, or a negative relationship between herring populations and cod populations (Fig. \ref{gfishBherrdiet}). Only dogfish spawning stock biomass had a positive relationship with the proporiton of herring in dogfish diet. Therefore, we selected dogfish as the groundfish predator for modeling.

```{r, fig.cap="Relationship of assessed groundfish spawning stock biomass (SSB) with the proportion of herring in diet \\label{gfishBherrdiet}", out.width='400pt',message = FALSE, warning = FALSE}
codGB<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GBcodpop.csv")
codGOM<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMcodpop.csv")
dogfish<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/dogfishpop.csv")
gfishdiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Gfishdiets.csv")
herring <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/herringassessout2015.csv")
herring <- herring[,1:6]
dogfishrec <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/dogfishrec.csv")

# add species column, tidy and combine datasets
codGB <- cbind(codGB, species=rep("GBcod", dim(codGB)[1]))
codGOM <- cbind(codGOM, species=rep("GOMcod", dim(codGOM)[1]))
dogfish <- cbind(dogfish, species=rep("dogfish", dim(dogfish)[1]))
dogfish <- full_join(dogfish, dogfishrec)

gfishdiet <- gfishdiet %>%
  gather(preypred, dietprop, 2:9) %>%
  separate(preypred, c("prey", "species"), sep="_")

codpops <- full_join(codGB, codGOM)
names(codpops)[1]<-"Yr"
names(dogfish)[2:12] <- c("LargeB", "MedB", "RecAge1t",
                          "TotJan1B","survSSB", "survSSB_MA","stocExpB",
                          "stocTotB", "stocSSB", "SSBkgtow", "SSB")
#put all in t not 1000t
dogfish <- dogfish %>%
  mutate_at(vars(LargeB, MedB, RecAge1t,TotJan1B, survSSB, survSSB_MA,
                 stocExpB, stocTotB, stocSSB, SSB),
            funs(.*1000))

dogfish <- dogfish %>%
  mutate(RecAge1 = RecAge1t/PupAvgWtkg) %>% #rec in 1000 fish matches cod units
  mutate(TotJan1B = replace(TotJan1B, TotJan1B==0, NA)) #2014 missed survey

gfishpops <- full_join(codpops, dogfish)

#merge fishpops and fish diets
gfishpopfood <- left_join(gfishpops, gfishdiet, by=c("Yr"="Year", "species"="species"))

gfishpopfood <- filter(gfishpopfood, !is.na(prey)) 

gfishSSBdiet <- ggplot(gfishpopfood, aes(x=dietprop, y=SSB, colour=prey)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) 
gfishSSBdiet + facet_wrap(~species)


```

### Dogfish population model
The dogfish model stock recruitment function, initial population, and annual natural mortality were adapted from information in [@rago_implications_1998; @rago_biological_2010; @bubley_reassessment_2012; @rago_update_2013]. Due to differential growth and fishing mortality by sex, our model best represents female dogfish (a split-sex delay difference model was not feasible within the time constraints of this MSE). Further, dogfish stock-recruit modeling to date based on Ricker functions [@rago_biological_2010] captures more nuances in productivity than the Beverton-Holt model we used. Our recruitment parameterization reflects a stock with generally low productivity and relatively high resilience, which we recognize is a rough approximation for a species such as dogfish. The annual fishing exploitation rate applied is average of the catch/adult female biomass from the most recent years of the 2016 data update provided to the Mid-Atlantic Fishery Management Council (Rago pers comm). 

### Herring-dogfish relationship model
There was a weak positive relationship between dogfish total biomass and herring total biomass from the respective stock assessments (Fig. \ref{pupherring}), but no clear relationship between dogfish weight or dogfish recruitment and herring population size. During the recent period of relatively low dogfish recruitment (1995-2007), there is a positive relationship between dogfish pup average weight and herring proportion in diet, suggesting a potential growth and or recruitment mechanism; however this relationship does not hold throughout the time series (Fig. \ref{pupherring}). 

```{r, fig.cap="Dogfish population relationships with herring total biomass (left) and herring proportion in diet (right) \\label{pupherring}", out.width='.49\\linewidth', fig.show='hold', fig.align='center', message = FALSE, warning = FALSE}
codGB<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GBcodpop.csv")
codGOM<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/GOMcodpop.csv")
dogfish<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/dogfishpop.csv")
gfishdiet<-read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/Gfishdiets.csv")
herring <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/herringassessout2015.csv")
herring <- herring[,1:6]
dogfishrec <- read.csv("/Users/sgaichas/Documents/0_Data/MSE/HerringMSE/dogfishrec.csv")

# add species column, tidy and combine datasets
codGB <- cbind(codGB, species=rep("GBcod", dim(codGB)[1]))
codGOM <- cbind(codGOM, species=rep("GOMcod", dim(codGOM)[1]))
dogfish <- cbind(dogfish, species=rep("dogfish", dim(dogfish)[1]))
dogfish <- full_join(dogfish, dogfishrec)

gfishdiet <- gfishdiet %>%
  gather(preypred, dietprop, 2:9) %>%
  separate(preypred, c("prey", "species"), sep="_")

codpops <- full_join(codGB, codGOM)
names(codpops)[1]<-"Yr"
names(dogfish)[2:12] <- c("LargeB", "MedB", "RecAge1t",
                          "TotJan1B","survSSB", "survSSB_MA","stocExpB",
                          "stocTotB", "stocSSB", "SSBkgtow", "SSB")
#put all in t not 1000t
dogfish <- dogfish %>%
  mutate_at(vars(LargeB, MedB, RecAge1t,TotJan1B, survSSB, survSSB_MA,
                 stocExpB, stocTotB, stocSSB, SSB),
            funs(.*1000))

dogfish <- dogfish %>%
  mutate(RecAge1 = RecAge1t/PupAvgWtkg) %>% #rec in 1000 fish matches cod units
  mutate(TotJan1B = replace(TotJan1B, TotJan1B==0, NA)) #2014 missed survey

gfishpops <- full_join(codpops, dogfish)

#merge fishpops and fish diets
gfishpopfood <- left_join(gfishpops, gfishdiet, by=c("Yr"="Year", "species"="species"))
gfishpopfood <- filter(gfishpopfood, !is.na(prey)) 
# merge herring assessment out with gfish
gfishpopfoodherr <- left_join(gfishpopfood, herring, by=c("Yr"="yr"))

#Dogfish population model
dogfishfoodherr <- gfishpopfoodherr %>%
  filter(species=="dogfish") 

#dogfishpoptrend <- ggplot(dogfishfoodherr, aes(x=Yr, y=TotJan1B)) + geom_point()
#dogfishdiet <- ggplot(dogfishfoodherr, aes(x=Yr, y=dietprop, colour=prey)) + geom_point()
pupwtdiet <- ggplot(dogfishfoodherr, aes(x=dietprop, y=PupAvgWtkg, colour=prey)) + geom_point() + ggtitle("Dogfish pup average weights 1980-2009")
#Fwtdiet <- ggplot(dogfishfoodherr, aes(x=dietprop, y=FAvgWtkg, colour=prey)) + geom_point()

dogfishlow <- dogfishfoodherr %>%
  filter(Yr>1994 & Yr<2008)

#lowdogfishdiet <- ggplot(dogfishlow, aes(x=Yr, y=dietprop, colour=prey)) + geom_point()
lowpupwtdiet <- ggplot(dogfishlow, aes(x=dietprop, y=PupAvgWtkg, colour=prey)) + geom_point() + ggtitle("Dogfish pup average weights 1995-2007")
#lowFwtdiet <- ggplot(dogfishlow, aes(x=dietprop, y=FAvgWtkg, colour=prey)) + geom_point()

#lowpupwtdiet #+ stat_smooth(method = "lm", formula = y ~ x, size = 1)
#pupwtdiet #+ stat_smooth(method = "lm", formula = y ~ x, size = 1)

#pup wt was higher earlier in the time series when herring low
#pupwttrend <- ggplot(dogfishfoodherr, aes(x=Yr, y=PupAvgWtkg)) + geom_point() + ggtitle("Dogfish average pup weight from NEFSC surveys")

#positive relationship between populations but weak
dogherrB <- ggplot(dogfishfoodherr, aes(x=HerringTotalB, y=TotJan1B)) +
  geom_point() +
  ggtitle("Dogfish total biomass 1968-2016")
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) + 
  #stat_smooth(method = "lm", formula = y ~ x, size = 1)

dogherrB
pupwtdiet
```

Therefore, to simulate a potential positive relationship between herring and dogfish, we assumed that dogfish survival increased (natural mortality was reduced) by an unspecified mechanism as herring abundance increased (Fig. \ref{herringdogfish}). Because dogfish are fully exploited by fisheries in this model, the impact of this change in natural mortality on total survival has small to moderate benefits to dogfish population numbers and biomass. Using a $\delta$ = 0.2 in equation \ref{varmort_equation} results in weak increases in dogfish biomass with herring abundance consistent with observations.

```{r, fig.cap="Modeled herring relative population size-dogfish natural mortality relationship \\label{herringdogfish}", out.width='300pt', message = FALSE, warning = FALSE}
# from MSE_GroundfishLink.R
PredannualA   <-  0.092  #0.6 #predator annual natural mortality rate (0-1); serves as max if time varying
Predexploit   <-  0.092  #0.3 #predator annual exploitation rate (0-1)
preyprednatm<-0.2 #strength of effect of prey N on predator annual natural mort (0-1) 0=no effect

preypopratio  <-  seq(0,2,by=0.01) # standing in for preyN/preyNzero
PredAnnualNatM   <-  c()

PredAnnualNatM<-PredannualA*exp(-(preypopratio)*preyprednatm) 

qplot(preypopratio, PredAnnualNatM, geom="line", xlab = "current herring abundance / unfished herring abundance", ylab = "dogfish annual natural mortality v") 

#TotalSurv<-(1-PredAnnualNatM)*(1-Predexploit) #total annual predator survival

```

## Marine mammals
Because no specific marine mammal was identified as a representative herring predator in the stakeholder workshop, as with groundfish, the first decision was which marine mammal to model. Diet information for a wide range of marine mammals on the Northeast US shelf suggests that minke whales, humpback whales, harbor seals, and harbor porpoises have the highest proportions of herring in their diets [@smith_consumption_2015], and therefore may show some reaction to changes in the herring ABC control rule. 

While some food habits data existed for marine mammals, consultation with marine mammal stock assessment scientists at the Northeast Fisheries Science Center confirmed that no data were available to parameterize a stock-recruitment relationship for any of these marine mammal species in the Northeast US region, and no such information was available in the literature for stocks in this region. Althouth it may be possible to develop stock-recruitment models for one or more of these species in the future, it was not possible within the time frame of the herring MSE. Therefore, we were unable to model marine mammals within the same framework as other predators. 

Potential effects of changes in herring production and/or biomass on marine mammals were instead evaluated using an updated version of an existing food web model for the Gulf of Maine [@link_northeast_2008; @link_response_2009; @link_documentation_2006] and incorporating food web model parameter uncertainty. Overall, food web modeling showed that a simulated increase in herring production in the Gulf of Maine may produce modest but uncertain benefits to marine mammal predators, primarily because increased herring was associated with decreases in other forage groups also preyed on by marine mammals. Please see Appendix 1 of this document for full analyses and results. 

## Summary of predator model input parameters
Table 2. Predator model input parameters
```{r pars, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'}
tabl <- "
| Parameter                | Tuna         | Tern         | Dogfish      |
|:-------------------------|:-------------|:-------------|:-------------|
| Numbers or Weight?       | Weight       | Numbers      | Weight       |
| Unfished spawning pop    | 6.69E+04     | 45000        | 300000       |
| Steepness *h*            | 1.0          | 0.26         | 0.97         |
| Annual nat. mortality *v*| 0.14         | 0.1          | 0.092        |
| Annual exploitation *u*  | 0.079        | 0.00         | 0.092        |
| Growth intercept *FWint* | 0.020605     | 0.00015      | 0.000278     |
| Growth slope *FWslope *  | 0.9675       | 0.0          | 0.9577       |
| Initial abundance        | 111864       | 3000         | 49629630     |
| Initial biomass          | 27966        | 1.5          | 134000       |
| Recruit delay (age) *a*  | 1            | 4            | 10           |
| Prey-recruitment link    | 1 (off)      | 1.09         | 1 (off)      |
| Prey-mortality link      | 0 (off)      | 0 (off)      | 0.2          |
| Prey-growth link         | 1.1          | 1 (off)      | 1 (off)      |
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion
```

