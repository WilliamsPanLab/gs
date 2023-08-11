Figure3
================
2023-05-06

``` r
# figure 3 - # be sure to run this one sequentially! relies on counters with same variable name across chunks.
library(mgcv)
```

    ## Loading required package: nlme

    ## This is mgcv 1.8-42. For overview type 'help("mgcv-package")'.

``` r
library(visreg)
library(gratia)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     collapse

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(stringr)
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
library(ggridges)
library(ggbeeswarm)
```

    ## Loading required package: ggplot2

``` r
library(ggExtra)
knitr::opts_chunk$set(fig.width=6, fig.height=2) 
```

``` r
# set functions
plot_bootstraps <- function(data,maxval,Name,maxValuePlot,BorderlineClinical,Clinical) {
  # Melt the data frame
  data_melt <- melt(t(data))
  data_melt$Var1 <- rep(seq(1, maxval), nrow(data))

  # Calculate percentiles
  percentiles <- data %>%
    summarise(across(everything(), quantile, probs = c(0.01, 0.99), na.rm = TRUE))
  
  percentiles_long <- tidyr::pivot_longer(percentiles, cols = everything(), names_to = "Percentile", values_to = "YValue")

  # Add CI column
  data_melt$CI <- 0
  
  # Prepare CIs for insertion
  CIs <- data.frame(rep(seq(1, maxval), 2), c(rep(10001, maxval), rep(10002, maxval)), percentiles_long$YValue, rep(1, (maxval*2)))
  colnames(CIs) <- colnames(data_melt)
  
  # Add CIs
  data_melt2 <- rbind(data_melt, CIs)
  
  # Convert CI column to factor
  data_melt2$CI <- as.factor(data_melt2$CI)
  
  # Plotting the lines
  ggplot(data = data_melt2, aes(x = Var1, y = value, group = Var2, color = Var2)) +
    geom_line(aes(alpha = CI), show.legend = FALSE) +
    scale_color_viridis_c(option = "inferno", direction = -1) +
    scale_alpha_manual(values = c(0.01, 1), guide = FALSE) + ylim(c(-1.5,1.5)) +
    theme_minimal(base_size=35) + 
    ylab(expression(italic(g)))+xlab(Name)+
    geom_vline(xintercept = BorderlineClinical, linetype = "dashed")+
    geom_vline(xintercept = Clinical, linetype = "dashed")+
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,maxValuePlot),expand = expansion(mult = c(0, 0)))
}

find_furthest_nonzero <- function(data) {
  numZeros=colSums(data==0)
  isZeroZeros=numZeros==0
  furthest_nonzero=sum(isZeroZeros)
}

# set colors
my_palette <- colorRampPalette(colors = c("#051099", "#1d5cb7", "white", "#e41a1c", "#a80009"))
```

``` r
# load in fits: ordered as follows:
# F_intFit,M_intFit,P_intFit,R_intFit,F_extFit,M_extFit,P_extFit,R_extFit,F_somFit,M_somFit,P_somFit,R_somFit,F_anxFit,M_anxFit,P_anxFit,R_anxFit,F_thoFit,M_thoFit,P_thoFit,R_thoFit,F_witdepFit,M_witdepFit,P_witdepFit,R_witdepFit,F_socFit,M_socFit,P_socFit,R_socFit,F_rulFit,M_rulFit,P_rulFit,R_rulFit,F_attFit,M_attFit,P_attFit,R_attFit,F_aggFit,M_aggFit,P_aggFit,R_aggFit
Fits=readRDS('~/Desktop/g_p/F3_gpSubscaleFits.rds')

# read in masterdf from sample construction
masterdf=readRDS('~/gp_masterdf.rds')

# pull clinical cutoff from master df: t scores > 65 = borderline clinical, 70 = clinical
masterdfP_bc<-masterdf[masterdf$cbcl_scr_syn_totprob_t==65,]
masterdfP_c<-masterdf[masterdf$cbcl_scr_syn_totprob_t==69,]

# borderline clinical and clinical cutoffs
Pbc=mean(masterdfP_bc$cbcl_scr_syn_totprob_r)
Pc=mean(masterdfP_c$cbcl_scr_syn_totprob_r)
```

``` r
# use max value derivation from slurm script to parse individual subscales
intMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_internal_r))
extMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_external_r))
somMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_somatic_r))
anxMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_anxdep_r))
thoMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_thought_r))
depMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_withdep_r))
socMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_social_r))
attMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_attention_r))
rulMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_rulebreak_r))
aggMaxVal=max(as.numeric(masterdf$cbcl_scr_syn_aggressive_r))

# and use clinical thresholds from regular t-scores
# pull clinical cutoff from master df: t scores > 65 = borderline clinical, 69 = clinical
masterdfP_bc<-masterdf[masterdf$cbcl_scr_syn_totprob_t==65,]
masterdfP_c<-masterdf[masterdf$cbcl_scr_syn_totprob_t==69,]
masterdfI_bc<-masterdf[masterdf$cbcl_scr_syn_internal_t==65,]
masterdfI_c<-masterdf[masterdf$cbcl_scr_syn_internal_t==69,]
masterdfE_bc<-masterdf[masterdf$cbcl_scr_syn_external_t==65,]
masterdfE_c<-masterdf[masterdf$cbcl_scr_syn_external_t==69,]
masterdfAnx_bc<-masterdf[masterdf$cbcl_scr_syn_anxdep_t==65,]
masterdfAnx_c<-masterdf[masterdf$cbcl_scr_syn_anxdep_t==69,]
# note no one has t==65 in this dataset for thought
masterdfTho_bc<-masterdf[masterdf$cbcl_scr_syn_thought_t==66,]
masterdfTho_c<-masterdf[masterdf$cbcl_scr_syn_thought_t==69,]
# note no one has t==65 in this dataset for withdrawn depression
masterdfWit_bc<-masterdf[masterdf$cbcl_scr_syn_withdep_t==66,]
masterdfWit_c<-masterdf[masterdf$cbcl_scr_syn_withdep_t==69,]
masterdfSom_bc<-masterdf[masterdf$cbcl_scr_syn_somatic_t==65,]
# no one has t=69
masterdfSom_c<-masterdf[masterdf$cbcl_scr_syn_somatic_t==70,]
masterdfSoc_bc<-masterdf[masterdf$cbcl_scr_syn_social_t==65,]
masterdfSoc_c<-masterdf[masterdf$cbcl_scr_syn_social_t==69,]
masterdfAtt_bc<-masterdf[masterdf$cbcl_scr_syn_attention_t==65,]
masterdfAtt_c<-masterdf[masterdf$cbcl_scr_syn_attention_t==69,]
masterdfRul_bc<-masterdf[masterdf$cbcl_scr_syn_rulebreak_t==65,]
masterdfRul_c<-masterdf[masterdf$cbcl_scr_syn_rulebreak_t==69,]
masterdfAgg_bc<-masterdf[masterdf$cbcl_scr_syn_aggressive_t==65,]
masterdfAgg_c<-masterdf[masterdf$cbcl_scr_syn_aggressive_t==69,]

# borderline clinical and clinical cutoffs
Pbc=mean(masterdfP_bc$cbcl_scr_syn_totprob_r)
Pc=mean(masterdfP_c$cbcl_scr_syn_totprob_r)
Ibc=mean(masterdfP_bc$cbcl_scr_syn_internal_r)
Ic=mean(masterdfP_c$cbcl_scr_syn_internal_r)
Ebc=mean(masterdfE_bc$cbcl_scr_syn_external_r)
Ec=mean(masterdfE_c$cbcl_scr_syn_external_r)
AnxBc=mean(as.numeric(masterdfAnx_bc$cbcl_scr_syn_anxdep_r))
AnxC=mean(as.numeric(masterdfAnx_c$cbcl_scr_syn_anxdep_r))
ThoBc=mean(as.numeric(masterdfTho_bc$cbcl_scr_syn_thought_r))
ThoC=mean(as.numeric(masterdfTho_c$cbcl_scr_syn_thought_r))
WitBc=mean(as.numeric(masterdfWit_bc$cbcl_scr_syn_withdep_r))
WitC=mean(as.numeric(masterdfWit_c$cbcl_scr_syn_withdep_r))
SomBc=mean(as.numeric(masterdfSom_bc$cbcl_scr_syn_somatic_r))
SomC=mean(as.numeric(masterdfSom_c$cbcl_scr_syn_somatic_r))
SocBc=mean(as.numeric(masterdfSom_bc$cbcl_scr_syn_social_r))
SocC=mean(as.numeric(masterdfSoc_c$cbcl_scr_syn_social_r))
AttBc=mean(as.numeric(masterdfAtt_bc$cbcl_scr_syn_attention_r))
AttC=mean(as.numeric(masterdfAtt_c$cbcl_scr_syn_attention_r))
RulBc=mean(as.numeric(masterdfRul_bc$cbcl_scr_syn_rulebreak_r))
RulC=mean(as.numeric(masterdfRul_c$cbcl_scr_syn_rulebreak_r))
AggBc=mean(as.numeric(masterdfAgg_bc$cbcl_scr_syn_aggressive_r))
AggC=mean(as.numeric(masterdfAgg_c$cbcl_scr_syn_aggressive_r))
```

``` r
# start with internalizing as 1 (corresponds to F_intFit,M_intFit,P_intFit,R_intFit) 
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=1
end=intMaxVal+1
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(intMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Internalizing', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = Ic, linetype = "dashed")+
  geom_vline(xintercept = Ibc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](Fig3_supp_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Int.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 17 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Int.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 17 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
plot_bootstraps(F_PFits[,1:MaxP],MaxP,"Internalizing",MaxP,Pbc,Pc)
```

    ## Warning in melt(t(data)): The melt generic in data.table has been passed a
    ## matrix and will attempt to redirect to the relevant reshape2 method; please
    ## note that reshape2 is deprecated, and this redirection is now deprecated as
    ## well. To continue using melt methods from reshape2 while both libraries are
    ## attached, e.g. melt.list, you can prepend the namespace like
    ## reshape2::melt(t(data)). In the next version, this warning will become an
    ## error.

    ## Warning: There was 1 warning in `summarise()`.
    ## ℹ In argument: `across(everything(), quantile, probs = c(0.01, 0.99), na.rm =
    ##   TRUE)`.
    ## Caused by warning:
    ## ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
    ## Supply arguments directly to `.fns` through an anonymous function instead.
    ## 
    ##   # Previously
    ##   across(a:b, mean, na.rm = TRUE)
    ## 
    ##   # Now
    ##   across(a:b, \(x) mean(x, na.rm = TRUE))

    ## Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in
    ## dplyr 1.1.0.
    ## ℹ Please use `reframe()` instead.
    ## ℹ When switching from `summarise()` to `reframe()`, remember that `reframe()`
    ##   always returns an ungrouped data frame and adjust accordingly.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 1 rows containing missing values (`geom_vline()`).
    ## Removed 1 rows containing missing values (`geom_vline()`).

    ## Warning: The `guide` argument in `scale_*()` cannot be `FALSE`. This was deprecated in
    ## ggplot2 3.3.4.
    ## ℹ Please use "none" instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](Fig3_supp_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

``` r
plot_bootstraps(M_PFits[,1:MaxP],MaxP,"Internalizing",MaxP,Pbc,Pc)
```

    ## Warning in melt(t(data)): The melt generic in data.table has been passed a
    ## matrix and will attempt to redirect to the relevant reshape2 method; please
    ## note that reshape2 is deprecated, and this redirection is now deprecated as
    ## well. To continue using melt methods from reshape2 while both libraries are
    ## attached, e.g. melt.list, you can prepend the namespace like
    ## reshape2::melt(t(data)). In the next version, this warning will become an
    ## error.

    ## Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in
    ## dplyr 1.1.0.
    ## ℹ Please use `reframe()` instead.
    ## ℹ When switching from `summarise()` to `reframe()`, remember that `reframe()`
    ##   always returns an ungrouped data frame and adjust accordingly.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 1 rows containing missing values (`geom_vline()`).
    ## Removed 1 rows containing missing values (`geom_vline()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(intMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(intMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Internalizing', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = Ec, linetype = "dashed")+
  geom_vline(xintercept = Ibc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Int.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 17 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Int.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 17 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
# end subscale 1
```

``` r
# now do externalizing as 2nd subscale: corresponds to F_extFit,M_extFit,P_extFit,R_extFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(extMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(extMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Externalizing', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = Ec, linetype = "dashed")+
  geom_vline(xintercept = Ebc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Ext.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 8 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Ext.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 8 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(extMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(extMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Externalizing', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = Ec, linetype = "dashed")+
  geom_vline(xintercept = Ebc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Ext.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 8 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Ext.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 8 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->

``` r
# end subscale 2
```

``` r
# now do somatic as 3rd subscale: corresponds to F_somFit,M_somFit,P_somFit,R_somFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(somMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(somMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Somatic', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = SomC, linetype = "dashed")+
  geom_vline(xintercept = SomBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Som.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 2 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Som.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 2 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(somMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(somMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Somatic', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = SomC, linetype = "dashed")+
  geom_vline(xintercept = SomBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Som.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 2 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Som.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 2 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->

``` r
# end subscale 3
```

``` r
# now do anxious depression as 4th subscale: corresponds to F_anxFit,M_anxFit,P_anxFit,R_anxFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(anxMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(anxMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Anxious Depression', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = AnxC, linetype = "dashed")+
  geom_vline(xintercept = AnxBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls AnxDep.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 6 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys AnxDep.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 6 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(anxMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(anxMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Anxious Depression', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = AnxC, linetype = "dashed")+
  geom_vline(xintercept = AnxBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. AnxDep.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 6 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. AnxDep.,')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 6 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->

``` r
# end subscale 4
```

``` r
# now do thought as 5th subscale: corresponds to F_thoFit,M_thoFit,P_thoFit,R_thoFit,
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(thoMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(thoMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Thought', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = ThoC, linetype = "dashed")+
  geom_vline(xintercept = ThoBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
# lil printout to see what max deriv is
paste(max(abs(dervPlotDf$sig_deriv)))
```

    ## [1] "0.111832973336162"

``` r
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    # Girls THOUGHT REQUIRES AN EXPANDED COLOR SCALE. NEGATIVE SLOPE IS BEYOND COLOR LIMITS OF ALL OTHER FIGURES BY .02
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.15),max(0.15)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Thought*')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Thought')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(thoMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(thoMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Thought', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = ThoC, linetype = "dashed")+
  geom_vline(xintercept = ThoBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Thought')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Thought')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-10-6.png)<!-- -->

``` r
# end subscale 5
```

``` r
# now do withdrawn depression as 6th subscale: corresponds to F_witdepFit,M_witdepFit,P_witdepFit,R_witdepFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(depMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(depMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Withdrawn Depression', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = WitC, linetype = "dashed")+
  geom_vline(xintercept = WitBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls WithDep.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys WithDep.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(depMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(depMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Withdrawn Depression', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = WitC, linetype = "dashed")+
  geom_vline(xintercept = WitBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. WithDep.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. WithDep.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->

``` r
# end subscale 6
```

``` r
# now do social as 7th subscale: corresponds to F_socFit,M_socFit,P_socFit,R_socFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(socMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(socMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Social', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = SocC, linetype = "dashed")+
  geom_vline(xintercept = SocBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Social')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Social')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(socMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(socMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Social', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = SocC, linetype = "dashed")+
  geom_vline(xintercept = SocBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Social')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Social')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->

``` r
# end subscale 7
```

``` r
# now do rule breaking as 8th subscale: corresponds to F_rulFit,M_rulFit,P_rulFit,R_rulFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(rulMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(rulMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Rule Breaking', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = RulC, linetype = "dashed")+
  geom_vline(xintercept = RulBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Rules')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 5 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Rules')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 5 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(rulMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(rulMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Rules', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = RulC, linetype = "dashed")+
  geom_vline(xintercept = RulBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Rules')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 5 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-13-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Rules')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 5 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-13-6.png)<!-- -->

``` r
# end subscale 8
```

``` r
# now do attention as 9th subscale: corresponds to F_attFit,M_attFit,P_attFit,R_attFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(attMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(attMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Attention', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = AttC, linetype = "dashed")+
  geom_vline(xintercept = AttBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Attention')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 1 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Attention')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 1 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(attMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(attMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Attention', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = AttC, linetype = "dashed")+
  geom_vline(xintercept = AttBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Attention')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 1 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-14-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Attention')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 1 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-14-6.png)<!-- -->

``` r
# end subscale 9
```

``` r
# now do aggression as 10th subscale: corresponds to F_aggFit,M_aggFit,P_aggFit,R_aggFit
# isolate Female, Male, Poor, Rich
# spans 0:maximum value of symptoms
start=end+1
end=end+(aggMaxVal+1)
F_PFits=Fits[,start:end]

# calculate some background info for plots
MaxP_F=find_furthest_nonzero(F_PFits)

# isolate male fits
start=end+1
end=end+(aggMaxVal+1)
M_PFits=Fits[start:end]

# calculate some background info for plots
MaxP_M=find_furthest_nonzero(M_PFits)
# use lowest common value for plots
MaxP <- min(MaxP_M, MaxP_F)

# get median value Girls
F_PFits_Coverage=F_PFits[,seq(1:MaxP)]
col_means=colMeans(F_PFits_Coverage)
FP_medians <- apply(F_PFits_Coverage, 2, median)

# get median value boy
M_PFits_Coverage=M_PFits[,seq(1:MaxP)]
col_means=colMeans(M_PFits_Coverage)
MP_medians <- apply(M_PFits_Coverage, 2, median)
MP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_girls = FP_medians,
  y_boys = MP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_boys)) +
  geom_line(aes(), color = "#fbad24", size = 3) +
  geom_line(aes(y=y_girls),color = "#923eb5", size = 3) +
  labs(x = 'Aggression', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = AggC, linetype = "dashed")+
  geom_vline(xintercept = AggBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
F_P_derivative_matrix <- matrix(0, nrow = nrow(F_PFits), ncol = ncol(F_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- F_PFits[, i + 1] - F_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  F_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
M_P_derivative_matrix <- matrix(0, nrow = nrow(M_PFits), ncol = ncol(M_PFits) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(F_PFits) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- M_PFits[, i + 1] - M_PFits[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  M_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(F_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(F_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(F_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Girls Aggr.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(M_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(M_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(M_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Boys Aggr.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
# poverty analysis: start with below poverty line group
start=end+1
end=end+(aggMaxVal+1)
povBoots=Fits[,start:end]
# and nonpov
start=end+1
end=end+(aggMaxVal+1)
nonpovBoots=Fits[,start:end]

# calculate some background info for plots
MaxP_P=find_furthest_nonzero(povBoots)
MaxP_R=find_furthest_nonzero(nonpovBoots)

# use lowest common value for plots
MaxP <- min(MaxP_P, MaxP_R)

# get median value pov
P_PFits_Coverage=povBoots[,seq(1:MaxP)]
col_means=colMeans(P_PFits_Coverage)
PP_medians <- apply(P_PFits_Coverage, 2, median)

# get median value nonpov
R_PFits_Coverage=nonpovBoots[,seq(1:MaxP)]
col_means=colMeans(R_PFits_Coverage)
RP_medians <- apply(M_PFits_Coverage, 2, median)
RP_medians=MP_medians[1:MaxP]

# merge data
data <- data.frame(
  x = 0:(MaxP-1),
  y_pov = PP_medians,
  y_nonpov = RP_medians
)

# Create the line plot for p
ggplot(data, aes(x = x, y = y_pov)) +
  geom_line(aes(), color = "#FF5003", size = 3) +
  geom_line(aes(y=y_nonpov),color = "#003F7D", size = 3) +
  labs(x = 'Aggression', y = expression(italic(g))) +
  theme_minimal(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"))+ylim(-1.5,1.5)+
  geom_vline(xintercept = AggC, linetype = "dashed")+
  geom_vline(xintercept = AggBc, linetype = "dashed")+
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
        scale_x_continuous(limits = c(0,MaxP-1),expand = expansion(mult = c(0, 0)))
```

![](Fig3_supp_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

``` r
# and derivatives
# Create an empty matrix to store the derivatives
P_P_derivative_matrix <- matrix(0, nrow = nrow(povBoots), ncol = ncol(povBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(povBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- povBoots[, i + 1] - povBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  P_P_derivative_matrix[, i] <- derivatives
}

# Create an empty matrix to store the derivatives
R_P_derivative_matrix <- matrix(0, nrow = nrow(nonpovBoots), ncol = ncol(nonpovBoots) - 1)

# Calculate the derivative for each column
for (i in 1:(ncol(nonpovBoots) - 1)) {
  # Calculate the differences in x (assuming a constant difference)
  dx <- 1
  # Calculate the differences in y (predicted values)
  dy <- nonpovBoots[, i + 1] - nonpovBoots[, i]
  # Calculate the derivatives (slopes)
  derivatives <- dy / dx
  # Store the derivatives in the derivative matrix
  R_P_derivative_matrix[, i] <- derivatives
}

# calc sig dervs
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(P_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(P_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(P_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5)) ###  ???
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('Pov. Aggr.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->

``` r
# get straightfoward of segment where 99% is over 0 or under
positive_counts <- colSums(R_P_derivative_matrix > 0, na.rm = TRUE)
negative_counts <- colSums(R_P_derivative_matrix < 0, na.rm = TRUE)
# find where each is 99% or greater
positive_countsSig=positive_counts>9900
negative_countsSig=negative_counts>9900
# make dataframe: 50th percentile of derivatives accompanied by posSig and NegSig vector
data <- apply(R_P_derivative_matrix, 2, function(x) quantile(x, probs = 0.5))
dervPlotDf<-data.frame(data,positive_countsSig,negative_countsSig)
# if either is sig at 99% plot
dervPlotDf$sig_derivMask=dervPlotDf[,2]+dervPlotDf[,3]>0
# use it to mask calculated derivs
dervPlotDf$sig_deriv=0
dervPlotDf$sig_deriv[dervPlotDf$sig_derivMask]=dervPlotDf$data[dervPlotDf$sig_derivMask]
dervPlotDf$seq=1:(dim(dervPlotDf)[1])
ggplot(data=dervPlotDf) + geom_raster(aes(x = seq, y = .5, fill = sig_deriv))+
    theme(panel.spacing = unit(-.01,"cm")) +
    scale_fill_gradientn(colors = my_palette(100),limits = c(min(-.1),max(0.1)))+theme_minimal(base_size = 35)+
    xlim(c(0,MaxP))+xlab('NonPov. Aggr.')+
    guides(fill=FALSE)+
    theme(axis.title.y = element_blank(),axis.text.y=element_blank())+theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
    scale_x_continuous(limits = c(0,MaxP),expand = expansion(mult = c(0, 0)))
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: Removed 4 rows containing missing values (`geom_raster()`).

![](Fig3_supp_files/figure-gfm/unnamed-chunk-15-6.png)<!-- -->

``` r
# end subscale 10
```
