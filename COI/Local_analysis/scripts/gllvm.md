GLLVM with richness data, COI
================
Emil Ellegaard Thomassen
\`r Sys.Date()

The following script performs Generalized linear latent variable models
of species richness for each order (gllvm).

The package gllvm is used.
(<https://jenniniku.github.io/gllvm/index.html>) (Niku *et al.*, 2019)

# Data:

First, all dung-associated species are selected, as other species would
represent random visits in dung pads, or contamination through airborne
DNA or carry-over by other organisms.

Abund: Abundance matrix. Matrix with samples as rows and orders as
columns. “Abundance values” are species counts (including ASVs
identified at higher taxonomic levels, e.g., “Aphodius sp.12
(sphacelatus/prodromus))

Envs: Environmental variables. Data frame with samples as rows, and the
following environmental variables as columns: - Source (values:
Exmoor/Galloway) - Season (values: Early/Mid/Late) - Habitat (values:
Forest/Open)

# Setup

First, R environment is set up, and packages are loaded:

# Load data

Data is loaded and re-organized to match input to the gllvm() function:

``` r
#Load dataset for extracting environmental variables
#dataset<-readRDS("COI/Local_analysis/output/RDS_files/dataset_final_all_dung")

dataset<-readRDS("../output/RDS_files/dataset_final_all_dung") #only when knitting


#Load richness data
#abund<-readRDS("COI/Local_analysis/output/RDS_files/richness_each_order")

abund<-readRDS("../output/RDS_files/richness_each_order") #only when knitting


#Extract environmental variables
envs<-data.frame("source"=dataset@samples$source, "season"=dataset@samples$season, "habitat"=dataset@samples$habitat, "month"=dataset@samples$month)

#Change season-column to factor with levels "Early", "Mid" & "Late" instead of spring, summer, autumn and winter

envs<-envs %>%
  dplyr::mutate(season=factor(season, levels=c("Early","Mid","Late")))

envs$season[which(envs$month=="January")]<-"Early"
envs$season[which(envs$month=="February")]<-"Early"
envs$season[which(envs$month=="March")]<-"Early"
envs$season[which(envs$month=="April")]<-"Early"
envs$season[which(envs$month=="May")]<-"Mid"
envs$season[which(envs$month=="June")]<-"Mid"
envs$season[which(envs$month=="July")]<-"Mid"
envs$season[which(envs$month=="August")]<-"Mid"
envs$season[which(envs$month=="September")]<-"Late"
envs$season[which(envs$month=="October")]<-"Late"
envs$season[which(envs$month=="November")]<-"Late"
envs$season[which(envs$month=="December")]<-"Late"

#Change class of abundance data to integers
abund<-as.data.frame(abund)
abund<-lapply(abund, as.integer)
abund<-as.data.frame(abund)

#Set rownames of both data frames to sample names
rownames(envs)<-dataset@samples$name
rownames(abund)<-dataset@samples$name

#Create list of abundance data and environmental variables
df_list<-list("abund"=abund,"env"=envs)

#Drop the month column in environmental variables, as it was only used to change season values
df_list$env<-df_list$env[,c(1:3)]

#Drop "control" as factor level for source column (remnant from the initial filtering steps)
df_list$env$source<-droplevels(df_list$env$source)


#Define inputs (y=abundance matrix, X=data frame of environmental variables)
y <- as.matrix(df_list$abund)
X <- df_list$env
```

# Modeling - GLLVM

First, modeling is performed with different distributions

## Gaussian distribution

``` r
# Gaussian
fitG <- gllvm(y = y, X = X, formula = y ~ source + season + habitat, family = "gaussian")

par(mfrow = c(2,2))
plot(fitG)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20gaussian-1.png)<!-- -->

``` r
AIC.fitG<-AIC(fitG)
AIC.fitG
```

    ## [1] 5817

![](gllvm_files/figure-gfm/GLLVM%20-%20gaussian-2.png)<!-- -->

AIC value is 5816.958 - diagnostics plots does not look good!

## Poisson distribution

``` r
# Poisson distribution
fitP <- gllvm(y = y, X = X, formula = y ~ source + season + habitat, family = poisson())

par(mfrow = c(2,2))
plot(fitP)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20poisson-1.png)<!-- -->

``` r
AIC.fitP<-AIC(fitP)
AIC.fitP
```

    ## [1] 6937

![](gllvm_files/figure-gfm/GLLVM%20-%20poisson-2.png)<!-- -->

AIC value is 6936.7299 - diagnostics plots looks better, but seems to be
an overdispersion problem!

## Negative binomial distibution

``` r
# Negative binomial
fitNB <- gllvm(y = y, X = X, formula = y ~ source + season + habitat, family = "negative.binomial")

par(mfrow = c(2,2))
plot(fitNB)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20negative%20binomial-1.png)<!-- -->

``` r
AIC.fitNB<-AIC(fitNB)
AIC.fitNB
```

    ## [1] 7171

![](gllvm_files/figure-gfm/GLLVM%20-%20negative%20binomial-2.png)<!-- -->

AIC value is 7170.7291 - diagnostics plots shows that NB model seems to
have taken care of the overdispersion problem!

# Test of number of latent variables

Now several models with increasing number of latent variables are tested
to find the optimum

``` r
AICs.lv.test<-NULL
LVs<-NULL
for (i in c(1:5)) {
  mod<-gllvm(y = y, X = X, formula = y ~ source + season + habitat, family = "negative.binomial", num.lv = i)
  AICs.lv.test<-append(AICs.lv.test, AIC(mod))
  LVs<-append(LVs, i)
}

par(mfrow=c(1,1))
plot(LVs,AICs.lv.test, main="AIC with increasing Latent variables", xlab="Number of LVs", ylab="AIC", col="red", pch=19)
```

![](gllvm_files/figure-gfm/Test%20of%20number%20of%20LVs-1.png)<!-- -->

It seems like AIC increases almost linearly with increasing latent
variables. I am not sure how to make sense of this, as i would have
expected AIC to decrease with increasing number of latent variables? But
seems like 1 would be the best model from AIC.

# Including interactions between environmental variables

Try to include interactions between predictors in the model

``` r
# Including interaction term
fitNB.int<-gllvm(y = y, X = X, formula = y ~ source * season * habitat, family = "negative.binomial", num.lv = 1)

par(mfrow = c(2,2))
plot(fitNB.int)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20NB%20and%20interactions-1.png)<!-- -->

``` r
AIC.fitNB.int<-AIC(fitNB.int)

#Compare to model without interactions
fitNB.noint<-gllvm(y = y, X = X, formula = y ~ source + season + habitat, family = "negative.binomial", num.lv = 1)

AIC.fitNB.noint<-AIC(fitNB.noint)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20NB%20and%20interactions-2.png)<!-- -->

AIC without interaction is: 7134.7299 AIC with interaction is: 7220.1617

Higher AIC indicate that including interactions does not improve the
model.

# Test if some orders migth be outliers

If some orders are only present in very few samples, they might
disproportionately drive the analysis. Try to investigate this, and
remove outlier orders.

``` r
#Calculate number of samples with presence of each order
nsamples<-NULL
for (i in 1:ncol(y)) {nsamples<-append(nsamples,length(y[which(y[,i]>0),i]))}

out<-cbind(colnames(y),nsamples,colSums(y))
out
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
nsamples
</th>
<th style="text-align:left;">
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Coleoptera
</td>
<td style="text-align:left;">
Coleoptera
</td>
<td style="text-align:left;">
184
</td>
<td style="text-align:left;">
618
</td>
</tr>
<tr>
<td style="text-align:left;">
Strongylida
</td>
<td style="text-align:left;">
Strongylida
</td>
<td style="text-align:left;">
228
</td>
<td style="text-align:left;">
2519
</td>
</tr>
<tr>
<td style="text-align:left;">
Chromadorea.unknown.order
</td>
<td style="text-align:left;">
Chromadorea.unknown.order
</td>
<td style="text-align:left;">
97
</td>
<td style="text-align:left;">
225
</td>
</tr>
<tr>
<td style="text-align:left;">
Diptera
</td>
<td style="text-align:left;">
Diptera
</td>
<td style="text-align:left;">
194
</td>
<td style="text-align:left;">
731
</td>
</tr>
<tr>
<td style="text-align:left;">
Rhabditida
</td>
<td style="text-align:left;">
Rhabditida
</td>
<td style="text-align:left;">
191
</td>
<td style="text-align:left;">
774
</td>
</tr>
<tr>
<td style="text-align:left;">
Plectida
</td>
<td style="text-align:left;">
Plectida
</td>
<td style="text-align:left;">
19
</td>
<td style="text-align:left;">
27
</td>
</tr>
<tr>
<td style="text-align:left;">
Entomobryomorpha
</td>
<td style="text-align:left;">
Entomobryomorpha
</td>
<td style="text-align:left;">
90
</td>
<td style="text-align:left;">
135
</td>
</tr>
<tr>
<td style="text-align:left;">
Poduromorpha
</td>
<td style="text-align:left;">
Poduromorpha
</td>
<td style="text-align:left;">
108
</td>
<td style="text-align:left;">
167
</td>
</tr>
<tr>
<td style="text-align:left;">
Symphypleona
</td>
<td style="text-align:left;">
Symphypleona
</td>
<td style="text-align:left;">
32
</td>
<td style="text-align:left;">
34
</td>
</tr>
<tr>
<td style="text-align:left;">
Phthiraptera
</td>
<td style="text-align:left;">
Phthiraptera
</td>
<td style="text-align:left;">
17
</td>
<td style="text-align:left;">
17
</td>
</tr>
<tr>
<td style="text-align:left;">
Sarcoptiformes
</td>
<td style="text-align:left;">
Sarcoptiformes
</td>
<td style="text-align:left;">
33
</td>
<td style="text-align:left;">
37
</td>
</tr>
<tr>
<td style="text-align:left;">
Nematoda.unknown.order
</td>
<td style="text-align:left;">
Nematoda.unknown.order
</td>
<td style="text-align:left;">
35
</td>
<td style="text-align:left;">
42
</td>
</tr>
<tr>
<td style="text-align:left;">
Mesostigmata
</td>
<td style="text-align:left;">
Mesostigmata
</td>
<td style="text-align:left;">
30
</td>
<td style="text-align:left;">
41
</td>
</tr>
<tr>
<td style="text-align:left;">
Polydesmida
</td>
<td style="text-align:left;">
Polydesmida
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Trombidiformes
</td>
<td style="text-align:left;">
Trombidiformes
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
Lithobiomorpha
</td>
<td style="text-align:left;">
Lithobiomorpha
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Neelipleona
</td>
<td style="text-align:left;">
Neelipleona
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Julida
</td>
<td style="text-align:left;">
Julida
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Insecta.unknown.order
</td>
<td style="text-align:left;">
Insecta.unknown.order
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
</tr>
</tbody>
</table>

There seems to be several orders that are only present in very few
samples. Orders present in fewer than 20 samples are removed

``` r
y.subs<-y[,which(nsamples>20)]
```

Now, modeling is performed with the dataset without the outlier orders

# Modeling - without outliers

## Poisson

``` r
# Poisson distribution
sfitP <- gllvm(y = y.subs, X = X, formula = y ~ source + season + habitat, family = poisson())

par(mfrow = c(2,2))
plot(sfitP)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20without%20outliers%20-%20poisson-1.png)<!-- -->

``` r
AIC.sfitP<-AIC(sfitP)
AIC.sfitP
```

    ## [1] 6314

![](gllvm_files/figure-gfm/GLLVM%20-%20without%20outliers%20-%20poisson-2.png)<!-- -->

Now a poisson distribution looks better than before - seems that the
zero-inflation might be caused by some of the underrepresented orders,
which are now removed as outliers.

AIC is 6314.4118

Still maybe a little overdispersion, so negative binomial is tried.

## Negative binomial

``` r
# Negative binomial
sfitNB <- gllvm(y = y.subs, X = X, formula = y ~ source + season + habitat, family = "negative.binomial")

par(mfrow = c(2,2))
plot(sfitNB)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20without%20outliers%20-%20negative%20binomial-1.png)<!-- -->

``` r
AIC.sfitNB<-AIC(sfitNB)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20without%20outliers%20-%20negative%20binomial-2.png)<!-- -->

Looks generally good from diagnostics plots. But AIC (6600.9281) is
higher than the poisson model. Thus the poisson model is chosen.

# Test number of latent variables

Models with increasing number of LVs are performed, and AICs are
plotted.

``` r
AICs.lv.test<-NULL
LVs<-NULL
for (i in c(1:5)) {
  mod<-gllvm(y = y.subs, X = X, formula = y ~ source + season + habitat, family = poisson(), num.lv = i)
  AICs.lv.test<-append(AICs.lv.test, AIC(mod))
  LVs<-append(LVs, i)
}

par(mfrow=c(1,1))
plot(LVs,AICs.lv.test, main="AIC with increasing Latent variables (without outliers)", xlab="Number of LVs", ylab="AIC", col="red", pch=19)
```

![](gllvm_files/figure-gfm/Test%20number%20of%20lvs%20-%20without%20outliers-1.png)<!-- -->

Including 3 latent variables seems to be the best model.

# Including interactions

Test if including interactions between predictors improve model
performance

``` r
# Including interaction term
sfitP.int<-gllvm(y = y.subs, X = X, formula = y ~ source * season * habitat, family = poisson(), num.lv = 3)

par(mfrow = c(2,2))
plot(sfitP.int)
```

![](gllvm_files/figure-gfm/Test%20-%20including%20interactions%20-%20without%20outliers-1.png)<!-- -->

``` r
AIC.sfitP.int<-AIC(sfitP.int)

sfitP.noint<-gllvm(y = y.subs, X = X, formula = y ~ source + season + habitat, family = poisson(), num.lv = 3)
AIC.sfitP.noint<-AIC(sfitP.noint)
```

![](gllvm_files/figure-gfm/Test%20-%20including%20interactions%20-%20without%20outliers-2.png)<!-- -->

AIC without interaction: 6315.2341 AIC with interaction: 6337.3332

Inclusion of interactions does not seem to improve model performance.

# Inclusion of random row effect

``` r
#Random row effect included
sfitP.lv3.RX<-gllvm(y = y.subs, X = X, formula = y ~ source + season + habitat, family = poisson(), num.lv = 3, row.eff = "random")

par(mfrow = c(2,2))
plot(sfitP.lv3.RX)                    
```

![](gllvm_files/figure-gfm/GLLVM%20-%20random%20row%20effect-1.png)<!-- -->

``` r
AIC.sfitP.lv3.RX<-AIC(sfitP.lv3.RX)
```

![](gllvm_files/figure-gfm/GLLVM%20-%20random%20row%20effect-2.png)<!-- -->

Inclusion of random row effect does not increase model performance:

AIC with row-effects: 6334.3531 AIC without: 6315.2341

Maybe we should include anyway, when our samples are not truly
independent?

So far, model with poisson distribution, without outlier orders,
including 2 latent variables and no random effects or interactions are
chosen as the best model

# Best model

``` r
#Best model

fitBEST<-gllvm(y = y.subs, X = X, formula = y ~ source + season + habitat, family = poisson(), num.lv = 3)

sumBEST<-summary(fitBEST)
sumBEST
```

    ## 
    ## Call:
    ## gllvm(y = y.subs, X = X, formula = y ~ source + season + habitat, 
    ##     num.lv = 3, family = poisson())
    ## 
    ## Family:  poisson 
    ## 
    ## AIC:  6454 AICc:  6460 BIC:  6749 LL:  -3142 df:  85 
    ## 
    ## Constrained LVs:  0 
    ## Reduced Ranks:  0 
    ## Unconstrained LVs:  3 
    ## Standard deviation of LVs:  0 0.457 0.684 
    ## 
    ## Formula:  ~ysource + season + habitat 
    ## LV formula:  ~yNULL 
    ## 
    ## Coefficients predictors:
    ##                                          Estimate Std. Error z value Pr(>|z|)    
    ## sourcegalloway:Coleoptera                  0.0729     0.0808    0.90   0.3664    
    ## sourcegalloway:Strongylida                -1.5590     0.0838  -18.61  < 2e-16 ***
    ## sourcegalloway:Chromadorea.unknown.order   1.1989     0.1981    6.05  1.4e-09 ***
    ## sourcegalloway:Diptera                     0.9058     0.0987    9.18  < 2e-16 ***
    ## sourcegalloway:Rhabditida                  0.6624     0.1055    6.28  3.4e-10 ***
    ## sourcegalloway:Entomobryomorpha           -0.2352     0.2017   -1.17   0.2435    
    ## sourcegalloway:Poduromorpha               -0.0992     0.1751   -0.57   0.5711    
    ## sourcegalloway:Symphypleona                0.8266     0.3674    2.25   0.0244 *  
    ## sourcegalloway:Sarcoptiformes              0.5961     0.3607    1.65   0.0984 .  
    ## sourcegalloway:Nematoda.unknown.order     -0.5210     0.3744   -1.39   0.1641    
    ## sourcegalloway:Mesostigmata                1.0960     0.3915    2.80   0.0051 ** 
    ## seasonMid:Coleoptera                       1.0851     0.1063   10.21  < 2e-16 ***
    ## seasonMid:Strongylida                      0.5830     0.0987    5.91  3.5e-09 ***
    ## seasonMid:Chromadorea.unknown.order        2.8518     0.3914    7.29  3.2e-13 ***
    ## seasonMid:Diptera                          1.7467     0.1424   12.27  < 2e-16 ***
    ## seasonMid:Rhabditida                       1.4775     0.1537    9.61  < 2e-16 ***
    ## seasonMid:Entomobryomorpha                 0.1975     0.2575    0.77   0.4430    
    ## seasonMid:Poduromorpha                     0.4076     0.2341    1.74   0.0816 .  
    ## seasonMid:Symphypleona                    -1.3781     0.5596   -2.46   0.0138 *  
    ## seasonMid:Sarcoptiformes                   0.7569     0.4570    1.66   0.0977 .  
    ## seasonMid:Nematoda.unknown.order           2.9412     1.0460    2.81   0.0049 ** 
    ## seasonMid:Mesostigmata                    14.3573     0.2540   56.53  < 2e-16 ***
    ## seasonLate:Coleoptera                      0.1649     0.1244    1.33   0.1850    
    ## seasonLate:Strongylida                     0.8386     0.0984    8.53  < 2e-16 ***
    ## seasonLate:Chromadorea.unknown.order       1.9952     0.4027    4.95  7.3e-07 ***
    ## seasonLate:Diptera                         0.9648     0.1508    6.40  1.6e-10 ***
    ## seasonLate:Rhabditida                      1.4400     0.1538    9.36  < 2e-16 ***
    ## seasonLate:Entomobryomorpha                0.4538     0.2479    1.83   0.0672 .  
    ## seasonLate:Poduromorpha                    0.7198     0.2253    3.19   0.0014 ** 
    ## seasonLate:Symphypleona                   -0.1842     0.3666   -0.50   0.6154    
    ## seasonLate:Sarcoptiformes                  0.4350     0.4836    0.90   0.3683    
    ## seasonLate:Nematoda.unknown.order          2.9262     1.0479    2.79   0.0052 ** 
    ## seasonLate:Mesostigmata                   12.2078     0.3651   33.44  < 2e-16 ***
    ## habitatopen:Coleoptera                     0.1734     0.0808    2.15   0.0318 *  
    ## habitatopen:Strongylida                    0.0759     0.0772    0.98   0.3252    
    ## habitatopen:Chromadorea.unknown.order      0.0382     0.1852    0.21   0.8367    
    ## habitatopen:Diptera                       -0.0832     0.0951   -0.87   0.3817    
    ## habitatopen:Rhabditida                     0.1592     0.1042    1.53   0.1265    
    ## habitatopen:Entomobryomorpha              -0.9615     0.2147   -4.48  7.5e-06 ***
    ## habitatopen:Poduromorpha                   0.7745     0.1841    4.21  2.6e-05 ***
    ## habitatopen:Symphypleona                   0.2464     0.3462    0.71   0.4765    
    ## habitatopen:Sarcoptiformes                -0.0550     0.3543   -0.16   0.8766    
    ## habitatopen:Nematoda.unknown.order        -0.5033     0.3676   -1.37   0.1709    
    ## habitatopen:Mesostigmata                   0.4805     0.3705    1.30   0.1946    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Create coefficient plots

First input is created from the summary of the model, and by extracting
confidence intervals

``` r
#Gather all parameters and confidence intervals in df

conf<-as.data.frame(confint(fitBEST))



df.params<-as.data.frame(cbind(fitBEST$params$beta0,fitBEST$params$Xcoef,
                               conf$`2.5 %`[which(grepl("galloway",rownames(conf)))],
                               conf$`97.5 %`[which(grepl("galloway",rownames(conf)))],
                               conf$`2.5 %`[which(grepl("Mid",rownames(conf)))],
                               conf$`97.5 %`[which(grepl("Mid",rownames(conf)))],
                               conf$`2.5 %`[which(grepl("Late",rownames(conf)))],
                               conf$`97.5 %`[which(grepl("Late",rownames(conf)))],
                               conf$`2.5 %`[which(grepl("open",rownames(conf)))],
                               conf$`97.5 %`[which(grepl("open",rownames(conf)))]
))
                               
                               
colnames(df.params)<-c("intercept", "coef.source.GAL","coef.season.MID","coef.season.LATE", "coef.habitat.OPEN", "conf.low.source.GAL","conf.high.source.GAL","conf.low.season.MID","conf.high.season.MID","conf.low.season.LATE","conf.high.season.LATE", "conf.low.habitat.OPEN", "conf.high.habitat.OPEN")

#Check if zero is included in 95% confint
signif.GAL<-NULL
for(i in 1:nrow(df.params)) {if (between(0, df.params$conf.low.source.GAL[i], df.params$conf.high.source.GAL[i])) {signif.GAL<-append(signif.GAL,"ns.")}
                             else {if (df.params$conf.high.source.GAL[i]<0) {signif.GAL<-append(signif.GAL,"higher richness horse")}
                               else {signif.GAL<-append(signif.GAL,"higher richness cattle")}
                               }
}

signif.MID<-NULL
for(i in 1:nrow(df.params)) {if (between(0, df.params$conf.low.season.MID[i], df.params$conf.high.season.MID[i])) {signif.MID<-append(signif.MID,"ns.")}
                             else {if (df.params$conf.high.season.MID[i]<0) {signif.MID<-append(signif.MID,"higher richness early")}
                               else {signif.MID<-append(signif.MID,"higher richness mid")}
                               }
}

signif.LATE<-NULL
for(i in 1:nrow(df.params)) {if (between(0, df.params$conf.low.season.LATE[i], df.params$conf.high.season.LATE[i])) {signif.LATE<-append(signif.LATE,"ns.")}
                             else {if (df.params$conf.high.season.LATE[i]<0) {signif.LATE<-append(signif.LATE,"higher richness early")}
                               else {signif.LATE<-append(signif.LATE,"higher richness late")}
                               }
}

signif.OPEN<-NULL
for(i in 1:nrow(df.params)) {if (between(0, df.params$conf.low.habitat.OPEN[i], df.params$conf.high.habitat.OPEN[i])) {signif.OPEN<-append(signif.OPEN,"ns.")}
                             else {if (df.params$conf.high.habitat.OPEN[i]<0) {signif.OPEN<-append(signif.OPEN,"higher richness forest")}
                               else {signif.OPEN<-append(signif.OPEN,"higher richness open")}
                               }
}


df.params<-cbind(df.params, "signif.GAL"=as.factor(signif.GAL),"signif.season.MID"=as.factor(signif.MID),"signif.season.LATE"=as.factor(signif.LATE),"signif.OPEN"=as.factor(signif.OPEN))
```

Then coefficient plots are created with ggplot2

``` r
#Plotting coefficient plot of source

df.source<-df.params[,which(grepl("GAL", colnames(df.params)))]

level_order<-rownames(df.source[order(df.source$coef.source.GAL, decreasing = FALSE),])

colnames(df.source)<-c("Coefficient","Conf.low","Conf.high","Significance")


ggcoefplot.source<-ggplot(df.source, aes(x=Coefficient, y=factor(rownames(df.source),level=level_order),col=Significance)) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  coord_cartesian(xlim=c(-2,2)) +
  geom_point() +
  scale_color_manual(values=c("darkseagreen4","dodgerblue3","grey")) +
  geom_segment(aes(x=Conf.low, xend=Conf.high, yend=factor(rownames(df.source)),colour=Significance)) +
  xlab("Coefficients") +
  ylab("Orders") +
  theme_bw()


ggcoefplot.source<-ggcoefplot.source +  theme(legend.position="top")
  
ggcoefplot.source
```

![](gllvm_files/figure-gfm/Plot%20-%20source-1.png)<!-- -->

``` r
#Plotting coefficient plot of source

df.season<-df.params[,which(grepl("season", colnames(df.params)))]

level_order1<-rownames(df.season[order(df.season$coef.season.MID, decreasing = FALSE),])
level_order2<-rownames(df.season[order(df.season$coef.season.LATE, decreasing = FALSE),])


colnames(df.season)<-c("Coefficient.MID","Coefficient.LATE","Conf.low.MID","Conf.high.MID","Conf.low.LATE","Conf.high.LATE","Significance", "Significance.LATE" )


ggcoefplot.season<-ggplot(df.season, aes(x=Coefficient.MID, y=factor(rownames(df.source),level=level_order1),col=Significance)) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  coord_cartesian(xlim=c(-3,20)) +
  geom_point(position = position_nudge(x=0, y=0.1)) +
  geom_segment(aes(x=Conf.low.MID, xend=Conf.high.MID, yend=factor(rownames(df.source)),colour=Significance), position = position_nudge(x=0, y=0.1)) +
  
  scale_color_manual(values=c("#4fe7c7","#83370d","#fb7264","grey")) +
  
  geom_point(data=df.season, aes(x=Coefficient.LATE, y=factor(rownames(df.source),level=level_order),col=Significance.LATE), position = position_nudge(x=0, y=-0.1)) +
  geom_segment(data=df.season,aes(x=Conf.low.LATE, xend=Conf.high.LATE, yend=factor(rownames(df.source)),colour=Significance.LATE), position = position_nudge(x=0, y=-0.1)) +
  
  
  xlab("Coefficients") +
  ylab("Orders") +
  theme(legend.position="top") +
  theme_bw()

ggcoefplot.season<-ggcoefplot.season +  theme(legend.position="top")
  
ggcoefplot.season
```

![](gllvm_files/figure-gfm/Plot%20-%20season-1.png)<!-- -->

``` r
#Plotting coefficient plot of source

df.habitat<-df.params[,which(grepl("OPEN", colnames(df.params)))]

level_order<-rownames(df.habitat[order(df.habitat$coef.habitat.OPEN, decreasing = FALSE),])

colnames(df.habitat)<-c("Coefficient","Conf.low","Conf.high","Significance")


ggcoefplot.habitat<-ggplot(df.habitat, aes(x=Coefficient, y=factor(rownames(df.habitat),level=level_order),col=Significance)) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  coord_cartesian(xlim=c(-2,2)) +
  geom_point() +
  scale_color_manual(values=c("#AED0A5","#F9F0BC","grey")) +
  geom_segment(aes(x=Conf.low, xend=Conf.high, yend=factor(rownames(df.habitat)),colour=Significance)) +
  xlab("Coefficients") +
  ylab("Orders") +
  theme(legend.position="top") +
  theme_bw()

ggcoefplot.habitat<-ggcoefplot.habitat +  theme(legend.position="top")
  
ggcoefplot.habitat
```

![](gllvm_files/figure-gfm/Plot%20-%20habitat-1.png)<!-- -->

# Arrange and export plots

``` r
gg.tmp<-ggarrange(ggcoefplot.source,ggcoefplot.habitat, labels = c("A","B"), nrow = 1, ncol=2)

ggarrange(gg.tmp,ggcoefplot.season, labels = c("","C"), nrow = 2, ncol=1)
```

![](gllvm_files/figure-gfm/Export%20plots-1.png)<!-- -->

``` r
#ggsave("COI/Local_analysis/output/plots/Coefficient_plots_gllvm.pdf", width=14, height=8)
```
