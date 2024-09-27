MA Fitness FA24
================
Rutuja
2024-09-23

# Setup functions

### This fits a spline to y=OD versus x=time, and then finds the maximum slope

``` r
spline.slope<-function(x, y, n=101, eps=1e-5, span=0.075){
  max(nderiv(loess(log(y) ~ x, degree=1, span=span), x), na.rm=TRUE)
}
```

### Used by the function above to get a local (linear) slope around a point

``` r
nderiv <- function(fit, x, eps=1e-5){
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}
```

## Spline fitting:

First we start by converting everything to log scale. During the
exponential growth phase, log(N) is proportional to t where N is the
population and t is time. There is noise in the original data. We remove
the noise and try to obtain the underlying curve using loess. It gives a
smooth curve. Then we use that curve to find local slopes. nderiv gives
precise estimate of the local linear slope at every point.

- this fits a spline to y=OD versus x=time, and then finds the
  time-stamp for maximum slope

``` r
spline.time<-function(x, y, n=101, eps=1e-5, span=0.075){
  estimates <- loess(log(y) ~ x, degree=1, span=span)
  slopes <- nderiv(estimates, x)
  return(which.max(slopes))
}
```

**Note that the unit of the time stamp is the number of 15 minute
intervals from the beginning of measurement.**

## Figuring out the parameters for the splines

``` r
d <- read.csv("data/Ref 04 28.csv")
assay.data <- read.delim("data/Rutuja 04 06.txt")
```

``` r
well <- assay.data$D4
time <- seq(1, length(well))


ggplot() + 
  geom_point(aes(x=time, y=well))
```

    ## Warning: Removed 43 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
## Now trying to smooth 
smoothed <- predict(loess(log(well) ~ time, degree=1, span=0.05), time)
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : k-d tree limited by memory. ncmax= 200

``` r
ggplot() + 
  geom_point(aes(x=time, y=smoothed), color="red") +
  geom_point(aes(x=time, y=log(well)))
```

    ## Warning: Removed 43 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 43 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
# slopes <- nderiv(log(well), time)
# fitted.slopes <- nderiv(log(well), time)
```

## Reading and cleaning the data

Load all files and setup variables

This is just some information from the experiment. Using the info to
make sure everything is labelled right and is categorized right
downstream.

``` r
dates.1 <- seq(mdy(04062023), mdy(04102023), 1)
dates.2 <- seq(mdy(04282023), mdy(11302023), 1)
dates.2 <- c(dates.2, mdy(11022023))

ancestors <- c('H1', 'H2', 'H3', 'D1', 'D2', 'D3')

# these are all the haploid ancestors including the diploids that are found to haploids
ancestors.haploid <- c('H1', 'H2', 'H3', 'D2', 'D3', 'C1', 'C3')

# these are haploids that were intended to be haploids
anc.hap.og <- c('H1', 'H2', 'H3', 'C1', 'C2')

# these are haploids that were intended to be diploids
anc.hap.fake <- c('D2', 'D3')

# this is the only actual diploid control that was diploid
ancestors.diploid <- c('D1')

# this MA line started out as a haploid but ended up as a diploid
fake.haploids <- c(81)

# these MA lines started out as diploids but ended up as haploids
fake.diploids <- c(20, 28, 48, 52, 54, 84, 100)

# these have an aneuploidy of chromosome 3
aneuploids <- c(45, 68, 98)
```

Adding additional information for categorization. Excluding MA lines
that changed ploidy. Excluding batch 1 due to lack of diploid controls.

``` r
d$batch[d$date %in% dates.1] <- 1
d$batch[d$date %in% dates.2] <- 2

d <- d %>% filter(batch == 2)

d$label <- d$treatment
d$label[d$treatment == 'Blank'] <- '0'
d$label[d$treatment == 'H1'] <- '101'
d$label[d$treatment == 'D1'] <- '102'
d$label[d$treatment == 'H2'] <- '103'
d$label[d$treatment == 'D2'] <- '104'
d$label[d$treatment == 'H3'] <- '105'
d$label[d$treatment == 'D3'] <- '106'
d$label[d$treatment == 'C1'] <- '107'
d$label[d$treatment == 'C3'] <- '109'
d$label <- as.numeric(d$label)

d <- d %>% mutate(category = case_when(label == 0 ~ 'Blank',
                              label > 100 & treatment %in% ancestors.haploid ~ 'Ctrl.H',
                              label > 100 & treatment %in% ancestors.diploid ~ 'Ctrl.D',
                              label %% 2 == 0 ~ 'MA.D',
                              label %% 2 == 1 ~ 'MA.H'))


d <- d %>% mutate(category = ifelse(label %in% fake.diploids, 'MA.H', category))
d <- d %>% mutate(category = ifelse(label %in% fake.haploids, 'MA.D', category))

# Labeling the dates
dates <- d %>% distinct(date)
alphabet <- c('A', 'B', 'C', 'D', 'E')
dates <- dates %>% mutate(day = alphabet)
d <- d %>% left_join(dates, by='date') %>% select(-date)

head(d)
```

    ##   well treatment      slope initial     final monotone time batch label
    ## 1   A1     Blank 0.01190268 0.19525 0.2028966        0   17     2     0
    ## 2   B1        D3 0.16205402 0.22850 0.6973506        0   62     2   106
    ## 3   C1        H1 0.24591168 0.25875 0.9610287        0   31     2   101
    ## 4   D1        H2 0.02848045 0.24225 0.2147241       72    8     2   103
    ## 5   E1        H3 0.47864236 0.23025 1.0269483        9   22     2   105
    ## 6   F1        D1 0.13472608 0.24650 0.5693218       15   62     2   102
    ##   category day
    ## 1    Blank   A
    ## 2   Ctrl.H   A
    ## 3   Ctrl.H   A
    ## 4   Ctrl.H   A
    ## 5   Ctrl.H   A
    ## 6   Ctrl.D   A

``` r
d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

There seems to be something wrong with day C. So I think I want to leave
it out of the analysis.

``` r
d <- d %>% filter(day != 'C')
```

### Additional information about the days

An additional assay was conducted that compared only the control lines
with frozen lab stocks. The ancestors used for the experiment were
replicates of the frozen lab stocks. However, they demonstrated
considerable variance from the lab stocks. This was a cause for concern
which led to the conduction of the supplemental assay. The supplemental
assay demonstrated that the control lines used for the original assay
behaved similar to the lab stocks. Thus it is reasonable to add these
values to the dataset before trimming for reasonable values.

# Preliminary exploration

Plot ancestors across all days to visually check for day effects
![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Plotting the distribution of the slope values

``` r
# ci <- c(quantile(d$slope, 0.25) - 1.5* IQR(d$slope),
        # quantile(d$slope, 0.75) + 1.5* IQR(d$slope))

ci <- c(0.05, 0.22)

d %>% ggplot() +
  geom_density(aes(x=slope)) +
  geom_vline(xintercept = ci, color='red', linetype='dashed')
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
blanks <-  d %>% filter(category == 'Blank')
```

visualize the blanks

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Setting cutoff for blanks around 0.05 and using the good blanks to make
blank predictions for each day.

``` r
bad.blanks <- blanks %>% filter(slope > 0.05)
good.blanks <- blanks %>% filter(slope < 0.05)

model <- lmer(slope~(1|day), data = good.blanks)
dates.predict <- data.frame(date=distinct(d, day))
dates.predict$null = predict(model, dates.predict)
```

Checking for effects of time (time taken to attain maximum growth rate).
Each time stamp is of 15 minutes.

Plotting the low values of time

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Samples with time \< 10 were identified to be erroneous. 20 timestamps
(4 hours) was a reasonable cutoff since the samples beyond that were
more likely to be good.

## Error removal

Accounting for experimental error by removing lines according to the
following rules:  
1. Remove unreasonable slope values  
2. Remove the bad blanks.  
3. Remove lines where the initial optical density was less than the
final optical density since the optical density should not decrease
unless there was an error.  
4. Remove the diploid ancestors that were found to be haploids.  
5. Remove the lines that reach saturation within the first 2 hours which
is too soon to reach saturation.

``` r
data <- d %>%
  filter(slope > ci[1] & slope < ci[2]) %>%
  anti_join(bad.blanks) %>%
  filter(initial <= final) %>%
  filter(!(treatment %in% anc.hap.fake)) %>%
  filter(time > 10) %>%
  filter(monotone == 0)
```

    ## Joining with `by = join_by(well, treatment, slope, initial, final, monotone,
    ## time, batch, label, category, day)`

Plotting the ancestors again
![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->

## Preparing for data analysis

``` r
df <- data %>%
  select(treatment, label, slope, initial, day, category, time) %>%
  filter(category != 'Blank')

df <- df %>% rename('lineid' = 'treatment')

df$ploidy <- case_when(
  df$category %in% c('Ctrl.H', 'MA.H') ~ 'Haploid',
  df$category %in% c('Ctrl.D', 'MA.D') ~ 'Diploid'
)

df$MA <- case_when(
  df$category %in% c('MA.H', 'MA.D') ~ 'MA',
  df$category %in% c('Ctrl.H', 'Ctrl.D') ~ 'Ctrl',
)

head(df)
```

    ##   lineid label     slope initial day category time  ploidy   MA
    ## 1     D1   102 0.1694542 0.20300   A   Ctrl.D   61 Diploid Ctrl
    ## 2     H1   101 0.2036260 0.32475   A   Ctrl.H   31 Haploid Ctrl
    ## 3     D1   102 0.1647003 0.20525   A   Ctrl.D   61 Diploid Ctrl
    ## 4     D1   102 0.1456489 0.22625   A   Ctrl.D   62 Diploid Ctrl
    ## 5      1     1 0.1201696 0.21475   A     MA.H   76 Haploid   MA
    ## 6      3     3 0.1485536 0.21100   A     MA.H   75 Haploid   MA

# Data Analysis

### Calculating relative fitness

``` r
# MA lines
trt <- df %>% filter(label > 0 & label <= 100)
# ancestor lines
ctrl <- df %>% filter(label > 100)

# model to predict the slope of the ancestors of each ploidy with date as a random effect
mod <- lmer(slope ~ ploidy + (1|day), ctrl)

# making predictions for each date
ctrl.predict <- data.frame(distinct(ctrl, day, ploidy))
ctrl.predict$ctrl <- predict(mod, ctrl.predict)

# combining with the MA line dateset
trt <- left_join(trt, ctrl.predict, by=c('day', 'ploidy'))

# calculating relative fitness as the difference between the slope for the MA line and its ancestor
trt$rel.fit <- trt$slope - trt$ctrl

# grouping all the readings for the lines and summarizing by calculating the mean relative fitness
trt <- trt %>% group_by(ploidy, label) %>%
  summarize(rel.fit = mean(rel.fit, na.rm=TRUE)) %>%
  ungroup()
```

    ## `summarise()` has grouped output by 'ploidy'. You can override using the
    ## `.groups` argument.

``` r
# mean and standard deviation for the diploids
mu.dip <- mean(trt$rel.fit[trt$ploidy == 'Diploid'], na.rm=TRUE)
sd.dip <- sd(trt$rel.fit[trt$ploidy == 'Diploid'])
paste0("Diploids: ", mu.dip, " Deviation: ", sd.dip)
```

    ## [1] "Diploids: -0.0235314170025826 Deviation: 0.019874181393014"

``` r
# mean and standard deviation for the haploids
mu.hap <- mean(trt$rel.fit[trt$ploidy == 'Haploid'], na.rm=TRUE)
sd.hap <- sd(trt$rel.fit[trt$ploidy == 'Haploid'])
paste0("Haploids: ", mu.hap, " Deviation: ", sd.dip)
```

    ## [1] "Haploids: -0.0293484693270281 Deviation: 0.019874181393014"

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 50 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Black line is the 0 line. The haploid MA lines performed worse than the
diploid MA lines.

### Using lmer models to test for any differences between the ancestors

``` r
null <- lmer(slope ~ 1 + (1|day), ctrl)
full <- lmer(slope ~ ploidy + (1|day), ctrl)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ ploidy + (1 | day)
    ##    Data: ctrl
    ## 
    ## REML criterion at convergence: -1492.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8537 -0.3869  0.1957  0.6135  2.7750 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  day      (Intercept) 0.0002532 0.01591 
    ##  Residual             0.0004053 0.02013 
    ## Number of obs: 306, groups:  day, 4
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)   0.151862   0.008210  18.498
    ## ploidyHaploid 0.003644   0.003078   1.184
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## ploidyHapld -0.110

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: ctrl
    ## Models:
    ## null: slope ~ 1 + (1 | day)
    ## full: slope ~ ploidy + (1 | day)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
    ## null    3 -1503.4 -1492.2 754.68  -1509.4                     
    ## full    4 -1502.6 -1487.7 755.29  -1510.6 1.2176  1     0.2698

Ploidy does not have a significant effect on fitness in the ancestors.

### Using lmer models with lineid and date as random effects.

##### Looking for MA-ploidy interaction:

``` r
null <- lmer(slope ~ MA + ploidy + (1|day) + (1|lineid), df)
summary(null)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ MA + ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2877
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8995 -0.4936  0.1717  0.6327  2.9114 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001494 0.01222 
    ##  day      (Intercept) 0.0001733 0.01316 
    ##  Residual             0.0004394 0.02096 
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)    0.161089   0.009055  17.790
    ## MAMA          -0.032883   0.006394  -5.143
    ## ploidyHaploid -0.007141   0.003388  -2.108
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA  
    ## MAMA        -0.616       
    ## ploidyHapld -0.302  0.144

``` r
full <- lmer(slope ~ MA*ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ MA * ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2871.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8738 -0.5002  0.1751  0.6394  2.9195 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001499 0.01224 
    ##  day      (Intercept) 0.0001793 0.01339 
    ##  Residual             0.0004392 0.02096 
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error t value
    ## (Intercept)         0.150169   0.014113  10.641
    ## MAMA               -0.021667   0.012775  -1.696
    ## ploidyHaploid       0.006471   0.013792   0.469
    ## MAMA:ploidyHaploid -0.014496   0.014237  -1.018
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA   pldyHp
    ## MAMA        -0.858              
    ## ploidyHapld -0.786  0.857       
    ## MAMA:pldyHp  0.762 -0.865 -0.969

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: df
    ## Models:
    ## null: slope ~ MA + ploidy + (1 | day) + (1 | lineid)
    ## full: slope ~ MA * ploidy + (1 | day) + (1 | lineid)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
    ## null    6 -2891.2 -2864.6 1451.6  -2903.2                     
    ## full    7 -2890.2 -2859.3 1452.1  -2904.2 1.0516  1     0.3051

No significant MA-ploidy interaction

##### Exploring effect of MA on the slope

``` r
# lmer for effect of MA
null <- lmer(slope ~ 1 + ploidy + (1|day) + (1|lineid), df)
full <- lmer(slope ~ MA + ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ MA + ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2877
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8995 -0.4936  0.1717  0.6327  2.9114 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001494 0.01222 
    ##  day      (Intercept) 0.0001733 0.01316 
    ##  Residual             0.0004394 0.02096 
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)    0.161089   0.009055  17.790
    ## MAMA          -0.032883   0.006394  -5.143
    ## ploidyHaploid -0.007141   0.003388  -2.108
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA  
    ## MAMA        -0.616       
    ## ploidyHapld -0.302  0.144

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: df
    ## Models:
    ## null: slope ~ 1 + ploidy + (1 | day) + (1 | lineid)
    ## full: slope ~ MA + ploidy + (1 | day) + (1 | lineid)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)    
    ## null    5 -2870.8 -2848.7 1440.4  -2880.8                         
    ## full    6 -2891.2 -2864.6 1451.6  -2903.2 22.314  1  2.315e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The MA lines have different slopes than the ancestors

##### Exploring effect of ploidy on the slope

``` r
null <- lmer(slope ~ 1 + MA + (1|day) + (1|lineid), df)
full <- lmer(slope ~ ploidy + MA + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ ploidy + MA + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2877
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8995 -0.4936  0.1717  0.6327  2.9114 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001494 0.01222 
    ##  day      (Intercept) 0.0001733 0.01316 
    ##  Residual             0.0004394 0.02096 
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)    0.161089   0.009055  17.790
    ## ploidyHaploid -0.007141   0.003388  -2.108
    ## MAMA          -0.032883   0.006394  -5.143
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) pldyHp
    ## ploidyHapld -0.302       
    ## MAMA        -0.616  0.144

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: df
    ## Models:
    ## null: slope ~ 1 + MA + (1 | day) + (1 | lineid)
    ## full: slope ~ ploidy + MA + (1 | day) + (1 | lineid)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)  
    ## null    5 -2888.7 -2866.6 1449.3  -2898.7                       
    ## full    6 -2891.2 -2864.6 1451.6  -2903.2 4.4603  1    0.03469 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Haploids and diploids may have different slopes

##### Combining information from the mutation rate dataset

``` r
# read the files
mut <- read_delim("pombe_MA_data.txt")
```

    ## Rows: 100 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): ploidy, ploidy_final
    ## dbl (8): line, transfers, shared_ancestry, generations, chr_gens1, chr_gens2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# combine with the relative fitness data
trt2 <- trt %>%
  left_join(mut, by=c('label'='line'))
trt2 <- trt2 %>% select(ploidy.x, label, rel.fit, ploidy.y, ploidy_final,
                        n.SNM, n.indel)
trt2 <- trt2 %>% mutate(mutations = n.SNM + n.indel)
```

    ## Warning: Removed 50 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

Some diploids have a relative fitness that is greater than 0. Diploids
have more mutations than haploids.

##### Using a linear model to explore mutation rate and ploidy interaction

``` r
mod <- lm(rel.fit ~ ploidy.x*mutations + ploidy.x + mutations, trt2)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = rel.fit ~ ploidy.x * mutations + ploidy.x + mutations, 
    ##     data = trt2)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.043783 -0.013695  0.001188  0.017205  0.033375 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -2.467e-02  6.265e-03  -3.937 0.000277 ***
    ## ploidy.xHaploid           -1.849e-02  1.676e-02  -1.103 0.275593    
    ## mutations                  6.204e-05  3.014e-04   0.206 0.837817    
    ## ploidy.xHaploid:mutations  2.796e-03  2.772e-03   1.009 0.318393    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01965 on 46 degrees of freedom
    ##   (50 observations deleted due to missingness)
    ## Multiple R-squared:  0.03323,    Adjusted R-squared:  -0.02982 
    ## F-statistic: 0.527 on 3 and 46 DF,  p-value: 0.6659

``` r
qqnorm(resid(mod))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

##### Linear model with no mutation rate and ploidy interactions

``` r
mod <- lm(rel.fit ~ ploidy.x + mutations, trt2)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = rel.fit ~ ploidy.x + mutations, data = trt2)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.043673 -0.013638  0.001264  0.016992  0.033353 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     -2.527e-02  6.237e-03  -4.052 0.000189 ***
    ## ploidy.xHaploid -4.535e-03  9.458e-03  -0.479 0.633832    
    ## mutations        9.509e-05  2.996e-04   0.317 0.752390    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01965 on 47 degrees of freedom
    ##   (50 observations deleted due to missingness)
    ## Multiple R-squared:  0.01185,    Adjusted R-squared:  -0.0302 
    ## F-statistic: 0.2817 on 2 and 47 DF,  p-value: 0.7558

``` r
qqnorm(resid(mod))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

##### Linear model with only the mutations

``` r
mod <- lm(rel.fit ~ mutations, trt2)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = rel.fit ~ mutations, data = trt2)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.043025 -0.015831  0.001851  0.017043  0.033756 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.0268425  0.0052669  -5.096  5.8e-06 ***
    ## mutations    0.0001565  0.0002687   0.582    0.563    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01949 on 48 degrees of freedom
    ##   (50 observations deleted due to missingness)
    ## Multiple R-squared:  0.007013,   Adjusted R-squared:  -0.01367 
    ## F-statistic: 0.339 on 1 and 48 DF,  p-value: 0.5631

``` r
qqnorm(resid(mod))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

##### Combining the mutation dataset with the entire assay data

``` r
mut_join <- trt2 %>% select(label, mutations)
df <- df %>% left_join(mut_join) %>%
  mutate(mutations = ifelse(is.na(mutations), 0, mutations))
```

    ## Joining with `by = join_by(label)`

##### Mutation rate - ploidy interaction

``` r
null <- lmer(slope ~ mutations + ploidy + (1|day) + (1|lineid), df)
summary(null)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ mutations + ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2848.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.0852 -0.4925  0.1377  0.6316  2.8137 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 2.237e-04 0.01496 
    ##  day      (Intercept) 9.683e-05 0.00984 
    ##  Residual             4.392e-04 0.02096 
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error t value
    ## (Intercept)    0.1357176  0.0069607  19.498
    ## mutations     -0.0002321  0.0002283  -1.017
    ## ploidyHaploid -0.0073251  0.0045261  -1.618
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns
    ## mutations   -0.555       
    ## ploidyHapld -0.567  0.553

``` r
full <- lmer(slope ~ mutations*ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ mutations * ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2837.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.0687 -0.4931  0.1392  0.6213  2.8298 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 2.180e-04 0.014766
    ##  day      (Intercept) 9.998e-05 0.009999
    ##  Residual             4.396e-04 0.020967
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                           Estimate Std. Error t value
    ## (Intercept)              1.310e-01  7.748e-03  16.902
    ## mutations                3.309e-05  2.924e-04   0.113
    ## ploidyHaploid           -1.050e-04  6.737e-03  -0.016
    ## mutations:ploidyHaploid -6.594e-04  4.591e-04  -1.436
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns pldyHp
    ## mutations   -0.656              
    ## ploidyHapld -0.658  0.757       
    ## mttns:pldyH  0.432 -0.633 -0.746

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: df
    ## Models:
    ## null: slope ~ mutations + ploidy + (1 | day) + (1 | lineid)
    ## full: slope ~ mutations * ploidy + (1 | day) + (1 | lineid)
    ##      npar     AIC     BIC logLik deviance Chisq Df Pr(>Chisq)
    ## null    6 -2869.8 -2843.3 1440.9  -2881.8                    
    ## full    7 -2869.9 -2839.0 1442.0  -2883.9 2.083  1     0.1489

No significant interaction in number of mutations and ploidy

##### Effect of mutations

``` r
null <- lmer(slope ~ 1 + ploidy + (1|day) + (1|lineid), df)
full <- lmer(slope ~ mutations + ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ mutations + ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2848.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.0852 -0.4925  0.1377  0.6316  2.8137 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 2.237e-04 0.01496 
    ##  day      (Intercept) 9.683e-05 0.00984 
    ##  Residual             4.392e-04 0.02096 
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error t value
    ## (Intercept)    0.1357176  0.0069607  19.498
    ## mutations     -0.0002321  0.0002283  -1.017
    ## ploidyHaploid -0.0073251  0.0045261  -1.618
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns
    ## mutations   -0.555       
    ## ploidyHapld -0.567  0.553

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: df
    ## Models:
    ## null: slope ~ 1 + ploidy + (1 | day) + (1 | lineid)
    ## full: slope ~ mutations + ploidy + (1 | day) + (1 | lineid)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
    ## null    5 -2870.8 -2848.7 1440.4  -2880.8                     
    ## full    6 -2869.8 -2843.3 1440.9  -2881.8 1.0176  1     0.3131

No significant fitness effects of number of mutations.

###### Effect of ploidy

``` r
null <- lmer(slope ~ 1 + mutations + (1|day) + (1|lineid), df)
full <- lmer(slope ~ ploidy + mutations + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ ploidy + mutations + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -2848.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.0852 -0.4925  0.1377  0.6316  2.8137 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 2.237e-04 0.01496 
    ##  day      (Intercept) 9.683e-05 0.00984 
    ##  Residual             4.392e-04 0.02096 
    ## Number of obs: 613, groups:  lineid, 106; day, 4
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error t value
    ## (Intercept)    0.1357176  0.0069607  19.498
    ## ploidyHaploid -0.0073251  0.0045261  -1.618
    ## mutations     -0.0002321  0.0002283  -1.017
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) pldyHp
    ## ploidyHapld -0.567       
    ## mutations   -0.555  0.553

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: df
    ## Models:
    ## null: slope ~ 1 + mutations + (1 | day) + (1 | lineid)
    ## full: slope ~ ploidy + mutations + (1 | day) + (1 | lineid)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
    ## null    5 -2869.2 -2847.1 1439.6  -2879.2                     
    ## full    6 -2869.8 -2843.3 1440.9  -2881.8 2.6498  1     0.1036
