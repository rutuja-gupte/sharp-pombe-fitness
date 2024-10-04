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

## Figuring out the parameters for the splines (I am just playing around)

``` r
d.test <- read.csv("data/Ref 04 10.csv")
assay.testdata <- read.delim("data/Rutuja 04 10.txt")
```

``` r
well <- assay.testdata$M22
time <- seq(1, length(well))


ggplot() + 
  geom_point(aes(x=time, y=well))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
## Now trying to smooth 
smoothed <- predict(loess(log(well) ~ time, degree=1, span=0.075), time)

ggplot() + 
  geom_point(aes(x=time, y=smoothed), color="red") +
  geom_point(aes(x=time, y=log(well)))
```

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
dates.3 <- c(mdy(11022023))
# dates.2 <- c(dates.2, mdy(11022023))

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
that changed ploidy. For now I am excluding the lines that changed
ploidy. Originally I was categorizing them as whatever their final
ploidy was but I need to rethink that.

``` r
d$batch[d$date %in% dates.1] <- 1
d$batch[d$date %in% dates.2] <- 2
d$batch[d$date %in% dates.3] <- 3


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


# 
# d <- d %>% mutate(category = ifelse(label %in% fake.diploids, 'MA.H', category))
# d <- d %>% mutate(category = ifelse(label %in% fake.haploids, 'MA.D', category))

# d <- d %>% filter(!(label %in% fake.diploids))
# d <- d %>% filter(!(label %in% fake.haploids))

# Labeling the dates
dates <- d %>% distinct(date)
alphabet <- c('1A', '1B', '1C', '1D', '2A', '2B', '2C', '2D', '3A')
dates <- dates %>% mutate(day = alphabet)
d <- d %>% left_join(dates, by='date') %>% select(-date)

head(d)
```

    ##   well treatment       slope initial     final monotone final_slope time batch
    ## 1   A1     Blank 0.006873357 0.16275 0.1664968        7  0.00146783   81     1
    ## 2   B1        H1 0.174031549 0.17475 0.6231032        0  0.04396843   57     1
    ## 3   C1        H2 0.180485167 0.17300 0.6569032        0  0.04575788   53     1
    ## 4   D1        H3 0.177495892 0.17825 0.5738387        0  0.04193429   57     1
    ## 5   E1        D1 0.160080758 0.17000 0.4803484        0  0.03960327   67     1
    ## 6   F1        D2 0.126360761 0.17175 0.3771548        2  0.03524252   67     1
    ##   label category day
    ## 1     0    Blank  1A
    ## 2   101   Ctrl.H  1A
    ## 3   103   Ctrl.H  1A
    ## 4   105   Ctrl.H  1A
    ## 5   102   Ctrl.D  1A
    ## 6   104   Ctrl.H  1A

``` r
d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Who are these outliers????

``` r
d %>% filter(slope > 1) %>% select(well, day)
```

    ## [1] well day 
    ## <0 rows> (or 0-length row.names)

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

# ci <- c(0.04, 0.27)

ci <- c(quantile(d$slope, 0.025), quantile(d$slope, 0.975))

d %>% ggplot() +
  geom_density(aes(x=slope)) +
  geom_vline(xintercept = ci, color='red', linetype='dashed')
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = ci, color="red")
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
blanks <-  d %>% filter(category == 'Blank')
```

Plotting the distribution of the final OD values

``` r
d %>% ggplot() +
  geom_density(aes(x=final))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

visualize the blanks

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Setting cutoff for blanks around 0.06 and using the good blanks to make
blank predictions for each day.

``` r
bad.blanks <- blanks %>% filter(slope > 0.06)
good.blanks <- blanks %>% filter(slope < 0.06)

model <- lmer(slope~(1|day), data = good.blanks)
dates.predict <- data.frame(date=distinct(d, day))
dates.predict$null = predict(model, dates.predict)
```

Checking for effects of time (time taken to attain maximum growth rate).
Each time stamp is of 15 minutes.

Plotting the low values of time

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Samples with time \< 10 were identified to be erroneous.

Looking to see if they reached saturation or not

``` r
# sat <- c(0.2, 1.1)
daily_cutoff <- d %>% group_by(day) %>%
  summarize(lower = quantile(final, 0.05))

d2 <- d %>% left_join(daily_cutoff, by=c("day"))

d2 %>% ggplot() +
  geom_density(aes(x=final)) +
  geom_vline(aes(xintercept=lower), color='red', linetype='dashed') +
  facet_wrap(vars(day))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
daily_cutoff <- d %>% group_by(day) %>%
  summarize(upper = quantile(initial, 0.95))

d2 <- d2 %>% left_join(daily_cutoff, by=c("day"))

d2 %>% ggplot() +
  geom_density(aes(x=initial)) +
  geom_vline(aes(xintercept=upper), color='red', linetype='dashed') +
  facet_wrap(vars(day))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

## Error removal

Accounting for experimental error by removing lines according to the
following rules:  
Rules need to be updated

1.  Remove unreasonable slope values  
2.  Remove the bad blanks.  
3.  Remove lines where the initial optical density was more than the
    final optical density since the optical density should not decrease
    unless there was an error.  
4.  Remove the diploid ancestors that were found to be haploids.  
5.  Remove the lines that reach saturation within the first 2 hours
    which is too soon to reach saturation.

These were not the current rules

``` r
data <- d2 %>%
  anti_join(bad.blanks) %>%
  filter(initial <= final) %>%
  filter(final >= lower) %>%
  filter(initial <= upper) %>%
  filter(monotone == 0) %>%
  filter(!(batch == 1 & treatment %in% c("D1", "D2", "D3"))) %>%
  filter(!(batch == 2 & treatment %in% c("H1", "H2", "H3", "D2", "D3"))) %>%
  filter(!(batch == 1 & category == "MA.D")) %>%
  filter(!(batch == 2 & category == "MA.H"))
```

    ## Joining with `by = join_by(well, treatment, slope, initial, final, monotone,
    ## final_slope, time, batch, label, category, day)`

``` r
data %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
data %>% filter(slope > 0.3) %>% select(well, day)
```

    ##   well day
    ## 1  K15  1A
    ## 2  M15  1A
    ## 3  M10  1B
    ## 4  I17  1B

``` r
data %>%
  filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=treatment, y=slope, color=treatment)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
cat(data$well[data$day=="1A"], file = "1A.txt")
cat(data$well[data$day=="1B"], file = "1B.txt")
cat(data$well[data$day=="1C"], file = "1C.txt")
cat(data$well[data$day=="1D"], file = "1D.txt")
cat(data$well[data$day=="2A"], file = "2A.txt")
cat(data$well[data$day=="2B"], file = "2B.txt")
cat(data$well[data$day=="2C"], file = "2C.txt")
cat(data$well[data$day=="2D"], file = "2D.txt")
cat(data$well[data$day=="3A"], file = "3A.txt")
```

## Preparing for data analysis

``` r
df <- data %>%
  select(treatment, slope, initial, time, batch, label, category, day) %>%
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

    ##   lineid     slope initial time batch label category day  ploidy   MA
    ## 1     H1 0.1740315 0.17475   57     1   101   Ctrl.H  1A Haploid Ctrl
    ## 2     H2 0.1804852 0.17300   53     1   103   Ctrl.H  1A Haploid Ctrl
    ## 3     H3 0.1774959 0.17825   57     1   105   Ctrl.H  1A Haploid Ctrl
    ## 4     H1 0.1261628 0.17575   57     1   101   Ctrl.H  1A Haploid Ctrl
    ## 5     H2 0.1369487 0.17250   60     1   103   Ctrl.H  1A Haploid Ctrl
    ## 6     H3 0.1497487 0.17400   55     1   105   Ctrl.H  1A Haploid Ctrl

# Data Analysis

### Calculating relative fitness

``` r
# MA lines
trt <- df %>% filter(MA == "MA")
# ancestor lines
ctrl <- df %>% filter(MA == "Ctrl")

# model to predict the slope of the ancestors of each ploidy with date as a random effect
mod <- lmer(slope ~ ploidy + (1|day), ctrl)

# making predictions for each date
ctrl.predict <- data.frame(distinct(ctrl, day, ploidy))
ctrl.predict$ctrl <- predict(mod, ctrl.predict)
ctrl.predict
```

    ##    day  ploidy      ctrl
    ## 1   1A Haploid 0.1601441
    ## 2   1B Haploid 0.1280030
    ## 3   1C Haploid 0.1738545
    ## 4   1D Haploid 0.1667350
    ## 5   2A Diploid 0.1311803
    ## 6   2B Diploid 0.1541905
    ## 7   2C Diploid 0.1639665
    ## 8   2D Diploid 0.1610390
    ## 9   3A Haploid 0.1272664
    ## 10  3A Diploid 0.1317625

``` r
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

    ## [1] "Diploids: -0.0187457138843519 Deviation: 0.0146348976314018"

``` r
# mean and standard deviation for the haploids
mu.hap <- mean(trt$rel.fit[trt$ploidy == 'Haploid'], na.rm=TRUE)
sd.hap <- sd(trt$rel.fit[trt$ploidy == 'Haploid'])
paste0("Haploids: ", mu.hap, " Deviation: ", sd.dip)
```

    ## [1] "Haploids: 0.0102180494669903 Deviation: 0.0146348976314018"

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Black line is the 0 line. The haploid MA lines performed worse than the
diploid MA lines.

### Using lmer models to test for any differences between the ancestors

``` r
null <- lmer(slope ~ 1 + (1|lineid) + (1|batch), ctrl)
full <- lmer(slope ~ ploidy + (1|lineid) + (1|batch), ctrl)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ ploidy + (1 | lineid) + (1 | batch)
    ##    Data: ctrl
    ## 
    ## REML criterion at convergence: -2365
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2571 -0.4819  0.2288  0.6257  2.4796 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 1.329e-05 0.003646
    ##  batch    (Intercept) 2.537e-04 0.015929
    ##  Residual             6.211e-04 0.024922
    ## Number of obs: 527, groups:  lineid, 6; batch, 3
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)    0.149765   0.010262  14.594
    ## ploidyHaploid -0.004459   0.005674  -0.786
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## ploidyHapld -0.397

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: ctrl
    ## Models:
    ## null: slope ~ 1 + (1 | lineid) + (1 | batch)
    ## full: slope ~ ploidy + (1 | lineid) + (1 | batch)
    ##      npar     AIC     BIC logLik deviance Chisq Df Pr(>Chisq)
    ## null    4 -2372.6 -2355.5 1190.3  -2380.6                    
    ## full    5 -2371.3 -2350.0 1190.7  -2381.3 0.783  1     0.3762

``` r
t.test(ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Haploid"],
       ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Diploid"])
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Haploid"] and ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Diploid"]
    ## t = -1.626, df = 91.016, p-value = 0.1074
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.010956863  0.001092968
    ## sample estimates:
    ## mean of x mean of y 
    ## 0.1270407 0.1319727

Ploidy does not have a significant effect on fitness in the ancestors.

### Ancestor-MA comparisons

At the moment, this is the analysis that I know I can do for sure.

``` r
hap <- df %>% filter(ploidy == "Haploid", batch == 1)
dip <- df %>% filter(ploidy == "Diploid", batch == 2)

## isSingular when including lineid

null <- lmer(slope ~ 1 + (1|day), hap)
summary(null)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ 1 + (1 | day)
    ##    Data: hap
    ## 
    ## REML criterion at convergence: -1880.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7726 -0.5113  0.0385  0.4516  9.6216 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  day      (Intercept) 0.0001921 0.01386 
    ##  Residual             0.0010814 0.03288 
    ## Number of obs: 476, groups:  day, 4
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   0.1636     0.0071   23.04

``` r
full <- lmer(slope ~ MA + (1|day), hap)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ MA + (1 | day)
    ##    Data: hap
    ## 
    ## REML criterion at convergence: -1883.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6184 -0.5449  0.0558  0.5105  9.5967 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  day      (Intercept) 0.0001868 0.01367 
    ##  Residual             0.0010547 0.03248 
    ## Number of obs: 476, groups:  day, 4
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.157478   0.007202   21.87
    ## MAMA        0.010880   0.003014    3.61
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## MAMA -0.234

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: hap
    ## Models:
    ## null: slope ~ 1 + (1 | day)
    ## full: slope ~ MA + (1 | day)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)    
    ## null    3 -1882.9 -1870.4 944.44  -1888.9                         
    ## full    4 -1893.8 -1877.1 950.90  -1901.8 12.909  1  0.0003271 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
null <- lmer(slope ~ 1 + (1|day), dip)
summary(null)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ 1 + (1 | day)
    ##    Data: dip
    ## 
    ## REML criterion at convergence: -1549.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.7029 -0.4939  0.1962  0.6342  2.2952 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  day      (Intercept) 0.0002216 0.01489 
    ##  Residual             0.0005470 0.02339 
    ## Number of obs: 336, groups:  day, 4
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.138945   0.007596   18.29

``` r
full <- lmer(slope ~ MA + (1|day), dip)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ MA + (1 | day)
    ##    Data: dip
    ## 
    ## REML criterion at convergence: -1589.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.7523 -0.4915  0.1571  0.6488  2.4279 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  day      (Intercept) 0.0003210 0.01792 
    ##  Residual             0.0004706 0.02169 
    ## Number of obs: 336, groups:  day, 4
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.153142   0.009276  16.509
    ## MAMA        -0.020448   0.002784  -7.346
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## MAMA -0.210

``` r
mod <- anova(null, full)
```

    ## refitting model(s) with ML (instead of REML)

``` r
mod
```

    ## Data: dip
    ## Models:
    ## null: slope ~ 1 + (1 | day)
    ## full: slope ~ MA + (1 | day)
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)    
    ## null    3 -1551.7 -1540.2 778.84  -1557.7                         
    ## full    4 -1599.2 -1583.9 803.59  -1607.2 49.505  1  1.978e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Pause reading here for now. Everything before this needs a sanity check first.

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
    ## REML criterion at convergence: -4521.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3143 -0.5020  0.1118  0.5224 10.9363 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001159 0.01077 
    ##  day      (Intercept) 0.0002605 0.01614 
    ##  Residual             0.0006853 0.02618 
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)   0.127594   0.008232  15.500
    ## MAMA          0.007503   0.005171   1.451
    ## ploidyHaploid 0.026237   0.005890   4.455
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA  
    ## MAMA        -0.624       
    ## ploidyHapld -0.491  0.185

``` r
full <- lmer(slope ~ MA*ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ MA * ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -4522.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2825 -0.5134  0.1042  0.5216 11.0291 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0000860 0.009273
    ##  day      (Intercept) 0.0002600 0.016123
    ##  Residual             0.0006895 0.026259
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error t value
    ## (Intercept)         0.152802   0.011080  13.790
    ## MAMA               -0.020121   0.009936  -2.025
    ## ploidyHaploid      -0.003758   0.010942  -0.343
    ## MAMA:ploidyHaploid  0.034595   0.011186   3.093
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA   pldyHp
    ## MAMA        -0.817              
    ## ploidyHapld -0.786  0.800       
    ## MAMA:pldyHp  0.707 -0.886 -0.851

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
    ## null    6 -4534.9 -4505.2 2273.5  -4546.9                        
    ## full    7 -4541.8 -4507.1 2277.9  -4555.8 8.8539  1   0.002925 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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
    ## REML criterion at convergence: -4521.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3143 -0.5020  0.1118  0.5224 10.9363 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001159 0.01077 
    ##  day      (Intercept) 0.0002605 0.01614 
    ##  Residual             0.0006853 0.02618 
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)   0.127594   0.008232  15.500
    ## MAMA          0.007503   0.005171   1.451
    ## ploidyHaploid 0.026237   0.005890   4.455
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA  
    ## MAMA        -0.624       
    ## ploidyHapld -0.491  0.185

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
    ## null    5 -4534.8 -4510.0 2272.4  -4544.8                     
    ## full    6 -4534.9 -4505.2 2273.5  -4546.9 2.1596  1     0.1417

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
    ## REML criterion at convergence: -4521.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3143 -0.5020  0.1118  0.5224 10.9363 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001159 0.01077 
    ##  day      (Intercept) 0.0002605 0.01614 
    ##  Residual             0.0006853 0.02618 
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)   0.127594   0.008232  15.500
    ## ploidyHaploid 0.026237   0.005890   4.455
    ## MAMA          0.007503   0.005171   1.451
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) pldyHp
    ## ploidyHapld -0.491       
    ## MAMA        -0.624  0.185

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
    ## null    5 -4517.8 -4493.0 2263.9  -4527.8                         
    ## full    6 -4534.9 -4505.2 2273.5  -4546.9 19.135  1  1.218e-05 ***
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

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

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
    ## -0.033949 -0.009221 -0.001507  0.010721  0.052226 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -0.0219748  0.0051283  -4.285 4.38e-05 ***
    ## ploidy.xHaploid            0.0377041  0.0069987   5.387 5.18e-07 ***
    ## mutations                  0.0001668  0.0002390   0.698    0.487    
    ## ploidy.xHaploid:mutations -0.0013617  0.0009429  -1.444    0.152    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01563 on 95 degrees of freedom
    ## Multiple R-squared:  0.4786, Adjusted R-squared:  0.4621 
    ## F-statistic: 29.06 on 3 and 95 DF,  p-value: 2.022e-13

``` r
qqnorm(resid(mod))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

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
    ## -0.034243 -0.009714 -0.001278  0.010378  0.049183 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     -2.028e-02  5.020e-03  -4.040 0.000108 ***
    ## ploidy.xHaploid  3.013e-02  4.663e-03   6.462 4.26e-09 ***
    ## mutations        7.928e-05  2.325e-04   0.341 0.733882    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01572 on 96 degrees of freedom
    ## Multiple R-squared:  0.4671, Adjusted R-squared:  0.456 
    ## F-statistic: 42.08 on 2 and 96 DF,  p-value: 7.552e-14

``` r
qqnorm(resid(mod))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

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
    ## -0.044804 -0.010467 -0.000912  0.011218  0.058809 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.0079620  0.0029449   2.704   0.0081 ** 
    ## mutations   -0.0010258  0.0001878  -5.463 3.61e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01873 on 97 degrees of freedom
    ## Multiple R-squared:  0.2353, Adjusted R-squared:  0.2274 
    ## F-statistic: 29.85 on 1 and 97 DF,  p-value: 3.61e-07

``` r
qqnorm(resid(mod))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

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
    ## REML criterion at convergence: -4512.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3318 -0.5010  0.1122  0.5151 10.9401 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001191 0.01091 
    ##  day      (Intercept) 0.0002697 0.01642 
    ##  Residual             0.0006856 0.02618 
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)   1.341e-01  7.746e-03  17.306
    ## mutations     5.085e-05  2.241e-04   0.227
    ## ploidyHaploid 2.549e-02  6.682e-03   3.814
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns
    ## mutations   -0.539       
    ## ploidyHapld -0.620  0.491

``` r
full <- lmer(slope ~ mutations*ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ mutations * ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -4500.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3265 -0.4979  0.1104  0.5162 10.9288 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001193 0.01092 
    ##  day      (Intercept) 0.0002677 0.01636 
    ##  Residual             0.0006861 0.02619 
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)             1.346e-01  7.811e-03  17.233
    ## mutations               1.518e-05  2.353e-04   0.064
    ## ploidyHaploid           2.372e-02  7.570e-03   3.133
    ## mutations:ploidyHaploid 3.868e-04  7.733e-04   0.500
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns pldyHp
    ## mutations   -0.553              
    ## ploidyHapld -0.609  0.556       
    ## mttns:pldyH  0.141 -0.303 -0.470

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
    ##      npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
    ## null    6 -4532.8 -4503.1 2272.4  -4544.8                     
    ## full    7 -4531.1 -4496.4 2272.6  -4545.1 0.2795  1      0.597

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
    ## REML criterion at convergence: -4512.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3318 -0.5010  0.1122  0.5151 10.9401 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001191 0.01091 
    ##  day      (Intercept) 0.0002697 0.01642 
    ##  Residual             0.0006856 0.02618 
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)   1.341e-01  7.746e-03  17.306
    ## mutations     5.085e-05  2.241e-04   0.227
    ## ploidyHaploid 2.549e-02  6.682e-03   3.814
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns
    ## mutations   -0.539       
    ## ploidyHapld -0.620  0.491

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
    ## null    5 -4534.8 -4510.0 2272.4  -4544.8                     
    ## full    6 -4532.8 -4503.1 2272.4  -4544.8 0.0512  1      0.821

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
    ## REML criterion at convergence: -4512.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3318 -0.5010  0.1122  0.5151 10.9401 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001191 0.01091 
    ##  day      (Intercept) 0.0002697 0.01642 
    ##  Residual             0.0006856 0.02618 
    ## Number of obs: 1046, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)   1.341e-01  7.746e-03  17.306
    ## ploidyHaploid 2.549e-02  6.682e-03   3.814
    ## mutations     5.085e-05  2.241e-04   0.227
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) pldyHp
    ## ploidyHapld -0.620       
    ## mutations   -0.539  0.491

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
    ## null    5 -4520.4 -4495.7 2265.2  -4530.4                         
    ## full    6 -4532.8 -4503.1 2272.4  -4544.8 14.383  1  0.0001491 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
