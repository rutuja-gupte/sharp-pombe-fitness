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

d.test$initial <- sapply(c(1:nrow(d.test)), function(r){
  mean(assay.testdata[,which(names(assay.testdata)==d.test$well[r])][1:4])
})

d.test$final <- sapply(c(1:nrow(d.test)), function(r){
  temp <- assay.testdata[,which(names(assay.testdata)==d.test$well[r])]
  temp <- temp[!is.na(temp)]
  n <- length(temp)
  mean(assay.testdata[,which(names(assay.testdata)==d.test$well[r])][n-3:n])
})
  
well <- assay.testdata$K12
time <- seq(1, length(well))

blanks <- d.test[d.test$treatment == "Blank",]
good.blanks <- blanks %>% filter(final - initial < 0.05)
blank.wells <- good.blanks$well
blank.val <- mean(unlist(assay.testdata[, blank.wells]), na.rm=TRUE)
# assay.testdata[, 3:ncol(assay.testdata)] <- assay.testdata[, 3:ncol(assay.testdata)] - blank.val
# assay.testdata <- assay.testdata %>% select(-blanks$well)
# d.test <- d.test %>% filter(!(well %in% blanks$well))
# assay.testdata[assay.testdata < blank.val] <- NA

ggplot() + 
  geom_point(aes(x=time, y=well)) +
  geom_point(aes(x=time, y=well-blank.val), color="blue")
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
## Now trying to smooth 
smoothed <- predict(loess(log(well) ~ time, degree=1, span=0.075), time)

derivs <- nderiv(loess(log(assay.testdata$K12) ~ time, degree=1), time)

ggplot() +
  geom_point(aes(x=time, y=log(assay.testdata$K12)))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
ggplot() + 
  geom_point(aes(x=time, y=derivs))
```

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggplot() + 
  geom_point(aes(x=time, y=smoothed), color="red") +
  geom_point(aes(x=time, y=log(well))) 
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

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

    ##   well treatment initial     final     slope batch label category day
    ## 1   B1        H1 0.17475 0.6231032 0.2883872     1   101   Ctrl.H  1A
    ## 2   C1        H2 0.17300 0.6569032 0.2966857     1   103   Ctrl.H  1A
    ## 3   D1        H3 0.17825 0.5738387 0.3083439     1   105   Ctrl.H  1A
    ## 4   E1        D1 0.17000 0.4803484 0.2835096     1   102   Ctrl.D  1A
    ## 5   F1        D2 0.17175 0.3771548 0.1335017     1   104   Ctrl.H  1A
    ## 6   G1        D3 0.16800 0.5019871 0.2989705     1   106   Ctrl.H  1A

``` r
d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

    ## Warning: Removed 581 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Who are these outliers????

``` r
# d %>% filter(slope > 1) %>% select(well, day)
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

    ## Warning: Removed 260 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Plotting the distribution of the slope values

``` r
# ci <- c(quantile(d$slope, 0.25) - 1.5* IQR(d$slope),
# quantile(d$slope, 0.75) + 1.5* IQR(d$slope))

ci <- c(0.04, 0.27)

# ci <- c(quantile(d$slope, 0.025, na.rm=TRUE), quantile(d$slope, 0.975, na.rm=TRUE))

d %>% ggplot() +
  geom_density(aes(x=slope)) +
  geom_vline(xintercept = ci, color='red', linetype='dashed')
```

    ## Warning: Removed 581 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = ci, color="red")
```

    ## Warning: Removed 581 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
blanks <-  d %>% filter(category == 'Blank')
```

Plotting the distribution of the final OD values

``` r
d %>% ggplot() +
  geom_density(aes(x=final))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

visualize the blanks
![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

Setting cutoff for blanks around 0.06 and using the good blanks to make
blank predictions for each day.

``` r
# bad.blanks <- blanks %>% filter(slope > 0.06)
# good.blanks <- blanks %>% filter(slope < 0.06)
# 
# model <- lmer(slope~(1|day), data = good.blanks)
# dates.predict <- data.frame(date=distinct(d, day))
# dates.predict$null = predict(model, dates.predict)

bad.blanks <- blanks %>% filter(final - initial >= 0.05)
good.blanks <- blanks %>% filter(final - initial < 0.05)
```

``` r
# media.od <- good.blanks %>% group_by(day) %>%
#   summarize(avg.init = mean(initial))
# 
# good.blanks %>% ggplot() +
#   geom_histogram(aes(x=initial)) +
#   facet_wrap(vars(day))
# 
# d3 <- d %>% left_join(media.od, by=c("day"))
# head(d3)
# # d3$slope.val <- d3$slope-
```

Checking for effects of time (time taken to attain maximum growth rate).
Each time stamp is of 15 minutes.

Plotting the low values of time

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
data <- d %>%
  # anti_join(bad.blanks) %>%
  filter(initial <= final) #%>%
  # filter(final >= lower) %>%
  # filter(initial <= upper) %>%
  # filter(monotone == 0) # %>%
  # filter(!(batch == 1 & treatment %in% c("D1", "D2", "D3"))) %>%
  # filter(!(batch == 2 & treatment %in% c("H1", "H2", "H3", "D2", "D3"))) %>%
  # filter(!(batch == 1 & category == "MA.D")) %>%
  # filter(!(batch == 2 & category == "MA.H"))

data %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

    ## Warning: Removed 482 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# data %>% filter(slope > 0.3) %>% select(well, day)

data %>%
  filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=treatment, y=slope, color=treatment)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

    ## Warning: Removed 167 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

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
# 
# df <- data  %>%  filter(category != 'Blank')
# # df <- data %>%
# #   select(treatment, slope, initial, time, batch, label, category, day) %>%
# #   filter(category != 'Blank')
# 
# df <- df %>% rename('lineid' = 'treatment')
# 
# df$ploidy <- case_when(
#   df$category %in% c('Ctrl.H', 'MA.H') ~ 'Haploid',
#   df$category %in% c('Ctrl.D', 'MA.D') ~ 'Diploid'
# )
# 
# df$MA <- case_when(
#   df$category %in% c('MA.H', 'MA.D') ~ 'MA',
#   df$category %in% c('Ctrl.H', 'Ctrl.D') ~ 'Ctrl',
# )
# 
# head(df)
```

# Data Analysis

### Calculating relative fitness

``` r
# # MA lines
# trt <- df %>% filter(MA == "MA")
# # ancestor lines
# ctrl <- df %>% filter(MA == "Ctrl")
# 
# # model to predict the slope of the ancestors of each ploidy with date as a random effect
# mod <- lmer(slope ~ ploidy + (1|day), ctrl)
# 
# # making predictions for each date
# ctrl.predict <- data.frame(distinct(ctrl, day, ploidy))
# ctrl.predict$ctrl <- predict(mod, ctrl.predict)
# ctrl.predict
# 
# # combining with the MA line dateset
# trt <- left_join(trt, ctrl.predict, by=c('day', 'ploidy'))
# 
# # calculating relative fitness as the difference between the slope for the MA line and its ancestor
# trt$rel.fit <- trt$slope - trt$ctrl
# 
# # grouping all the readings for the lines and summarizing by calculating the mean relative fitness
# trt <- trt %>% group_by(ploidy, label) %>%
#   summarize(rel.fit = mean(rel.fit, na.rm=TRUE)) %>%
#   ungroup()
# 
# # mean and standard deviation for the diploids
# mu.dip <- mean(trt$rel.fit[trt$ploidy == 'Diploid'], na.rm=TRUE)
# sd.dip <- sd(trt$rel.fit[trt$ploidy == 'Diploid'])
# paste0("Diploids: ", mu.dip, " Deviation: ", sd.dip)
# 
# # mean and standard deviation for the haploids
# mu.hap <- mean(trt$rel.fit[trt$ploidy == 'Haploid'], na.rm=TRUE)
# sd.hap <- sd(trt$rel.fit[trt$ploidy == 'Haploid'])
# paste0("Haploids: ", mu.hap, " Deviation: ", sd.dip)
```

Black line is the 0 line. The haploid MA lines performed worse than the
diploid MA lines.

### Using lmer models to test for any differences between the ancestors

``` r
# null <- lmer(slope ~ 1 + (1|lineid) + (1|batch), ctrl)
# full <- lmer(slope ~ ploidy + (1|lineid) + (1|batch), ctrl)
# summary(full)
# 
# mod <- anova(null, full)
# mod
# 
# t.test(ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Haploid"],
#        ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Diploid"])
```

Ploidy does not have a significant effect on fitness in the ancestors.

### Ancestor-MA comparisons

At the moment, this is the analysis that I know I can do for sure.

``` r
# hap <- df %>% filter(ploidy == "Haploid", batch == 1)
# dip <- df %>% filter(ploidy == "Diploid", batch == 2)
# 
# ## isSingular when including lineid
# 
# null <- lmer(slope ~ 1 + (1|day), hap)
# summary(null)
# full <- lmer(slope ~ MA + (1|day), hap)
# summary(full)
# 
# mod <- anova(null, full)
# mod
# 
# null <- lmer(slope ~ 1 + (1|day), dip)
# summary(null)
# full <- lmer(slope ~ MA + (1|day), dip)
# summary(full)
# 
# mod <- anova(null, full)
# mod
```

# Pause reading here for now. Everything before this needs a sanity check first.

### Using lmer models with lineid and date as random effects.

##### Looking for MA-ploidy interaction:

``` r
# null <- lmer(slope ~ MA + ploidy + (1|day) + (1|lineid), df)
# summary(null)
# full <- lmer(slope ~ MA*ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod
```

No significant MA-ploidy interaction

##### Exploring effect of MA on the slope

``` r
# # lmer for effect of MA
# null <- lmer(slope ~ 1 + ploidy + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ MA + ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod
```

The MA lines have different slopes than the ancestors

##### Exploring effect of ploidy on the slope

``` r
# null <- lmer(slope ~ 1 + MA + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ ploidy + MA + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod
```

Haploids and diploids may have different slopes

##### Combining information from the mutation rate dataset

``` r
# # read the files
# mut <- read_delim("pombe_MA_data.txt")
# 
# # combine with the relative fitness data
# trt2 <- trt %>%
#   left_join(mut, by=c('label'='line'))
# trt2 <- trt2 %>% select(ploidy.x, label, rel.fit, ploidy.y, ploidy_final,
#                         n.SNM, n.indel)
# trt2 <- trt2 %>% mutate(mutations = n.SNM + n.indel)
```

Some diploids have a relative fitness that is greater than 0. Diploids
have more mutations than haploids.

##### Using a linear model to explore mutation rate and ploidy interaction

``` r
# mod <- lm(rel.fit ~ ploidy.x*mutations + ploidy.x + mutations, trt2)
# summary(mod)
# qqnorm(resid(mod))
```

##### Linear model with no mutation rate and ploidy interactions

``` r
# mod <- lm(rel.fit ~ ploidy.x + mutations, trt2)
# summary(mod)
# qqnorm(resid(mod))
```

##### Linear model with only the mutations

``` r
# mod <- lm(rel.fit ~ mutations, trt2)
# summary(mod)
# qqnorm(resid(mod))
```

##### Combining the mutation dataset with the entire assay data

``` r
# mut_join <- trt2 %>% select(label, mutations)
# df <- df %>% left_join(mut_join) %>%
#   mutate(mutations = ifelse(is.na(mutations), 0, mutations))
```

##### Mutation rate - ploidy interaction

``` r
# null <- lmer(slope ~ mutations + ploidy + (1|day) + (1|lineid), df)
# summary(null)
# full <- lmer(slope ~ mutations*ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod
```

No significant interaction in number of mutations and ploidy

##### Effect of mutations

``` r
# null <- lmer(slope ~ 1 + ploidy + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ mutations + ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod
```

No significant fitness effects of number of mutations.

###### Effect of ploidy

``` r
# null <- lmer(slope ~ 1 + mutations + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ ploidy + mutations + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod
```
