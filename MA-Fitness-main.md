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
d.test <- read.csv("data/Ref 04 06.csv")
assay.testdata <- read.delim("data/Rutuja 04 06.txt")
```

``` r
well <- assay.testdata$A23
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
    ## 1   A1     Blank 0.008556741 0.16275 0.1664968        1  0.00143824   22     1
    ## 2   B1        H1 0.177044571 0.17475 0.6231032        0  0.04391984   56     1
    ## 3   C1        H2 0.181338943 0.17300 0.6569032        0  0.04570547   52     1
    ## 4   D1        H3 0.179039043 0.17825 0.5738387        0  0.04193645   56     1
    ## 5   E1        D1 0.163659261 0.17000 0.4803484        0  0.03961091   67     1
    ## 6   F1        D2 0.132832780 0.17175 0.3771548        1  0.03524914   67     1
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
![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Plotting the distribution of the slope values

``` r
# ci <- c(quantile(d$slope, 0.25) - 1.5* IQR(d$slope),
# quantile(d$slope, 0.75) + 1.5* IQR(d$slope))

ci <- c(0.04, 0.27)

d %>% ggplot() +
  geom_density(aes(x=slope)) +
  geom_vline(xintercept = ci, color='red', linetype='dashed')
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = ci, color="red")
```

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

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

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

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Samples with time \< 10 were identified to be erroneous.

Looking to see if they reached saturation or not

``` r
# sat <- c(0.2, 1.1)
daily_cutoff <- d %>% group_by(day) %>%
  summarize(lower = quantile(final, 0.05))

d <- d %>% left_join(daily_cutoff, by=c("day"))

d %>% ggplot() +
  geom_density(aes(x=final)) +
  geom_vline(aes(xintercept=lower), color='red', linetype='dashed') +
  # geom_vline(aes(xintercept=upper), color='red', linetype='dashed') +
  facet_wrap(vars(day))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## Error removal

Accounting for experimental error by removing lines according to the
following rules:  
1. Remove unreasonable slope values  
2. Remove the bad blanks.  
3. Remove lines where the initial optical density was more than the
final optical density since the optical density should not decrease
unless there was an error.  
4. Remove the diploid ancestors that were found to be haploids.  
5. Remove the lines that reach saturation within the first 2 hours which
is too soon to reach saturation.

``` r
data <- d %>%
  anti_join(bad.blanks) %>%
  filter(initial <= final) %>%
  filter(final >= lower) %>%
  filter(monotone == 0) %>%
  filter(!(batch == 1 & treatment %in% c("D1", "D2", "D3"))) %>%
  filter(!(batch == 2 & treatment %in% ancestors.haploid)) %>%
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

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
data %>%
  filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=treatment, y=slope, color=treatment)) + facet_grid(cols=vars(day))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
cat(data$well[data$day=="1A"])
```

    ## B1 C1 D1 H1 I1 J1 N1 E2 F2 G2 K2 L2 M2 B3 F3 G3 H3 L3 M3 N3 O3 B4 C4 E4 G4 M4 B5 C5 E5 G5 I5 K5 M5 P5 C6 E6 G6 I6 K6 M6 P6 A7 C7 E7 G7 I7 K7 M7 O7 C8 E8 G8 I8 K8 M8 O8 B9 C9 E9 G9 I9 K9 M9 O9 B10 C10 E10 G10 I10 K10 M10 A11 B11 C11 G11 K11 M11 P11 A12 C12 E12 G12 H12 I12 J12 K12 L12 M12 N12 K13 M13 O13 C14 E14 M14 O14 B15 C15 E15 G15 I15 K15 M15 O15 B16 C16 E16 G16 I16 K16 M16 A17 B17 C17 E17 G17 I17 K17 M17 P17 A18 C18 E18 G18 I18 K18 M18 P18 A19 C19 E19 G19 I19 K19 M19 O19 P19 C20 E20 G20 I20 K20 O20 C21 G21 M21 O21 D22 I22 J22 K22 A23 B23 C23 G23 H23 D24 E24 F24 J24 K24 L24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="1B"])
```

    ## C1 D1 E1 A2 B2 L2 B3 D3 O3 B4 I4 M4 A5 C5 E5 G5 I5 A6 C6 E6 G6 I6 K6 M6 A7 C7 E7 G7 I7 K7 M7 O7 P7 C8 G8 I8 K8 O8 B9 C9 E9 G9 K9 O9 C10 E10 G10 K10 M10 A11 C11 E11 K11 M11 A12 C12 E12 N12 C13 O13 C14 E14 O14 B15 C15 E15 O15 B16 C16 E16 A17 B17 I17 K17 M17 A18 A19 E19 K19 P19 E20 G20 O20 B21 B22 F22 G22 B23 D23 E23 K23 F24 G24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="1C"])
```

    ## C1 D1 E1 I1 J1 K1 O1 A2 B2 F2 G2 H2 L2 M2 N2 B3 C3 D3 E3 I3 J3 K3 O3 B4 C4 E4 G4 K4 M4 A5 B5 C5 E5 G5 I5 K5 M5 C6 E6 G6 I6 K6 M6 P6 A7 C7 G7 I7 K7 M7 O7 P7 C8 E8 G8 I8 K8 M8 O8 P8 B9 E9 G9 I9 K9 M9 O9 B10 G10 I10 K10 A11 B11 G11 M11 A12 C12 G12 H12 I12 M12 N12 P12 A13 C13 G13 H13 K13 O13 P13 C14 E14 G14 K14 M14 O14 P14 B15 C15 E15 G15 I15 K15 M15 O15 B16 C16 E16 G16 I16 K16 M16 A17 B17 C17 G17 I17 K17 M17 A18 C18 E18 G18 K18 M18 P18 A19 C19 E19 G19 K19 O19 P19 C20 E20 G20 K20 M20 O20 P20 B21 C21 G21 K21 M21 O21 B22 F22 G22 H22 L22 M22 N22 O22 B23 C23 D23 E23 I23 J23 K23 B24 F24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="1D"])
```

    ## D1 E1 K1 O1 B2 H2 M2 B3 C3 E3 I3 J3 K3 O3 B4 G4 I4 K4 M4 A5 B5 E5 G5 I5 K5 M5 A6 C6 E6 G6 I6 K6 M6 P6 A7 C7 E7 G7 I7 M7 P7 C8 G8 I8 K8 M8 B9 C9 E9 G9 I9 K9 M9 O9 B10 C10 E10 G10 I10 K10 M10 A11 B11 C11 G11 K11 M11 C12 G12 H12 I12 C13 G13 H13 I13 K13 M13 O13 C14 E14 G14 I14 M14 B15 C15 E15 G15 I15 K15 M15 O15 B16 C16 E16 G16 I16 K16 M16 A17 B17 C17 E17 G17 I17 M17 A18 C18 E18 G18 I18 K18 M18 A19 C19 E19 G19 I19 M19 O19 C20 E20 G20 O20 B21 G21 O21 B22 F22 G22 H22 L22 M22 O22 B23 C23 D23 J23 K23 F24 M24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="2A"])
```

    ## I2 F3 P3 F4 H4 O4 F5 H5 J5 L5 B6 F6 H6 L6 N6 F7 H7 J7 A8 L8 N8 F9 N9 P9 H10 N10 N11 B12 D13 J13 L14 N14 P15 L16 F17 J17 L17 D19 F19 F20 C22 I22 F23 L23 C24 I24 P24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="2B"])
```

    ## F1 L1 O2 F3 L3 P3 D4 F4 H4 J4 L4 N4 O4 D5 F5 J5 N5 B6 D6 F6 H6 J6 L6 N6 F7 H7 J7 L7 N7 A8 D8 H8 J8 L8 N8 D9 F9 H9 N9 P9 D10 F10 H10 J10 L10 N10 O10 D11 F11 H11 J11 L11 N11 D12 F12 J12 D13 J13 L13 N13 A14 F14 L14 D15 F15 H15 J15 L15 N15 P15 D16 H16 J16 L16 N16 O16 F17 H17 J17 L17 B18 F18 J18 L18 N18 A19 F19 H19 L19 N19 F20 H20 J20 L20 N20 D21 F21 H21 J21 L21 N21 P21 C22 I22 F23 O23 I24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="2C"])
```

    ## F1 L1 C2 O2 F3 L3 D4 F4 H4 J4 L4 N4 O4 D5 F5 H5 J5 L5 N5 D6 F6 H6 J6 L6 N6 D7 F7 H7 J7 L7 N7 A8 D8 F8 H8 J8 L8 N8 D9 F9 H9 J9 L9 N9 P9 D10 F10 H10 J10 L10 N10 O10 D11 F11 H11 J11 L11 N11 B12 D12 F12 J12 J13 L13 N13 D14 F14 H14 J14 L14 N14 D15 F15 H15 J15 L15 N15 P15 D16 F16 H16 J16 L16 N16 O16 D17 F17 H17 L17 D18 F18 H18 J18 L18 N18 D19 F19 H19 N19 D20 F20 H20 J20 L20 N20 H21 J21 N21 P21 F23 L23 O23 I24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="2D"])
```

    ## F1 L1 C2 I2 O2 F3 L3 P3 D4 F4 H4 L4 N4 O4 D5 F5 H5 J5 L5 N5 B6 D6 F6 H6 J6 L6 N6 D7 F7 H7 J7 L7 N7 A8 H8 L8 D9 F9 H9 J9 L9 N9 P9 D10 F10 H10 J10 L10 N10 O10 D11 F11 H11 J11 L11 N11 B12 D12 F12 J13 L13 N13 A14 D14 F14 H14 J14 N14 D15 F15 H15 J15 L15 N15 P15 F16 H16 J16 L16 N16 O16 D17 H17 J17 N17 H18 J18 L18 D19 F19 H19 J19 L19 N19 A20 F20 H20 J20 L20 N20 F21 H21 J21 L21 N21 P21 C22 I22 F23 L23 O23 C24 I24

``` r
print("")
```

    ## [1] ""

``` r
cat(data$well[data$day=="3A"])
```

    ## B1 C1 D1 E1 F1 G1 H1 I1 K1 L1 M1 N1 O1 B2 C2 E2 F2 G2 H2 I2 L2 M2 N2 O2 B3 C3 D3 E3 F3 H3 K3 L3 M3 N3 O3 P3 A4 B4 C4 D4 E4 F4 I4 J4 K4 L4 M4 N4 O4 P4 A5 B5 C5 D5 E5 F5 G5 H5 I5 J5 K5 L5 M5 N5 O5 P5 A6 B6 C6 K6 M6 N6 O6 P6 A7 B7 C7 N7 O7 P7 A8 B8 C8 N8 O8 P8 D9 E9 F9 G9 H9 I9 J9 K9 L9 P9 A10 B10 C10 D10 E10 F10 G10 H10 I10 J10 K10 L10 M10 N10 O10 P10 A11 B11 C11 D11 E11 F11 G11 I11 J11 K11 M11 N11 O11 P11 B12 C12 D12 E12 F12 G12 H12 I12 J12 K12 L12 M12 N12 O12 P12 B13 C13 D13 E13 F13 G13 H13 I13 J13 K13 L13 M13 N13 O13 P13 B14 C14 D14 E14 F14 G14 H14 I14 J14 L14 N14 O14 P14 A15 B15 C15 D15 E15 F15 K15 L15 M15 N15 O15 A16 B16 C16 N16 O16 P16 B17 C17 E17 H17 M17 N17 O17 P17 A18 B18 C18 N18 O18 P18 A19 B19 C19 D19 E19 F19 G19 H19 I19 J19 K19 L19 M19 O19 P19 A20 B20 C20 D20 E20 F20 G20 I20 J20 K20 L20 M20 N20 O20 P20 A22 B22 C22 D22 E22 F22 H22 I22 K22 O22 A23 B23 C23 D23 E23 F23 H23 I23 M23 O23 H24 J24

``` r
print("")
```

    ## [1] ""

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
    ## 1     H1   101 0.1770446 0.17475  1A   Ctrl.H   56 Haploid Ctrl
    ## 2     H2   103 0.1813389 0.17300  1A   Ctrl.H   52 Haploid Ctrl
    ## 3     H3   105 0.1790390 0.17825  1A   Ctrl.H   56 Haploid Ctrl
    ## 4     H1   101 0.1262669 0.17575  1A   Ctrl.H   55 Haploid Ctrl
    ## 5     H2   103 0.1403413 0.17250  1A   Ctrl.H   60 Haploid Ctrl
    ## 6     H3   105 0.1540794 0.17400  1A   Ctrl.H   54 Haploid Ctrl

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
ctrl.predict
```

    ##    day  ploidy      ctrl
    ## 1   1A Haploid 0.1654282
    ## 2   1B Haploid 0.1374850
    ## 3   1C Haploid 0.1856451
    ## 4   1D Haploid 0.1729962
    ## 5   2A Diploid 0.1345879
    ## 6   2B Diploid 0.1579472
    ## 7   2C Diploid 0.1661826
    ## 8   2D Diploid 0.1633252
    ## 9   3A Haploid 0.1293319
    ## 10  3A Diploid 0.1354217

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

    ## [1] "Diploids: -0.0205261449139338 Deviation: 0.0153652158371832"

``` r
# mean and standard deviation for the haploids
mu.hap <- mean(trt$rel.fit[trt$ploidy == 'Haploid'], na.rm=TRUE)
sd.hap <- sd(trt$rel.fit[trt$ploidy == 'Haploid'])
paste0("Haploids: ", mu.hap, " Deviation: ", sd.dip)
```

    ## [1] "Haploids: 0.0106648055111777 Deviation: 0.0153652158371832"

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

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
    ## REML criterion at convergence: -2710.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6234 -0.3793  0.0585  0.5134  5.3879 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  day      (Intercept) 0.0004020 0.02005 
    ##  Residual             0.0005332 0.02309 
    ## Number of obs: 587, groups:  day, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)    0.160375   0.007082   22.64
    ## ploidyHaploid -0.006090   0.003604   -1.69
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## ploidyHapld -0.276

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
    ## null    3 -2719.3 -2706.2 1362.7  -2725.3                       
    ## full    4 -2720.1 -2702.6 1364.0  -2728.1 2.7564  1    0.09687 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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
    ## REML criterion at convergence: -4703.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8607 -0.4622  0.0725  0.4270 13.6502 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001356 0.01164 
    ##  day      (Intercept) 0.0003378 0.01838 
    ##  Residual             0.0009650 0.03106 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)   0.132935   0.009160  14.513
    ## MAMA          0.006239   0.005617   1.111
    ## ploidyHaploid 0.026321   0.006656   3.955
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA  
    ## MAMA        -0.609       
    ## ploidyHapld -0.487  0.178

``` r
full <- lmer(slope ~ MA*ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ MA * ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -4704.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8683 -0.4815  0.0680  0.4190 13.7110 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001019 0.01009 
    ##  day      (Intercept) 0.0003256 0.01805 
    ##  Residual             0.0009704 0.03115 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error t value
    ## (Intercept)         0.159062   0.012186  13.053
    ## MAMA               -0.022494   0.010857  -2.072
    ## ploidyHaploid      -0.004656   0.012027  -0.387
    ## MAMA:ploidyHaploid  0.035998   0.012218   2.946
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA   pldyHp
    ## MAMA        -0.811              
    ## ploidyHapld -0.781  0.791       
    ## MAMA:pldyHp  0.702 -0.886 -0.842

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
    ## null    6 -4716.8 -4686.3 2364.4  -4728.8                        
    ## full    7 -4722.9 -4687.4 2368.4  -4736.9 8.1513  1   0.004303 **
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
    ## REML criterion at convergence: -4703.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8607 -0.4622  0.0725  0.4270 13.6502 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001356 0.01164 
    ##  day      (Intercept) 0.0003378 0.01838 
    ##  Residual             0.0009650 0.03106 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)   0.132935   0.009160  14.513
    ## MAMA          0.006239   0.005617   1.111
    ## ploidyHaploid 0.026321   0.006656   3.955
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) MAMA  
    ## MAMA        -0.609       
    ## ploidyHapld -0.487  0.178

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
    ## null    5 -4717.5 -4692.1 2363.7  -4727.5                     
    ## full    6 -4716.8 -4686.3 2364.4  -4728.8 1.2785  1     0.2582

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
    ## REML criterion at convergence: -4703.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8607 -0.4622  0.0725  0.4270 13.6502 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001356 0.01164 
    ##  day      (Intercept) 0.0003378 0.01838 
    ##  Residual             0.0009650 0.03106 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)   0.132935   0.009160  14.513
    ## ploidyHaploid 0.026321   0.006656   3.955
    ## MAMA          0.006239   0.005617   1.111
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) pldyHp
    ## ploidyHapld -0.487       
    ## MAMA        -0.609  0.178

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
    ## null    5 -4704.1 -4678.7 2357.1  -4714.1                         
    ## full    6 -4716.8 -4686.3 2364.4  -4728.8 14.659  1  0.0001288 ***
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
    ## -0.035092 -0.012297  0.000974  0.011070  0.051601 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -0.0233731  0.0056782  -4.116 8.20e-05 ***
    ## ploidy.xHaploid            0.0356016  0.0077492   4.594 1.33e-05 ***
    ## mutations                  0.0001471  0.0002647   0.556    0.580    
    ## ploidy.xHaploid:mutations -0.0004861  0.0010440  -0.466    0.643    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01731 on 95 degrees of freedom
    ## Multiple R-squared:  0.4597, Adjusted R-squared:  0.4426 
    ## F-statistic: 26.94 on 3 and 95 DF,  p-value: 1.077e-12

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
    ## -0.035228 -0.012499  0.001083  0.011846  0.050515 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     -0.0227684  0.0055051  -4.136 7.58e-05 ***
    ## ploidy.xHaploid  0.0328990  0.0051128   6.435 4.85e-09 ***
    ## mutations        0.0001158  0.0002550   0.454    0.651    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01723 on 96 degrees of freedom
    ## Multiple R-squared:  0.4584, Adjusted R-squared:  0.4472 
    ## F-statistic: 40.63 on 2 and 96 DF,  p-value: 1.642e-13

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
    ## -0.047965 -0.011927 -0.001316  0.011294  0.061024 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.0080667  0.0032249   2.501    0.014 *  
    ## mutations   -0.0010907  0.0002056  -5.305 7.12e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02051 on 97 degrees of freedom
    ## Multiple R-squared:  0.2249, Adjusted R-squared:  0.2169 
    ## F-statistic: 28.14 on 1 and 97 DF,  p-value: 7.118e-07

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
    ## REML criterion at convergence: -4696.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8543 -0.4552  0.0745  0.4178 13.6586 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001365 0.01168 
    ##  day      (Intercept) 0.0003484 0.01867 
    ##  Residual             0.0009655 0.03107 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)   1.388e-01  8.594e-03  16.153
    ## mutations     1.796e-05  2.460e-04   0.073
    ## ploidyHaploid 2.522e-02  7.426e-03   3.396
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns
    ## mutations   -0.518       
    ## ploidyHapld -0.605  0.466

``` r
full <- lmer(slope ~ mutations*ploidy + (1|day) + (1|lineid), df)
summary(full)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: slope ~ mutations * ploidy + (1 | day) + (1 | lineid)
    ##    Data: df
    ## 
    ## REML criterion at convergence: -4684.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8580 -0.4619  0.0676  0.4204 13.6255 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001319 0.01149 
    ##  day      (Intercept) 0.0003421 0.01850 
    ##  Residual             0.0009669 0.03110 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                           Estimate Std. Error t value
    ## (Intercept)              1.401e-01  8.601e-03  16.292
    ## mutations               -6.823e-05  2.563e-04  -0.266
    ## ploidyHaploid            2.102e-02  8.300e-03   2.533
    ## mutations:ploidyHaploid  8.965e-04  8.342e-04   1.075
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns pldyHp
    ## mutations   -0.529              
    ## ploidyHapld -0.595  0.533       
    ## mttns:pldyH  0.134 -0.306 -0.456

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
    ## null    6 -4715.5 -4685.1 2363.7  -4727.5                     
    ## full    7 -4714.7 -4679.2 2364.3  -4728.7 1.2067  1      0.272

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
    ## REML criterion at convergence: -4696.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8543 -0.4552  0.0745  0.4178 13.6586 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001365 0.01168 
    ##  day      (Intercept) 0.0003484 0.01867 
    ##  Residual             0.0009655 0.03107 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)   1.388e-01  8.594e-03  16.153
    ## mutations     1.796e-05  2.460e-04   0.073
    ## ploidyHaploid 2.522e-02  7.426e-03   3.396
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) muttns
    ## mutations   -0.518       
    ## ploidyHapld -0.605  0.466

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
    ##      npar     AIC     BIC logLik deviance Chisq Df Pr(>Chisq)
    ## null    5 -4717.5 -4692.1 2363.7  -4727.5                    
    ## full    6 -4715.5 -4685.1 2363.7  -4727.5 0.005  1     0.9434

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
    ## REML criterion at convergence: -4696.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8543 -0.4552  0.0745  0.4178 13.6586 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  lineid   (Intercept) 0.0001365 0.01168 
    ##  day      (Intercept) 0.0003484 0.01867 
    ##  Residual             0.0009655 0.03107 
    ## Number of obs: 1177, groups:  lineid, 105; day, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)   1.388e-01  8.594e-03  16.153
    ## ploidyHaploid 2.522e-02  7.426e-03   3.396
    ## mutations     1.796e-05  2.460e-04   0.073
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) pldyHp
    ## ploidyHapld -0.605       
    ## mutations   -0.518  0.466

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
    ## null    5 -4706.3 -4681.0 2358.2  -4716.3                         
    ## full    6 -4715.5 -4685.1 2363.7  -4727.5 11.162  1   0.000835 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
