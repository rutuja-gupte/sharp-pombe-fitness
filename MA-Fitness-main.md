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
d.test <- read.csv("data/Ref 04 08.csv")
assay.testdata <- read.delim("data/Rutuja 04 08.txt")

# d.test <- read.csv("test/Ref 10 11.csv")
# assay.testdata <- read.delim("test/Rutuja 10 11.txt")

d.test$initial <- sapply(c(1:nrow(d.test)), function(r){
  mean(assay.testdata[,which(names(assay.testdata)==d.test$well[r])][1:4])
})

d.test$final <- sapply(c(1:nrow(d.test)), function(r){
  temp <- assay.testdata[,which(names(assay.testdata)==d.test$well[r])]
  temp <- temp[!is.na(temp)]
  n <- length(temp)
  # print(assay.testdata[,which(names(assay.testdata)==d.test$well[r])][(n-3):n])
  return(mean(assay.testdata[,which(names(assay.testdata)==d.test$well[r])][(n-3):n]))
})

blanks <- d.test[d.test$treatment == "Blank",]
good.blanks <- blanks %>% filter(final - initial < 0.05)
blank.wells <- good.blanks$well
blank.val <- mean(unlist(assay.testdata[, blank.wells]), na.rm=TRUE)
# testcols <- colnames(assay.testdata)
assay.testdata2 <- bind_cols(assay.testdata[,1:2], data.frame(sapply(c(1:nrow(d.test)), function(r){
  temp <- (assay.testdata[,which(names(assay.testdata)==d.test$well[r])]-blank.val)
  temp <- replace(temp, which(temp<0), NA)
  if (d.test$well[r] ==  "I18"){
    print(which(names(assay.testdata)==d.test$well[r]))
    plot(assay.testdata[,"I18"])
  }
  return(temp)
  })))
```

    ## [1] 212

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
colnames(assay.testdata2)
```

    ##   [1] "Time"  "T.500" "X1"    "X2"    "X3"    "X4"    "X5"    "X6"    "X7"   
    ##  [10] "X8"    "X9"    "X10"   "X11"   "X12"   "X13"   "X14"   "X15"   "X16"  
    ##  [19] "X17"   "X18"   "X19"   "X20"   "X21"   "X22"   "X23"   "X24"   "X25"  
    ##  [28] "X26"   "X27"   "X28"   "X29"   "X30"   "X31"   "X32"   "X33"   "X34"  
    ##  [37] "X35"   "X36"   "X37"   "X38"   "X39"   "X40"   "X41"   "X42"   "X43"  
    ##  [46] "X44"   "X45"   "X46"   "X47"   "X48"   "X49"   "X50"   "X51"   "X52"  
    ##  [55] "X53"   "X54"   "X55"   "X56"   "X57"   "X58"   "X59"   "X60"   "X61"  
    ##  [64] "X62"   "X63"   "X64"   "X65"   "X66"   "X67"   "X68"   "X69"   "X70"  
    ##  [73] "X71"   "X72"   "X73"   "X74"   "X75"   "X76"   "X77"   "X78"   "X79"  
    ##  [82] "X80"   "X81"   "X82"   "X83"   "X84"   "X85"   "X86"   "X87"   "X88"  
    ##  [91] "X89"   "X90"   "X91"   "X92"   "X93"   "X94"   "X95"   "X96"   "X97"  
    ## [100] "X98"   "X99"   "X100"  "X101"  "X102"  "X103"  "X104"  "X105"  "X106" 
    ## [109] "X107"  "X108"  "X109"  "X110"  "X111"  "X112"  "X113"  "X114"  "X115" 
    ## [118] "X116"  "X117"  "X118"  "X119"  "X120"  "X121"  "X122"  "X123"  "X124" 
    ## [127] "X125"  "X126"  "X127"  "X128"  "X129"  "X130"  "X131"  "X132"  "X133" 
    ## [136] "X134"  "X135"  "X136"  "X137"  "X138"  "X139"  "X140"  "X141"  "X142" 
    ## [145] "X143"  "X144"  "X145"  "X146"  "X147"  "X148"  "X149"  "X150"  "X151" 
    ## [154] "X152"  "X153"  "X154"  "X155"  "X156"  "X157"  "X158"  "X159"  "X160" 
    ## [163] "X161"  "X162"  "X163"  "X164"  "X165"  "X166"  "X167"  "X168"  "X169" 
    ## [172] "X170"  "X171"  "X172"  "X173"  "X174"  "X175"  "X176"  "X177"  "X178" 
    ## [181] "X179"  "X180"  "X181"  "X182"  "X183"  "X184"  "X185"  "X186"  "X187" 
    ## [190] "X188"  "X189"  "X190"  "X191"  "X192"  "X193"  "X194"  "X195"  "X196" 
    ## [199] "X197"  "X198"  "X199"  "X200"  "X201"  "X202"  "X203"  "X204"  "X205" 
    ## [208] "X206"  "X207"  "X208"  "X209"  "X210"  "X211"  "X212"  "X213"  "X214" 
    ## [217] "X215"  "X216"  "X217"  "X218"  "X219"  "X220"  "X221"  "X222"  "X223" 
    ## [226] "X224"  "X225"  "X226"  "X227"  "X228"  "X229"  "X230"  "X231"  "X232" 
    ## [235] "X233"  "X234"  "X235"  "X236"  "X237"  "X238"  "X239"  "X240"  "X241" 
    ## [244] "X242"  "X243"  "X244"  "X245"  "X246"  "X247"  "X248"  "X249"  "X250" 
    ## [253] "X251"  "X252"  "X253"  "X254"  "X255"  "X256"  "X257"  "X258"  "X259" 
    ## [262] "X260"  "X261"  "X262"  "X263"  "X264"  "X265"  "X266"  "X267"  "X268" 
    ## [271] "X269"  "X270"  "X271"  "X272"  "X273"  "X274"  "X275"  "X276"  "X277" 
    ## [280] "X278"  "X279"  "X280"  "X281"  "X282"  "X283"  "X284"  "X285"  "X286" 
    ## [289] "X287"  "X288"  "X289"  "X290"  "X291"  "X292"  "X293"  "X294"  "X295" 
    ## [298] "X296"  "X297"  "X298"  "X299"  "X300"  "X301"  "X302"  "X303"  "X304" 
    ## [307] "X305"  "X306"  "X307"  "X308"  "X309"  "X310"  "X311"  "X312"  "X313" 
    ## [316] "X314"  "X315"  "X316"  "X317"  "X318"  "X319"  "X320"  "X321"  "X322" 
    ## [325] "X323"  "X324"  "X325"  "X326"  "X327"  "X328"  "X329"  "X330"  "X331" 
    ## [334] "X332"  "X333"  "X334"  "X335"  "X336"  "X337"  "X338"  "X339"  "X340" 
    ## [343] "X341"  "X342"  "X343"  "X344"  "X345"  "X346"  "X347"  "X348"  "X349" 
    ## [352] "X350"  "X351"  "X352"  "X353"  "X354"  "X355"  "X356"  "X357"  "X358" 
    ## [361] "X359"  "X360"  "X361"  "X362"  "X363"  "X364"  "X365"  "X366"  "X367" 
    ## [370] "X368"  "X369"  "X370"  "X371"  "X372"  "X373"  "X374"  "X375"  "X376" 
    ## [379] "X377"  "X378"  "X379"  "X380"  "X381"  "X382"  "X383"  "X384"

``` r
colnames(assay.testdata2) = c(names(assay.testdata[1:2]), d.test$well)

well <- assay.testdata2$I2
well.og <-  assay.testdata$I2
time <- seq(1, length(well))/4

ggplot() +
  geom_point(aes(x=time, y=well), color="green") +
  geom_point(aes(x=time, y=well.og), color="red")
```

    ## Warning: Removed 99 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 99 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
## Now trying to smooth

y = replace(well, which(well<=0), NA)
# print(assay.data$Time)
# print(y)
if (sum(is.na(y)) < (length(y)-50)) smoothed = predict(loess(log(y) ~ time, degree=1, span=0.075), time)
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : k-d tree limited by memory. ncmax= 200

``` r
# smoothed <- predict(loess(log(well) ~ time, degree=1, span=0.075), time)

derivs <- nderiv(loess(log(y) ~ time, degree=1, span=0.075), time)
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : k-d tree limited by memory. ncmax= 200

``` r
ggplot() +
  geom_point(aes(x=time, y=log(y)))
```

    ## Warning: Removed 99 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggplot() +
  geom_point(aes(x=time, y=derivs))
```

    ## Warning: Removed 101 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

``` r
ggplot() +
  geom_point(aes(x=time, y=smoothed), color="red") +
  geom_point(aes(x=time, y=log(well)))
```

    ## Warning: Removed 99 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 99 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

``` r
max(derivs, na.rm=TRUE)
```

    ## [1] 0.5335337

``` r
derivs
```

    ##   [1]            NA  7.262536e-02  9.223580e-02  1.163979e-01  2.257661e-06
    ##   [6] -2.257661e-06 -1.164002e-01 -2.655263e-01 -1.491283e-01  1.491261e-01
    ##  [11]  1.491283e-01 -3.272810e-02  2.257661e-06  1.491283e-01  1.163980e-01
    ##  [16]  4.152501e-06  9.769603e-02  2.228628e-01  3.070401e-01  2.897239e-01
    ##  [21]  2.659813e-01  3.188200e-01  4.328295e-01  5.335337e-01  4.848870e-01
    ##  [26]  3.750090e-01  3.298419e-01  3.891355e-01  4.080447e-01  3.653383e-01
    ##  [31]  3.538425e-01  3.927831e-01  4.112597e-01  3.931340e-01  3.821451e-01
    ##  [36]  3.693363e-01  3.379043e-01  3.249010e-01  3.310898e-01  3.363353e-01
    ##  [41]  3.443309e-01  3.452484e-01  3.391983e-01  3.218630e-01  3.002877e-01
    ##  [46]  2.923927e-01  2.766818e-01  2.706180e-01  2.789068e-01  2.973028e-01
    ##  [51]  2.760176e-01  2.239119e-01  2.043478e-01  2.137715e-01  2.289644e-01
    ##  [56]  2.399678e-01  2.397737e-01  2.189811e-01  1.892217e-01  1.725525e-01
    ##  [61]  1.705237e-01  1.819323e-01  1.945618e-01  1.941166e-01  2.014116e-01
    ##  [66]  2.013711e-01  1.920351e-01  1.824068e-01  1.814408e-01  1.793118e-01
    ##  [71]  1.659906e-01  1.525540e-01  1.595679e-01  1.623049e-01  1.424503e-01
    ##  [76]  1.088657e-01  1.000404e-01  1.003164e-01  8.801192e-02  8.505201e-02
    ##  [81]  8.430373e-02  7.933532e-02  6.928644e-02  7.028984e-02  8.199254e-02
    ##  [86]  9.253363e-02  9.867060e-02  1.015723e-01  9.307107e-02  8.245061e-02
    ##  [91]  7.971039e-02  8.084335e-02  7.405038e-02  6.884475e-02  6.747276e-02
    ##  [96]  6.058062e-02  4.800336e-02  4.097239e-02  4.409588e-02  5.086926e-02
    ## [101]  5.019120e-02            NA            NA            NA            NA
    ## [106]            NA            NA            NA            NA            NA
    ## [111]            NA            NA            NA            NA            NA
    ## [116]            NA            NA            NA            NA            NA
    ## [121]            NA            NA            NA            NA            NA
    ## [126]            NA            NA            NA            NA            NA
    ## [131]            NA            NA            NA            NA            NA
    ## [136]            NA            NA            NA            NA            NA
    ## [141]            NA            NA            NA            NA            NA
    ## [146]            NA            NA            NA            NA            NA
    ## [151]            NA            NA            NA            NA            NA
    ## [156]            NA            NA            NA            NA            NA
    ## [161]            NA            NA            NA            NA            NA
    ## [166]            NA            NA            NA            NA            NA
    ## [171]            NA            NA            NA            NA            NA
    ## [176]            NA            NA            NA            NA            NA
    ## [181]            NA            NA            NA            NA            NA
    ## [186]            NA            NA            NA            NA            NA
    ## [191]            NA            NA            NA            NA            NA
    ## [196]            NA            NA            NA            NA            NA
    ## [201]            NA

``` r
spline.slope(time, well)
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : k-d tree limited by memory. ncmax= 200

    ## [1] 0.5335337

``` r
# 0.7*(1-1/1000)
# # slopes <- nderiv(log(well), time)
# # fitted.slopes <- nderiv(log(well), time)
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



# d <- d %>% mutate(category = ifelse(label %in% fake.diploids, 'MA.H', category))
# d <- d %>% mutate(category = ifelse(label %in% fake.haploids, 'MA.D', category))
# 
# d <- d %>% filter(!(label %in% fake.diploids))
# d <- d %>% filter(!(label %in% fake.haploids))

# Labeling the dates
dates <- d %>% distinct(date)
# alphabet <- c('A', 'B')
alphabet <- c('1A', '1B', '1C', '1D', '2A', '2B', '2C', '2D', '3A')
dates <- dates %>% mutate(day = alphabet)
d <- d %>% left_join(dates, by='date') %>% select(-date)

head(d)
```

    ##   well treatment      initial     final    slope double_time monotone
    ## 1   A1     Blank           NA        NA       NA          NA      100
    ## 2   B1        H1 0.0020535098 0.7138035 2.640040   0.2625517        0
    ## 3   C1        H2 0.0020535098 0.7620535 2.797476   0.2477759        0
    ## 4   D1        H3 0.0003035098 0.6668035 2.941911   0.2356112        0
    ## 5   E1        D1 0.0015535098 0.5585535 2.865809   0.2418679        0
    ## 6   F1        D2 0.0023035098 0.4658035 3.025559   0.2290972        0
    ##   final_slope batch label category day
    ## 1          NA     1     0    Blank  1A
    ## 2 0.005953000     1   101   Ctrl.H  1A
    ## 3 0.003223910     1   103   Ctrl.H  1A
    ## 4 0.005398506     1   105   Ctrl.H  1A
    ## 5 0.003648250     1   102   Ctrl.D  1A
    ## 6 0.015899532     1   104   Ctrl.H  1A

``` r
d %>% head()
```

    ##   well treatment      initial     final    slope double_time monotone
    ## 1   A1     Blank           NA        NA       NA          NA      100
    ## 2   B1        H1 0.0020535098 0.7138035 2.640040   0.2625517        0
    ## 3   C1        H2 0.0020535098 0.7620535 2.797476   0.2477759        0
    ## 4   D1        H3 0.0003035098 0.6668035 2.941911   0.2356112        0
    ## 5   E1        D1 0.0015535098 0.5585535 2.865809   0.2418679        0
    ## 6   F1        D2 0.0023035098 0.4658035 3.025559   0.2290972        0
    ##   final_slope batch label category day
    ## 1          NA     1     0    Blank  1A
    ## 2 0.005953000     1   101   Ctrl.H  1A
    ## 3 0.003223910     1   103   Ctrl.H  1A
    ## 4 0.005398506     1   105   Ctrl.H  1A
    ## 5 0.003648250     1   102   Ctrl.D  1A
    ## 6 0.015899532     1   104   Ctrl.H  1A

``` r
d %>%
  filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>%
  filter(double_time < 5) %>%
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Who are these outliers????

``` r
d %>% filter(double_time > 100) %>% select(well, day, treatment)
```

    ##    well day treatment
    ## 1   K13  1B        97
    ## 2   L14  1B        10
    ## 3   G15  1B        17
    ## 4   J15  1B        20
    ## 5   I18  1B        55
    ## 6    C1  2C        H1
    ## 7    I1  2C        H1
    ## 8    L2  2C        H1
    ## 9    P6  2C        H1
    ## 10  G13  2C        H1
    ## 11  O13  2C        H1
    ## 12  O22  2C        H1
    ## 13   D3  2D        H2
    ## 14   A6  2D        H2

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

ci <- c(0.05, 2)

# ci <- c(quantile(d$slope, 0.025, na.rm=TRUE), quantile(d$slope, 0.975, na.rm=TRUE))

d %>% ggplot() +
  geom_density(aes(x=slope)) +
  geom_vline(xintercept = ci, color='red', linetype='dashed')
```

    ## Warning: Removed 153 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = ci, color="red")
```

    ## Warning: Removed 153 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
blanks <-  d %>% filter(category == 'Blank')
```

Plotting the distribution of the final OD values

``` r
d %>% ggplot() +
  geom_density(aes(x=final))
```

    ## Warning: Removed 70 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

visualize the blanks

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 60 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 45 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 45 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 45 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->

Setting cutoff for blanks around 0.05 and using the good blanks to make
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

Checking for effects of time (time taken to attain maximum growth rate).
Each time stamp is of 15 minutes.

Plotting the low values of time

Samples with time \< 10 were identified to be erroneous.

Looking to see if they reached saturation or not

``` r
ci_sat <- c(-0.015, 0.03)
ci_sat
```

    ## [1] -0.015  0.030

``` r
d %>% 
  filter(final_slope<0.25) %>%
  filter(final_slope >-0.25)  %>% 
  ggplot() +
  geom_histogram(aes(x=final_slope)) +
  geom_vline(xintercept = ci_sat, color='red', linetype='dashed')
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
# # sat <- c(0.2, 1.1)
# daily_cutoff <- d %>% group_by(day) %>%
#   summarize(lower = quantile(final, 0.05))
# 
# d2 <- d %>% left_join(daily_cutoff, by=c("day"))
# 
# d2 %>% ggplot() +
#   geom_density(aes(x=final)) +
#   geom_vline(aes(xintercept=lower), color='red', linetype='dashed') +
#   facet_wrap(vars(day))
# 
# daily_cutoff <- d %>% group_by(day) %>%
#   summarize(upper = quantile(initial, 0.95))
# 
# d3 <- d2 %>% left_join(daily_cutoff, by=c("day"))
# 
# d3 %>% ggplot() +
#   geom_density(aes(x=initial)) +
#   geom_vline(aes(xintercept=upper), color='red', linetype='dashed') +
#   facet_wrap(vars(day))
```

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
  anti_join(bad.blanks) %>%
  filter(initial <= final) %>%
  # filter(final >= lower) %>%
  # filter(initial <= upper)  %>%
  filter(final_slope < ci_sat[2]) %>%
  filter(final_slope > ci_sat[1]) %>%
  # filter(slope < ci[2]) %>%
  # filter(slope > ci[1])
  filter(monotone == 0) # %>%
```

    ## Joining with `by = join_by(well, treatment, initial, final, slope, double_time,
    ## monotone, final_slope, batch, label, category, day)`

``` r
  # filter(!(batch == 1 & treatment %in% c("D1", "D2", "D3"))) %>%
  # filter(!(batch == 2 & treatment %in% c("H1", "H2", "H3", "D2", "D3"))) %>%
  # filter(!(batch == 1 & category == "MA.D")) %>%
  # filter(!(batch == 2 & category == "MA.H"))

data %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](MA-Fitness-main_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# data %>% filter(slope > 0.3) %>% select(well, day)

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
