## ----setup, include=FALSE-------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----include=FALSE--------------------------------------------------------------------
library(dplyr)
library(stringr)
library(tidyverse)
library(lme4)
library(vioplot)


## -------------------------------------------------------------------------------------
spline.slope<-function(x, y, n=101, eps=1e-5, span=0.075){
  max(nderiv(loess(log(y) ~ x, degree=1, span=span), x), na.rm=TRUE)
}


## -------------------------------------------------------------------------------------
nderiv <- function(fit, x, eps=1e-5){
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}


## -------------------------------------------------------------------------------------
spline.time<-function(x, y, n=101, eps=1e-5, span=0.075){
  estimates <- loess(log(y) ~ x, degree=1, span=span)
  slopes <- nderiv(estimates, x)
  return(which.max(slopes))
}


## -------------------------------------------------------------------------------------
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
colnames(assay.testdata2)
colnames(assay.testdata2) = c(names(assay.testdata[1:2]), d.test$well)

well <- assay.testdata2$I2
well.og <-  assay.testdata$I2
time <- seq(1, length(well))/4

ggplot() +
  geom_point(aes(x=time, y=well), color="green") +
  geom_point(aes(x=time, y=well.og), color="red")

## Now trying to smooth

y = replace(well, which(well<=0), NA)
# print(assay.data$Time)
# print(y)
if (sum(is.na(y)) < (length(y)-50)) smoothed = predict(loess(log(y) ~ time, degree=1, span=0.075), time)
    
# smoothed <- predict(loess(log(well) ~ time, degree=1, span=0.075), time)

derivs <- nderiv(loess(log(y) ~ time, degree=1, span=0.075), time)


ggplot() +
  geom_point(aes(x=time, y=log(y)))

ggplot() +
  geom_point(aes(x=time, y=derivs))

ggplot() +
  geom_point(aes(x=time, y=smoothed), color="red") +
  geom_point(aes(x=time, y=log(well)))

max(derivs, na.rm=TRUE)
derivs
spline.slope(time, well)

# 0.7*(1-1/1000)
# # slopes <- nderiv(log(well), time)
# # fitted.slopes <- nderiv(log(well), time)


## ----include=FALSE--------------------------------------------------------------------
folder <-"data/"
l <- list.files(folder)
span = 0.075

dfs <- lapply(l[endsWith(l, '.csv')], function(r){
  d <- read.csv(paste0(folder,r))
  
  assay.name <- str_replace(str_replace(r, 'Ref', 'Rutuja'), '.csv', '.txt')
  assay.data <- read.delim(paste0(folder,assay.name))
  assay.data$Time<-seq(from=0.25,by=0.25,length.out=nrow(assay.data))
  
  ## Check and note the initial OD

  d$initial <- sapply(c(1:nrow(d)), function(r){
    mean(assay.data[,which(names(assay.data)==d$well[r])][1:4])
  })
  
  ## Check and note the final OD

  d$final <- sapply(c(1:nrow(d)), function(r){
    temp <- assay.data[,which(names(assay.data)==d$well[r])]
    temp <- temp[!is.na(temp)]
    n <- length(temp)
    mean(temp[(n-3):n])
  })

  ## Subtracting the blanks

  blanks <- d[d$treatment == "Blank",]
  good.blanks <- blanks %>% filter(final - initial < 0.05)
  blank.wells <- good.blanks$well
  blank.val <- mean(unlist(assay.data[, blank.wells]), na.rm=TRUE)
  testcols <- colnames(assay.data)
  assay.data <- bind_cols(assay.data[,1:2], data.frame(sapply(c(1:nrow(d)), function(r){
    temp <- assay.data[,which(names(assay.data)==d$well[r])]-blank.val
    temp <- replace(temp, which(temp<=0), NA)
    temp
    })))
  colnames(assay.data) = c(testcols[1:2], d$well)
  

  ## Now doing the initial and final calculations again

  ## Check and note the initial OD

  d$initial <- sapply(c(1:nrow(d)), function(r){
    temp <- assay.data[,which(names(assay.data)==d$well[r])]
    temp <- temp[!is.na(temp)]
    n <- length(temp)
    mean(temp[1:4])
  })

  ## Check and note the final OD

  d$final <- sapply(c(1:nrow(d)), function(r){
    temp <- assay.data[,which(names(assay.data)==d$well[r])]
    temp <- temp[!is.na(temp)]
    n <- length(temp)
    if (n<4) return(NA)
    mean(temp[(n-3):n])
  })

  ## Calculating the maximum slope mu_max
  ## We use this as a proxy for mu which is the inherent growth rate

  d$slope <- sapply(c(1:nrow(d)), function(r){
    y = assay.data[,which(names(assay.data)==d$well[r])]
    # y = replace(y, which(y<=0), NA)
    if (sum(is.na(y)) >= (length(y)-50)) return(NA)
    else spline.slope(assay.data$Time, y, span=span)
  })
  
  d$double_time <- log(2)/d$slope

  ## Dates in a good format are nice

  d$date <- mdy(str_replace(str_replace(r, 'Ref ', ''), '.csv', ' 2023'))

  d$monotone <- sapply(c(1:nrow(d)), function(r){
    y = assay.data[,which(names(assay.data)==d$well[r])]
    # y = replace(y, which(y<=0), NA)
    x = assay.data$Time
    if (sum(is.na(y)) >= (length(y)-50)) return(100)
    # y <- assay.data[9:nrow(assay.data),which(names(assay.data)==d$well[r])]
    # x <- assay.data$Time[9:nrow(assay.data)]
    return(sum(nderiv(loess(log(y) ~ x, degree=1, span=span), x) < 0, na.rm=TRUE))
  })

  d$final_slope <- sapply(c(1:nrow(d)), function(r){
    y <- assay.data[,which(names(assay.data)==d$well[r])]
    # y = replace(y, which(y<=0), NA)
    x <- assay.data$Time
    if (sum(is.na(y)) >= (length(y)-50)) return(NA)

    temp <- nderiv(loess(log(y) ~ x, degree=1, span=span), x)
    temp <- temp[!is.na(temp)]
    n = length(temp)
    if ((d$well[r]=="I18")&(d$date[r]==mdy("04/08/2023"))){
      print(blank.val)
      print(n)
      print(y)
    }
    return(mean(temp[(n-4):n]))
  })

  # # find the time stamps for the maximum growth rate
  # d$time <- sapply(c(1:nrow(d)), function(r){
  #   spline.time(assay.data$Time, assay.data[,which(names(assay.data)==d$well[r])], span=span)
  # })

  return(d)
})

d <- do.call("rbind", dfs)

head(d)


## -------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------
d %>% head()


## -------------------------------------------------------------------------------------
d %>%
  filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>%
  filter(double_time < 5) %>%
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## -------------------------------------------------------------------------------------
d %>% filter(double_time > 100) %>% select(well, day, treatment)


## ----echo=FALSE-----------------------------------------------------------------------
d %>%
  filter(double_time < 5) %>%
  filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=treatment, y=double_time, color=treatment)) + facet_wrap(vars(day))


## -------------------------------------------------------------------------------------
# ci <- c(quantile(d$slope, 0.25) - 1.5* IQR(d$slope),
# quantile(d$slope, 0.75) + 1.5* IQR(d$slope))

ci <- c(0.05, 2)

# ci <- c(quantile(d$slope, 0.025, na.rm=TRUE), quantile(d$slope, 0.975, na.rm=TRUE))

d %>% ggplot() +
  geom_density(aes(x=slope)) +
  geom_vline(xintercept = ci, color='red', linetype='dashed')

d %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = ci, color="red")


## -------------------------------------------------------------------------------------
blanks <-  d %>% filter(category == 'Blank')


## -------------------------------------------------------------------------------------
d %>% ggplot() +
  geom_density(aes(x=final))


## ----echo=FALSE-----------------------------------------------------------------------
blanks %>% ggplot() + geom_histogram(aes(x=slope)) +
  geom_vline(xintercept=0.06, color="red", linetype="dotted") +
  scale_x_log10()

blanks %>% ggplot() + geom_histogram(aes(x=initial)) +
  # geom_vline(xintercept=0.06, color="red", linetype="dotted") +
  scale_x_log10()

blanks %>% ggplot() + geom_histogram(aes(x=final)) # +
  # geom_vline(xintercept=0.06, color="red", linetype="dotted") +
  # scale_x_log10()

blanks %>% ggplot() + geom_histogram(aes(x=final-initial))


## -------------------------------------------------------------------------------------
# bad.blanks <- blanks %>% filter(slope > 0.06)
# good.blanks <- blanks %>% filter(slope < 0.06)
# 
# model <- lmer(slope~(1|day), data = good.blanks)
# dates.predict <- data.frame(date=distinct(d, day))
# dates.predict$null = predict(model, dates.predict)

bad.blanks <- blanks %>% filter(final - initial >= 0.05)
good.blanks <- blanks %>% filter(final - initial < 0.05)


## ----echo=FALSE-----------------------------------------------------------------------
# plotting the extremely low values of time
# d %>%
#   filter(time<50) %>%
#   ggplot() +
#   geom_histogram(aes(x=time))


## -------------------------------------------------------------------------------------
ci_sat <- c(-0.015, 0.03)
ci_sat
d %>% 
  filter(final_slope<0.25) %>%
  filter(final_slope >-0.25)  %>% 
  ggplot() +
  geom_histogram(aes(x=final_slope)) +
  geom_vline(xintercept = ci_sat, color='red', linetype='dashed')


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


## -------------------------------------------------------------------------------------
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
  # filter(!(batch == 1 & treatment %in% c("D1", "D2", "D3"))) %>%
  # filter(!(batch == 2 & treatment %in% c("H1", "H2", "H3", "D2", "D3"))) %>%
  # filter(!(batch == 1 & category == "MA.D")) %>%
  # filter(!(batch == 2 & category == "MA.H"))

data %>%
  # filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=category, y=slope, color=category)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# data %>% filter(slope > 0.3) %>% select(well, day)

data %>%
  filter(category == 'Ctrl.H' | category == 'Ctrl.D') %>% 
  ggplot() + geom_point(aes(x=treatment, y=slope, color=treatment)) + facet_grid(cols=vars(day)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cat(data$well[data$day=="1A"], file = "1A.txt")
cat(data$well[data$day=="1B"], file = "1B.txt")
cat(data$well[data$day=="1C"], file = "1C.txt")
cat(data$well[data$day=="1D"], file = "1D.txt")
cat(data$well[data$day=="2A"], file = "2A.txt")
cat(data$well[data$day=="2B"], file = "2B.txt")
cat(data$well[data$day=="2C"], file = "2C.txt")
cat(data$well[data$day=="2D"], file = "2D.txt")
cat(data$well[data$day=="3A"], file = "3A.txt")


## -------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------
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


## ----echo=FALSE-----------------------------------------------------------------------
# trt %>% ggplot() +
#   geom_histogram(aes(x=rel.fit), fill='black') +
#   facet_grid(rows = vars(ploidy)) +
#   geom_vline(data = data.frame(ploidy = c('Haploid', 'Diploid'),
#                                mu = c(mu.hap, mu.dip)),
#              aes(xintercept = mu, color = ploidy)) +
#   geom_vline(xintercept = 0, color = 'black') +
#   xlab("Relative Fitness") +
#   ylab("Number of Lines") +
#   guides(color = guide_legend(title = "Ploidy"))


## -------------------------------------------------------------------------------------
# null <- lmer(slope ~ 1 + (1|lineid) + (1|batch), ctrl)
# full <- lmer(slope ~ ploidy + (1|lineid) + (1|batch), ctrl)
# summary(full)
# 
# mod <- anova(null, full)
# mod
# 
# t.test(ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Haploid"],
#        ctrl$slope[ctrl$batch == 3 & ctrl$ploidy == "Diploid"])


## -------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------
# null <- lmer(slope ~ MA + ploidy + (1|day) + (1|lineid), df)
# summary(null)
# full <- lmer(slope ~ MA*ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod


## -------------------------------------------------------------------------------------
# # lmer for effect of MA
# null <- lmer(slope ~ 1 + ploidy + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ MA + ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod


## ----include=FALSE--------------------------------------------------------------------
# # Looking at the predictions made by the lmer model (understanding the model)
# 
# predictions <- predict(full, df[, c('MA', 'ploidy', 'day', 'lineid')])
# ggplot() +
#   geom_point(aes(x=df$slope, y=predictions)) +
#   geom_abline(slope = 1, intercept = 0, color='red') +
#   xlim(0, NA) +
#   ylim(0, NA)
# 
# ggplot() +
#   geom_point(aes(x=df$slope, y=fitted(full), color=df$day)) +
#   geom_abline(slope = 1, intercept = 0, color='red') +
#   xlim(0, NA) +
#   ylim(0, NA)
# 
# 
# ggplot() + geom_point(aes(x=fitted(full), y= resid(full),
#                           color=df$day)) #+ theme(legend.position = "none")
# 
# plot(full)
# summary(full)
# 
# ggplot() + geom_histogram(aes(x=resid(full)), fill='black')
# 
# qqnorm(resid(full))


## -------------------------------------------------------------------------------------
# null <- lmer(slope ~ 1 + MA + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ ploidy + MA + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod


## -------------------------------------------------------------------------------------
# # read the files
# mut <- read_delim("pombe_MA_data.txt")
# 
# # combine with the relative fitness data
# trt2 <- trt %>%
#   left_join(mut, by=c('label'='line'))
# trt2 <- trt2 %>% select(ploidy.x, label, rel.fit, ploidy.y, ploidy_final,
#                         n.SNM, n.indel)
# trt2 <- trt2 %>% mutate(mutations = n.SNM + n.indel)


## ----echo=FALSE-----------------------------------------------------------------------
# trt2 %>% ggplot() +
#   geom_point(aes(x=mutations , y=rel.fit, color=ploidy.x)) +
#   geom_hline(yintercept = 0) +
#   geom_point(data = trt2 %>% filter(label %in% aneuploids) %>%
#                filter(ploidy.x == 'Diploid'),
#              aes(x=mutations, y=rel.fit), color = 'red', shape='triangle', size = 3) +
#   geom_point(data = trt2 %>% filter(label %in% aneuploids) %>%
#                filter(ploidy.x == 'Haploid'),
#              aes(x=mutations, y=rel.fit), color = 'blue', shape='triangle', size = 3)


## -------------------------------------------------------------------------------------
# mod <- lm(rel.fit ~ ploidy.x*mutations + ploidy.x + mutations, trt2)
# summary(mod)
# qqnorm(resid(mod))


## -------------------------------------------------------------------------------------
# mod <- lm(rel.fit ~ ploidy.x + mutations, trt2)
# summary(mod)
# qqnorm(resid(mod))


## -------------------------------------------------------------------------------------
# mod <- lm(rel.fit ~ mutations, trt2)
# summary(mod)
# qqnorm(resid(mod))



## -------------------------------------------------------------------------------------
# mut_join <- trt2 %>% select(label, mutations)
# df <- df %>% left_join(mut_join) %>%
#   mutate(mutations = ifelse(is.na(mutations), 0, mutations))



## -------------------------------------------------------------------------------------
# null <- lmer(slope ~ mutations + ploidy + (1|day) + (1|lineid), df)
# summary(null)
# full <- lmer(slope ~ mutations*ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod


## -------------------------------------------------------------------------------------
# null <- lmer(slope ~ 1 + ploidy + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ mutations + ploidy + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod


## -------------------------------------------------------------------------------------
# null <- lmer(slope ~ 1 + mutations + (1|day) + (1|lineid), df)
# full <- lmer(slope ~ ploidy + mutations + (1|day) + (1|lineid), df)
# summary(full)
# 
# mod <- anova(null, full)
# mod

