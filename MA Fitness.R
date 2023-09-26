# Author: Dr Nathaniel Sharp

# this fits a spline to y=OD versus x=time, and then finds the maximum slope
# code by 
spline.slope<-function(x, y, n=101, eps=1e-5, span=0.2){
  max(nderiv(loess(log(y) ~ x, degree=1, span=span), seq(min(x), max(x), length=n)), na.rm=TRUE)
}

# used by the function above to get a local (linear) slope around a point
nderiv <- function(fit, x, eps=1e-5){
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}

#Author: Rutuja Gupte

library(dplyr)
library(stringr)
library(tidyverse)

l <- list.files()
dfs <- lapply(l[endsWith(l, '.csv')], function(r){
  d <- read.csv(r)
  
  assay.name <- str_replace(str_replace(r, 'Ref', 'Rutuja'), '.csv', '.txt')
  assay.data <- read.delim(assay.name)
  assay.data$Time<-seq(from=0.25,by=0.25,length.out=nrow(assay.data))
  # calculate slopes
  d$slope <- sapply(c(1:nrow(d)), function(r){
    spline.slope(assay.data$Time, assay.data[,which(names(assay.data)==d$well[r])])
  })
  
  d$fname <- r
  return(d)
})

d <- dfs %>% reduce(bind_rows)

df.raw <- d %>%
  group_by(treatment) %>% 
  summarise(avg = mean(slope))

dates.1 <- c('Ref 04 06.csv', 'Ref 04 08.csv', 'Ref 04 09.csv', 'Ref 04 10.csv')
dates.2 <- c('Ref 04 28.csv', 'Ref 04 30.csv', 'Ref 05 03.csv', 'Ref 05 05.csv')

ancestors <- c('H1', 'H2', 'H3', 'D1', 'D2', 'D3')
ancestors.haploid <- c('H1', 'H2', 'H3')
ancestors.diploid <- c('D1', 'D2', 'D3')
# avg.haploid <- as.numeric(df.raw %>% filter(treatment %in% ancestors.haploid) %>% summarize(avg = mean(avg)))
# avg.diploid <- as.numeric(df.raw %>% filter(treatment %in% ancestors.diploid) %>% summarize(avg = mean(avg)))

df <- df.raw
df$treatment[df$treatment == 'Blank'] <- '0'
df$treatment[df$treatment == 'H1'] <- '101'
df$treatment[df$treatment == 'D1'] <- '102'
df$treatment[df$treatment == 'H2'] <- '103'
df$treatment[df$treatment == 'D2'] <- '104'
df$treatment[df$treatment == 'H3'] <- '105'
df$treatment[df$treatment == 'D3'] <- '106'
df$treatment <- as.numeric(df$treatment)


df <- df %>%
  mutate(dates = ifelse(fname %in% dates.1, 1, 2)) %>%
  mutate(fit = case_when(treatment == 0 ~ 0,
                                    treatment > 100 ~ 1,
                                    treatment %% 2 == 0 ~ avg/avg.diploid,
                                    treatment %% 2 == 1 ~ avg/avg.haploid)) %>%
  mutate(category = case_when(treatment == 0 ~ 'Blank',
                              treatment > 100 & treatment %% 2 == 0 ~ 'Diploid',
                              treatment > 100 & treatment %% 2 == 1 ~ 'Haploid',
                              treatment %% 2 == 0 ~ 'D',
                              treatment %% 2 == 1 ~ 'H'))

df %>% filter(category %in% c('H', 'D')) %>% ggplot() +
  geom_histogram(aes(x=fit), color='black') +
  facet_grid(rows = vars(category)) + 
  geom_vline(xintercept=1, color='red')

#D1 from second set might be diploids. Filter for that.
dip.hunt <- d
dip.hunt <- d %>% filter(treatment %in% ancestors) %>%
  mutate(dip = ifelse(fname %in% 
                        dates.2 &
                        treatment == 'D1', 'Yes', 'No'))

dip.hunt %>% ggplot() +
  geom_point(aes(x=treatment, y=slope)) + 
  facet_grid(cols=vars(dip))

# Noticed that H3 from second set were doing funny stuff.
dip.hunt <- d
dip.hunt <- d %>% filter(treatment %in% ancestors) %>%
  mutate(dip = ifelse(fname %in% 
                        dates.2 &
                        treatment == 'H3', 'Yes', 'No'))

dip.hunt %>% ggplot() +
  geom_point(aes(x=treatment, y=slope)) + 
  facet_grid(cols=vars(dip))

#ANOVA for controls
ctrl <- d %>% mutate(dates = ifelse(fname %in% dates.1, 1, 2)) %>%
  filter(treatment %in% ancestors) %>%
  mutate(group = paste(treatment, dates, sep = ":"))

mod <- aov(ctrl$slope ~ ctrl$group)
summary(mod)
plot(fitted(mod), resid(mod))
qqnorm(resid(mod))

#tukey for all combinations
tukey.multiplier <- qtukey(0.95, nmeans = 12, df = 1418)/sqrt(2)
counts <- count(ctrl, group)
means <- ctrl %>% group_by(group) %>% summarize(average = mean(slope))
MS_e <- 0.00215

for (i in seq(1, 12)){
  for (j in seq(1,12)){
    moe_ij <- tukey.multiplier * sqrt(MS_e*(1/counts$n[i] + 1/counts$n[j]))
    pe_ij <- means$average[i] - means$average[j]
    ci_ij <- c(pe_ij - moe_ij, pe_ij + moe_ij)
    if (ci_ij[1] < 0 && ci_ij[2] > 0) {
      next
    }
    else {
      print(paste(counts$group[i], counts$group[j], sep=' - '))
    }
  }
}

#t test D1 vs everybody
d1 <- ctrl %>% filter(treatment == 'D1')
not.d1 <- ctrl %>% filter(treatment %in% ancestors & treatment != 'D1')
t.test(d1$slope, not.d1$slope, var.equal = F)

#t test D1:2 vs everybody
d1 <- ctrl %>% filter(group == 'D1:2')
not.d1 <- ctrl %>% filter(treatment %in% ancestors & group != 'D1:2')
t.test(d1$slope, not.d1$slope, var.equal = F)

#t test D2 vs everybody
d2 <- ctrl %>% filter(treatment == 'D2')
not.d2 <- ctrl %>% filter(treatment %in% ancestors & treatment != 'D2')
t.test(d2$slope, not.d2$slope, var.equal = F)

#t test D2 vs diploids
d2 <- ctrl %>% filter(treatment == 'D2')
not.d2 <- ctrl %>% filter(treatment %in% ancestors.diploid & treatment != 'D2')
t.test(d2$slope, not.d2$slope, var.equal = F)

#t test D3 vs everybody
d3 <- ctrl %>% filter(treatment == 'D3')
not.d3 <- ctrl %>% filter(treatment %in% ancestors & treatment != 'D3')
t.test(d3$slope, not.d3$slope, var.equal = T)
t.test(d3$slope, not.d3$slope, var.equal = F)


#t test haploids vs diploids
dip <- ctrl %>% filter(group == 'D1:2')
hap <- ctrl %>% filter(treatment %in% ancestors.haploid)
t.test(dip$slope, hap$slope, var.equal = F)

#t test D1:2 vs everybody
dip <- ctrl %>% filter(group == 'D1:2')
hap <- ctrl %>% filter(treatment %in% ancestors.haploid)
t.test(dip$slope, hap$slope, var.equal = F)



#t test mutant dip vs mutant hap

###########
# NEXT STEP: FOR REL FITNESS, USE CONTROLS FROM THE SAME ASSAY




trial %>% filter(category %in% c('H', 'D')) %>% ggplot() +
  geom_point(aes(x=as.numeric(treatment), y=fit, color=category)) +
  geom_smooth(aes(x=as.numeric(treatment), y=fit, color=category), method='lm', se=FALSE)

hap.mut <- trial %>% filter(category == 'H')
dip.mut <- trial %>% filter(category == 'D')
mut <- merge(hap.mut, dip.mut, by='row.names')

mut %>% ggplot() + geom_point(aes(x=fit.x, y=fit.y))
