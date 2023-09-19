# Author: Dr Nathaniel Sharp

# this fits a spline to y=OD versus x=time, and then finds the maximum slope
# code by 
spline.slope<-function(x, y, n=101, eps=1e-5, span=0.2){
  max(nderiv(loess(log(y) ~ x, degree=1, span=span), seq(min(x), max(x), length=n)), na.rm=TRUE)
}

# used by the function above to get a local (linear) slope around a point
nderiv <- function(fit, x, eps=1e-5){
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}


# load the file that indicates treatments in different wells
d <- read.csv("Ref 04 06.csv")

# load the OD data
assay.data <- read.delim("Rutuja 04 10.txt")

# changing the format of time to be quarter hour increments
assay.data$Time<-seq(from=0.25,by=0.25,length.out=nrow(assay.data))

# calculate slopes
d$slope <- sapply(c(1:nrow(d)), function(r){
  spline.slope(assay.data$Time, assay.data[,which(names(assay.data)==d$well[r])])
})

#Author: Rutuja Gupte

library(dplyr)

df <- d %>%
  group_by(treatment) %>% 
  summarise(avg = mean(slope))

ancestors <- c('H1', 'H2', 'H3', 'D1', 'D2', 'D3')
ancestors.haploid <- c('H1', 'H2', 'H3')
ancestors.diploid <- c('D1', 'D2', 'D3')
avg.haploid <- as.numeric(df %>% filter(treatment %in% ancestors.haploid) %>% summarize(avg = mean(avg)))
avg.diploid <- as.numeric(df %>% filter(treatment %in% ancestors.diploid) %>% summarize(avg = mean(avg)))


trial <- df %>% mutate(fit = ifelse(treatment == 'Blank', 0, ifelse(treatment %in% ancestors, 1,
                ifelse(as.numeric(treatment) %% 2 == 0, avg/avg.diploid,
                avg/avg.haploid)))) %>% mutate(category = ifelse(treatment == 'Blank', 0, 
                  ifelse(treatment %in% ancestors.haploid, 'Haploid',
                  ifelse(treatment %in% ancestors.haploid, 'Diploid',
                  ifelse(as.numeric(treatment) %% 2 == 0, 'H', 'D')))))

ggplot(trial %>% filter(category %in% c('H', 'D'))) +
  geom_histogram(aes(x=fit), color='black') +
  facet_grid(rows = vars(category)) + 
  geom_vline(xintercept=1, color='red')

trial %>% filter(category %in% c('H', 'D')) %>% summarise(mean(fit))
