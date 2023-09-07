# this fits a spline to y=OD versus x=time, and then finds the maximum slope
# code by 
spline.slope<-function(x, y, n=101, eps=1e-5, span=0.2){
  max(nderiv(loess(log(y) ~ x, degree=1, span=span), seq(min(x), max(x), length=n)), na.rm=TRUE)
}

# used by the function above to get a local (linear) slope around a point
nderiv <- function(fit, x, eps=1e-5){
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}


# load the file that indicates treatments in different wells
d96 <- read.csv("Trial96 2 ref.csv")
d384 <- read.csv("Trial384 2 ref.csv")

# load the OD data
assay96.data <- read.delim("Trial96 2 - Copy (2).txt")
assay384.data <- read.delim("Trial384 3 - Copy.txt")

# changing the format of time to be quarter hour increments
assay96.data$Time<-seq(from=0.25,by=0.25,length.out=nrow(assay96.data))
assay384.data$Time<-seq(from=0.25,by=0.25,length.out=nrow(assay384.data))

# calculate slopes
d96$slope <- sapply(c(1:nrow(d96)),function(r){
  spline.slope(assay96.data$Time, assay96.data[,which(names(assay96.data)==d96$Well[r])])
})
d384$slope <- sapply(c(1:nrow(d384)),function(r){
  spline.slope(assay384.data$Time, assay384.data[,which(names(assay384.data)==d384$Well[r])])
})

