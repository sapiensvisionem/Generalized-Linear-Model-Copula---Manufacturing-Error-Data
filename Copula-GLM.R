# Business analytics group of a manufacturing company is asked to investigate causes of malfunctions in technological process at one of the plants that result in significant increase of cost for the end product of the business.
# One of suspected reasons for malfunctions is deviation of temperature in burning zone from optimal levels. The sample in the provided file contains times of malfunctions in seconds since the start of measurement and records of temperature.

# Import Data 
# File contains time stamps of malfunction events expressed in seconds and temperature readings at the time of each event.
dataPath <- "/Users/jihun/Desktop"
# Temperature sensor updates readings once a minute.
Course.Project.Data<-read.csv(file=paste(dataPath,"projectdata.csv",sep="/"))
# Counting process is a step function that jumps by 1 at every time when event occurs.
Counting.Process <- as.data.frame(cbind(Time=Course.Project.Data$Time, Count=1:length(Course.Project.Data$Time)))
plot(counting_process$Time,counting_process$Count,type="s")
# The counting process trajectory looks pretty smooth and grows steadily. This may mean that malfunctions are caused by a persistent rather than episodic problem.

library(MASS)
library(pscl)
library(AER)

# Use the method based on comparison of the deviance and the number of degrees of freedom shown in the output of generalized linear model for Poisson regression. Calculate approximate confidence interval centered at the theoretical mean value of the deviance statistic and check if observed deviance belongs to it.
test.dat<-read.csv(file=paste(dataPath,"MScA_LinearNonlinear_step2_sample.csv",sep="/"))
Test.Deviance.Overdispersion.Poisson <- function(test.dat,Parameter.Lambda){
  my.Sample <- rpois(test.dat,Parameter.Lambda)
  Model <- glm(my.Sample~1,family=poisson)
  Dev <- Model$deviance
  Deg.Fred <- Model$df.residual
  (((Dev/Deg.Fred-1)/sqrt(2/Deg.Fred)>-1.96)&((Dev/Deg.Fred-1)/sqrt(2/Deg.Fred)<=1.96))*1
}
Test.Deviance.Overdispersion.Poisson(100,1)
sum(replicate(300,Test.Deviance.Overdispersion.Poisson(100,1)))
# 253
exp(glm(rpois(1000,2)~1,family=poisson)$coeff)
# 2.027
Test.Deviance.Overdispersion.NBinom <- function(test.dat,Parameter.prob){
  my.Sample <- rnbinom(test.dat,2,Parameter.prob)
  Model <- glm(my.Sample~1,family=poisson)
  Dev <- Model$deviance
  Deg.Fred <- Model$df.residual
  (((Dev/Deg.Fred-1)/sqrt(2/Deg.Fred)>-1.96)&((Dev/Deg.Fred-1)/sqrt(2/Deg.Fred)<=1.96))*1
}
sum(replicate(300,Test.Deviance.Overdispersion.NBinom(100,.2)))
# 0
GLM.model <- glm(test.dat$count~1,family=poisson)
summary(GLM.model)
# 
-1.96*sqrt(2*499) + 499
1.96*sqrt(2*499)+499

# Use regression-based tests for over-dispersion by Cameron-Trivedi implemented in dispersiontest() from library AER.

Apply dispersiontest() to the glm model fitted in the previous section.
Disp.Test <- dispersiontest(GLM.model,alternative="two.sided")
Disp.Test
# t-test = 0.041314
# p-value - 0.967
# underdispersion

# Use method based on negative binomial regression.

Apply glm.nb() from MASS to one-minute counts to fit a negative binomial model.
Then use odTest() from library pscl to test with level of 5% if the data can be described by Poisson distribution (no over-dispersion) or not (over-dispersion).
NB.model <- glm.nb(test.dat$count~1)
NB.model
odTest(NB.model)
# 0.0016 0.484 

###########################################
# Kolmogorov-Smirnov test is used to test hypotheses of equivalence between two empirical distributions or equivalence between one empirical distribution and one theoretical distribution.
# check if it is a Poisson flow of events with parameter equal to average intensity (number of events per minute).
Minute.times <- matrix(data=NA,nrow=250,ncol=1)
Minute.times[,1] <- row(Minute.times)[,1]*60-30
colnames(Minute.times) <- "Minute.times"
Course.Project.Data$Minute.times <- ceiling(Course.Project.Data$Time/60)*60-30
Course.Project.Data$One <- 1
Minute.counts <- aggregate(One~Temperature+Minute.times, data=Course.Project.Data, sum)
One.Minute.Counts.Temps <- merge(Minute.times,Minute.counts,by="Minute.times",all.x=TRUE)
One.Minute.Counts.Temps <- One.Minute.Counts.Temps[,c(1,3,2)]
colnames(One.Minute.Counts.Temps) <- c("Minute.times","Minute.counts","Minute.Temps")
One.Minute.Counts.Temps$Minute.counts[c(is.na(One.Minute.Counts.Temps$Minute.counts))] <- 0
test.dist<-function(data){

# enter Kolmogorov-Smirnov statistics and their p-values for all 5 intensity distribution candidates: normal, exponential, gamma, lognormal, logistic. Select distribution that is most consistent with the sample according to Kolmogorov-Smirnov test.  dist=c('normal', 'exponential', 'gamma', 'lognormal', 'logistic')
cdf=c('pnorm','pexp','pgamma', 'plnorm', 'plogis')
df<-data.frame(dist,cdf)
  
data.mean<-mean(data)
data.sd<-sd(data)
  
dist <- list()
cdf <- list()
stat<-list()
pval<-list()
Parameter1<-list()
Parameter2<-list()
  
for (i in 1:dim(df)[1]){
    
  icdf<-toString(df$cdf[i])
  idist<-toString(df$dist[i])
    
  d<-fitdistr(data, idist)#, data.mean, data.sd)
  if(length(d$estimate)==1) {
    d.ks<-ks.test(data, icdf, d$estimate)
    Parameter1[[i]] <- d$estimate
    Parameter2[[i]] <- NULL
      
  } else {
    d.ks<-ks.test(data, icdf, d$estimate[1], d$estimate[2])
    Parameter1[[i]] <- d$estimate[1]
    Parameter2[[i]] <- d$estimate[2]
  }
    
    
  dist[[i]]<-idist
  cdf[[i]]<-icdf
  stat[[i]]<-d.ks$statistic
  pval[[i]]<-d.ks$p.value
    
}
  
res<-data.frame(cbind(dist,cdf,stat,pval, Parameter1, Parameter2))
  
return(res)
}
(countParam <- test.dist(One.Minute.Counts.Temps$Minute.counts))
(tempParam <- test.dist(One.Minute.Counts.Temps$Minute.Temps))
write.csv(One.Minute.Counts.Temps,file=paste(dataPath,"OneMinuteCountsTemps.csv",sep="/"),row.names=FALSE)
write.table(distrParam,file=paste(dataPath,"DistrParameters.csv",sep="/"))

suppressWarnings((library(MASS)))
suppressWarnings(library(copula))
Part2.Data <- read.csv(file=paste(dataPath,"OneMinuteCountsTemps.csv",sep="/"))
dim(Part2.Data)
Part2.Data <- as.data.frame(cbind(Part2.Data,Part2.Data[,2]/60)) #per minute intensity
colnames(Part2.Data) <- c("Times","Counts","Temperatures","Intensities")

# Explore possible types of dependence between one-minute counts and temperatures.

# create a data frame with one-minute breaks counts and temperature measurements.

dat <- readRDS(paste(dataPath, 'nonlinear_models_course_assignment_step_5_data.rds', sep = '/'))
head(dat)
plot(dat$predictor, dat$output)
y <- dat$output
x <- dat$predictor
u <- rank(y)/length(y)
v <- rank(x)/length(x)
plot(v,u)

library(copula)


# First step : 
copulaFitData <- cbind(dat$predictor,dat$output)

copula.type <- c("normal")
Copula.Fit.Object <- normalCopula(param = 0, dim = 2)
Copula.fit.Gaussian  <- fitCopula(Copula.Fit.Object, 
                                  pobs(copulaFitData, ties.method = "average"), 
                                  method = "ml",
                                  optim.method = "BFGS", 
                                  optim.control = list(maxit = 1000))

Copula.fit.Gaussian@loglik
# 255.3742


copula.type <- c("Clayton")
Copula.Fit.Object <- claytonCopula(param = 5, dim = 2)
Copula.fit.Clayton <- fitCopula(Copula.Fit.Object, 
                                pobs(copulaFitData, ties.method = "average"), 
                                method = "ml",
                                optim.method = "BFGS", 
                                optim.control = list(maxit = 1000))
Copula.fit.Clayton@loglik # 92.16895


copula.type <- c("Frank")
Copula.Fit.Object <- frankCopula(param = 5, dim = 2)
Copula.fit.Frank <- fitCopula(Copula.Fit.Object, 
                              pobs(copulaFitData, ties.method = "average"), 
                              method = "ml",
                              optim.method = "BFGS", 
                              optim.control = list(maxit = 1000))
Copula.fit.Frank@loglik # 234.7638

copulaType <- "normal"

dat$predictor_DistrType; dat$predictor_DistrParameters; dat$output_DistrType; dat$output_DistrParameters

# Calculate copula image of the predictor
normal.mean <- 90
normal.variance <- 2
u_vector <- pobs(dat$predictor,ties.method="average")

# Estimate 5%, 95% and 50% image quantiles
theta <- Copula.fit.Gaussian@estimate
lowMidHigh <- as.data.frame(cbind(dat$predictor, dat$output))
colnames(lowMidHigh) <- c('predictor', 'output')

alpha<-.50
lowMidHigh$intensMid <- sapply(u_vector,
                               function(z)
                                 pnorm(qnorm(alpha)*sqrt(1-theta^2) + theta*qnorm(z)))
alpha<-.05
lowMidHigh$intensLow <- sapply(u_vector,
                               function(z)
                                 pnorm(qnorm(alpha)*sqrt(1-theta^2) + theta*qnorm(z)))
alpha<-.95
lowMidHigh$intensHigh <- sapply(u_vector,
                                function(z)
                                  pnorm(qnorm(alpha)*sqrt(1-theta^2) + theta*qnorm(z)))
quantileLow <- lowMidHigh$intensLow
quantileMid <- lowMidHigh$intensMid
quantileHigh <- lowMidHigh$intensHigh

# identify minutes during which intensities were anomalously high relative to the 95% level 
# of conditional distribution of intensities for a given temperature level

# sample, identify and fit one parametric copula of the three (normal, Frank, Clayton) and calculate 3 outputs of quantile regression for levels 5%,50% and 95%, for the given marginal distributions of both variables.

library(copula)
copula_data <- as.data.frame(copula_data)
Clayton.fit <- fitCopula(
  copula = claytonCopula(),
  data = 1-copula_data,
  method = "ml",
  optim.method = "BFGS")
head(copula_data)
# now find the anomalies
theta <- Clayton.fit@estimate
alpha<-0.95
imgY <- -pgamma(copula_data$Intensities,1.73511795,0.12705900)+1
imgX <- -pnorm(copula_data$Temperatures,mean=100.2008143,sd=5.2210835)+1

lowBound <- sapply(imgX, 
                   function(z)
                     ((0.05^(-theta/(1+theta))-1)*z^(-theta)+1)^(-1/theta))
midBound <- sapply(imgX, 
                   function(z)
                     ((0.5^(-theta/(1+theta))-1)*z^(-theta)+1)^(-1/theta))
highBound <- sapply(imgX, 
                    function(z)
                      ((alpha^(-theta/(1+theta))-1)*z^(-theta)+1)^(-1/theta))

copula <- -pobs(copula_data)+1
anomHighIdx <- (imgY<lowBound)
anomHighIdx
anomalies <- copula[anomHighIdx,]

intensity <- copula_data$intensity
temp <- copula_data$mean
anomalies_intensity <- intensity[anomHighIdx]
anomalies_temperature <- temp[anomHighIdx]
# Check the graph
plot(part2_data$mean,part2_data$intensity)
plot(anomalies_temperature,anomalies_intensity )

dataPath <- "/Users/jihun/Desktop"
part2_data<-read.csv(file=paste(dataPath,"OneMinuteCountsTemps.csv",sep="/"))
#head(part2_data)
library(copula)
part2_data<-part2_data[complete.cases(part2_data),]
part2_data$intensity<-part2_data[,"Minute.counts"]
copula_data <- copula::pobs(part2_data[,c("Minute.Temps","intensity")])
Clayton.fit <- fitCopula(
  copula = claytonCopula(),
  data = 1-copula_data,
  method = "ml",
  optim.method = "BFGS")
head(copula_data)
head(part2_data)
# now find the anomalies
theta <- Clayton.fit@estimate
alpha<-0.95
imgY <- -pgamma(part2_data$intensity,1.73511795,0.12705900)+1
imgX <- -pnorm(part2_data$Minute.Temps,mean=100.2008143,sd=5.2210835)+1

lowBound <- sapply(imgX, 
                   function(z)
                     ((0.05^(-theta/(1+theta))-1)*z^(-theta)+1)^(-1/theta))
midBound <- sapply(imgX, 
                   function(z)
                     ((0.5^(-theta/(1+theta))-1)*z^(-theta)+1)^(-1/theta))
highBound <- sapply(imgX, 
                    function(z)
                      ((alpha^(-theta/(1+theta))-1)*z^(-theta)+1)^(-1/theta))

copula2 <- -pobs(part2_data)+1
anomHighIdx <- (imgY<lowBound)
anomHighIdx
anomalies <- copula2[anomHighIdx,]

intensity <- part2_data$intensity
temp <- part2_data$Minute.Temps
anomalies_intensity <- intensity[anomHighIdx]
anomalies_temperature <- temp[anomHighIdx]
# Check the graph
plot(part2_data$Minute.Temps,part2_data$intensity)
plot(anomalies_temperature,anomalies_intensity )


dataPath <- "/Users/jihun/Desktop"
dat <- readRDS(paste(dataPath, 'nonlinear_models_course_assignment_step_5_data.rds', sep = '/'))
dat$predictor_DistrType; dat$predictor_DistrParameters; dat$output_DistrType; dat$output_DistrParameters
# predictor: gamma 4,4  
# output: normal 90, 2
#View(dat2)
library(copula)
Y <- dat$output
X <- dat$predictor

CopulaFitData <- cbind(X,Y)

# normal
Copula.Fit.Object <- normalCopula(param = 0, dim = 2)
Copula.fit.Gaussian  <- fitCopula(Copula.Fit.Object, 
                                  pobs(CopulaFitData, ties.method = "average"), 
                                  method = "ml",
                                  optim.method = "BFGS", 
                                  optim.control = list(maxit = 1000))

Copula.fit.Gaussian@loglik # 255.3742

# Frank copula
Copula.Fit.Object <- frankCopula(param = 5, dim = 2)
Copula.fit.Frank <- fitCopula(Copula.Fit.Object, 
                              pobs(CopulaFitData, ties.method = "average"), 
                              method = "ml",
                              optim.method = "BFGS", 
                              optim.control = list(maxit = 1000))
Copula.fit.Frank@loglik # 234.7638

# Fit Clayton copula
Copula.Fit.Object <- claytonCopula(param = 5, dim = 2)
Copula.fit.Clayton <- fitCopula(Copula.Fit.Object, 
                                pobs(CopulaFitData, ties.method = "average"), 
                                method = "ml",
                                optim.method = "BFGS", 
                                optim.control = list(maxit = 1000))
Copula.fit.Clayton@loglik # 92.16895

# copulaType - looks normal, acts normal
CopulaType <- "normal"

# Calculate copula image of the predictor
gamma.shape <- 4
gamma.rate <- 4
imgX <- sapply(dat$predictor,
               function(z)
                 (pgamma(z, gamma.shape, gamma.rate, lower.tail = TRUE, log.p = FALSE)))


# Estimate 5%, 95% and 50% image quantiles
theta <- Copula.fit.Gaussian@estimate
lowMidHigh <- as.data.frame(cbind(dat$predictor, dat$output))
colnames(lowMidHigh) <- c('predictor', 'output')

alpha<-.50
lowMidHigh$intensMid <- sapply(u_vector,
                               function(z)
                                 pnorm(qnorm(alpha)*sqrt(1-theta^2) + theta*qnorm(z)))
alpha<-.05
lowMidHigh$intensLow <- sapply(u_vector,
                               function(z)
                                 pnorm(qnorm(alpha)*sqrt(1-theta^2) + theta*qnorm(z)))
alpha<-.95
lowMidHigh$intensHigh <- sapply(u_vector,
                                function(z)
                                  pnorm(qnorm(alpha)*sqrt(1-theta^2) + theta*qnorm(z)))

# For Clayton
# alpha<-.50
# lowMidHigh$intensMid <- sapply(imgX,
#                              function(z)
#                                 (((alpha^(-theta/(1+theta)))-1)*(z^(-theta)) + 1)^(-1/theta))
# alpha<-.05
# lowMidHigh$intensLow <- sapply(imgX,
#                               function(z)
#                                 (((alpha^(-theta/(1+theta)))-1)*(z^(-theta)) + 1)^(-1/theta))
# alpha<-.95
# lowMidHigh$intensHigh <- sapply(imgX,
#                                function(z)
#                                  (((alpha^(-theta/(1+theta)))-1)*(z^(-theta)) + 1)^(-1/theta))

# For Frank

# alpha<-.50
# lowMidHigh$intensMid <- sapply(imgX,
#                               function(z)
#                                 -log(1-alpha*(1-exp(-theta))/(exp(-theta*z)+alpha*(1-exp(-theta*z))))/theta)
# alpha<-.05
# lowMidHigh$intensLow <- sapply(imgX,
#                               function(z)
#                                 -log(1-alpha*(1-exp(-theta))/(exp(-theta*z)+alpha*(1-exp(-theta*z))))/theta)
# alpha<-.95
# lowMidHigh$intensHigh <- sapply(imgX,
#                                function(z)
#                                  -log(1-alpha*(1-exp(-theta))/(exp(-theta*z)+alpha*(1-exp(-theta*z))))/theta)

# Reconstruct quantiles for the original data units from the image quantiles
normal.mean <- 90
normal.sd <- 2
lowMidHigh$intensMid <- sapply(lowMidHigh$intensMid,
                               function(z)
                                 (qlnorm(z, mean = normal.mean, sd = normal.sd, lower.tail = TRUE, log.p = FALSE)))
lowMidHigh$intensLow <- sapply(lowMidHigh$intensLow,
                               function(z)
                                 (qlnorm(z, mean = normal.mean, sd = normal.sd, lower.tail = TRUE, log.p = FALSE)))
lowMidHigh$intensHigh <- sapply(lowMidHigh$intensHigh,
                                function(z)
                                  (qlnorm(z, mean = normal.mean, sd = normal.sd, lower.tail = TRUE, log.p = FALSE)))


# Find anomalies
anomaly_intensity_idx <- (dat$output > lowMidHigh$intensHigh)|(dat$output < lowMidHigh$intensLow)

# plot
plot(dat$predictor, dat$output)
points(dat$predictor, lowMidHigh[, "intensLow"], col = "orange")
points(dat$predictor, lowMidHigh[, "intensHigh"], col = "orange")
points(dat$predictor, lowMidHigh[, "intensMid"], col = "cyan")
points(dat$predictor[anomaly_intensity_idx],
       dat$output[anomaly_intensity_idx], col = "magenta", pch = 16)
legend("topleft", legend = c("Median", "Tail", "Anomalies"), col = c("cyan", "orange", "magenta"), pch = c(1,1,16))

# set up output
quantileLow <- lowMidHigh$intensLow
quantileMid <- lowMidHigh$intensMid
quantileHigh <- lowMidHigh$intensHigh

# Count Regression: temperature-intensity data
poission.model <- glm(intensity ~ temperature, family=poisson(link='log'), data=part2_data)
summary(poission.model)$deviance
# 272.786, 248

nb.model <- MASS::glm.nb(intensity ~ temperature, data=part2_data)
summary(nb.model)$deviance
# 272.761, 248, 97524






