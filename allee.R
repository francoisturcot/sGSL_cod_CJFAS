library(tidyverse)
library(segmented)

#read SCA data
SCA_data <- readRDS("SCA_data.RDS")

catch = SCA_data$catch
ssb = SCA_data$ssb
b = SCA_data$b2to11
year = SCA_data$year

#calculate production
p = b #prepare vector

for (t in 1:length(1950:2017)){
    
    p[t] = catch[t] + b[t+1] - b[t]
    
}

#build df
dat = data.frame(year = year, biomass = b, ssb = ssb, production = p, production.rate = p/b)
dat = dat %>% filter(year<2018)

#find a breakpoint using peicewise regression
model = lm(production.rate~biomass, data = dat)
segmented <- segmented(model, seg.Z = ~biomass)
summary(segmented)

# Plot the original data with the fitted model
seg_preds <- predict(segmented)
seg_res <- dat$production.rate - seg_preds

plot(
    dat$biomass, dat$production.rate,
    main = "Piecewise Regression Fit",
    xlab = "Independent Variable (x)",
    ylab = "Dependent Variable (y)",
    col = "blue"
)
lines(dat$biomass, seg_preds,col = "red", lwd = 2)

#linear regressions for each period
dat1 = dat %>% filter(year<1992)
m1 = lm(production.rate~biomass, data = dat1)
summary(m1)
plot(m1)

dat2 = dat %>% filter(year>1991)
m2 = lm(production.rate~biomass, data = dat2)
summary(m2)
plot(m2)

