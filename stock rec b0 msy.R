library(tidyverse)
library(FSA)
library(reshape2)
library(fishmethods)

#read data 
data = readRDS("SCA_data.RDS")

#lag recruitment by 2 years
df = data.frame(stock = data$ssb[1:67]/1000,recruits = data$recruits[3:69]/1000000)
df$logR = log(df$recruits)

#starting values using FSA
bh1s <- srStarts(recruits~stock,data=df,type="BevertonHolt",param=1)
unlist(bh1s) 

# Fit BH model with lognormal error
beverton_holt_lognorm <- nls(
    logR ~ log((a * stock) / (1 + b * stock)),
    data = df,
    start=bh1s
)

summary(beverton_holt_lognorm)

# Plot residuals vs. fitted
resid <- residuals(beverton_holt_lognorm)
plot(fitted(beverton_holt_lognorm), resid, xlab = "Fitted recruitment", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

#SSB at 50%Rmax
S_50_lognorm <- 1 / coef(beverton_holt_lognorm)["b"]
S_50_lognorm

# Fit Ricker model with lognormal error
ricker_lognorm <- nls(
    logR ~ log(a) + log(stock) - b * stock,
    data = df,
    start=list(a = 1, b = 0.001)
)
summary(ricker_lognorm)

# Calculate residuals
residuals_ricker <- residuals(ricker_lognorm)

# Plot residuals vs fitted values
plot(fitted(ricker_lognorm), residuals_ricker,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

# Get coefficients
b_lognorm <- coef(ricker_lognorm)["b"]
a_lognorm <- coef(ricker_lognorm)["a"]

#find SSB at 50%Rmax
SSB_at_max <- 1 / b_lognorm
R_max_ricker <- a_lognorm * SSB_at_max * exp(-1)
RK50R <- 0.5 * R_max_ricker
ricker_fn <- function(S) {
    a_lognorm * S * exp(-b_lognorm * S) - RK50R
}
RK50 <- uniroot(ricker_fn, lower = 1, upper = SSB_at_max)$root
RK50

#replacement

# survivorship function 
survivorship_F <- function(f = 0, M, n_ages, sel, message = TRUE) {
    l_age <- rep(NA, n_ages)
    l_age[1] <- 1
    for (a in 2:(n_ages - 1)) {
        l_age[a] <- l_age[a - 1] * exp(-(M[a - 1] + f * sel[a - 1]))
    }
    l_age[n_ages] <- l_age[n_ages - 1] * exp(-(M[a - 1] + f * sel[a - 1])) /
        (1 - exp(-(M[a] + f * sel[a])))
    return(l_age)
}

# Prepare data
sel <- data$sel
M <- data$M
n_ages <- length(2:12)
waa = data$waa
mat = data$mat

# Compute survivorship
surv <- t(apply(cbind(1:69), 1, function(i) {
    survivorship_F(f = 0, M = M[i, ], n_ages = n_ages, sel = sel[i, ], message = FALSE)
}))

#annual phi0
phi0 = (rowSums(surv*waa*mat))*1000
phi0

# SRR with 1/phi0 lines
# Extract BH coefficients
a <- coef(beverton_holt_lognorm)["a"]
b <- coef(beverton_holt_lognorm)["b"]

# Create a sequence of stock values
stock_seq <- seq(0, max(df$stock*5), length.out = 200)

# Compute predicted log(R), then transform to R
pred_logR <- log((a * stock_seq) / (1 + b * stock_seq))
pred_R <- exp(pred_logR)

# Define plot limits
x_max <- max(df$stock) * 1.2
y_max <- max(df$recruits) * 1.2

# Initial plot without y-axis
plot(df$recruits ~ df$stock, 
     xaxs = "i", yaxs = "i", 
     xlim = c(0, x_max), 
     ylim = c(0, y_max),
     xlab = "SSB (kt)",
     ylab = "Recruits (millions)")

# Draw all lines and highlight groupings as you already did
for(i in 1:39) { abline(a = 0, b = (1/phi0[i]), col = "grey") }
for(i in 40:44) { abline(a = 0, b = (1/phi0[i]), col = "red") }
for(i in 45:53) { abline(a = 0, b = (1/phi0[i]), col = "grey") }
for(i in 54:69) { abline(a = 0, b = (1/phi0[i]), col = "red") }

# Add BH curve 
lines(stock_seq, pred_R, lty = 5, lwd = 2)


# Extract BH coefficients
a <- coef(beverton_holt_lognorm)["a"]
b <- coef(beverton_holt_lognorm)["b"]

# Create a sequence of stock values and predicted recruitment
stock_seq <- seq(0, max(df$stock) * 1.2, length.out = 500)
pred_R <- (a * stock_seq) / (1 + b * stock_seq)

# Function to check intersection for a given slope (1/phi0)
check_intersection <- function(slope) {
    # Recruitment line: R = slope * SSB
    # Find difference between BH curve and line
    diff <- pred_R - slope * stock_seq
    # Check if sign changes (indicating intersection)
    any(diff[-1] * diff[-length(diff)] < 0)
}

# Apply to all phi0 values
intersections <- sapply(1/phi0, check_intersection)

year = 1950:2018

year[intersections]
#1950 to 1985 dont intersect but below the curve, not over
#first year post 1985 that dont intersect is 1990



## ------------------ MSY
#Dynamic

#fit BH without scaling SSB and recruits
#lag recruitment by 2 years
df = data.frame(stock = data$ssb[1:67],recruits = data$recruits[3:69])
df$logR = log(df$recruits)

#starting values using FSA
bh1s <- srStarts(recruits~stock,data=df,type="BevertonHolt",param=1)
unlist(bh1s) 

# Fit BH model with lognormal error
beverton_holt_lognorm <- nls(
    logR ~ log((a * stock) / (1 + b * stock)),
    data = df,
    start=bh1s
)

# Extract BH coefficients
a <- coef(beverton_holt_lognorm)["a"]
b <- coef(beverton_holt_lognorm)["b"]

#prepare data
ny = length(year)
Fmsy <- rep(NA, ny)
msy <- rep(NA, ny)
SSBmsy <- rep(NA,ny)
ypr = rep(NA,ny)
ypr <- matrix(rep(ypr, each = 501), nrow = 501, ncol = ny, byrow = TRUE)
eq_rec_f = rep(NA,ny)
eq_rec_f <- matrix(rep(eq_rec_f, each = 501), nrow = 501, ncol = ny, byrow = TRUE)
yield = rep(NA,ny)
yield <- matrix(rep(yield, each = 501), nrow = 501, ncol = ny, byrow = TRUE)
f <- seq(0,5,0.01)
phi_f <- rep(NA,length(f))
YPR_a <- SURV_a <- matrix(NA,ncol=length(waa[1,]),nrow=length(f))
rownames(SURV_a)<-rownames(YPR_a)<-f

for (k in 1:ny){
    
    M_i = M[k,]
    waa_i = waa[k,]
    mat_i = mat[k,]
    sel_i = sel[k,]
    
    for(i in 1:length(f)){
        SURV_a[i,] <- survivorship_F(f=f[i],M=M_i,n_ages=length(sel_i),sel=sel_i)
        for(j in 1:ncol(YPR_a)){
            YPR_a[i,j] <- SURV_a[i,j]*waa_i[j]*(1-exp(-(M_i[j]+f[i]*sel_i[j])))*f[i]*sel_i[j]/(M_i[j]+f[i]*sel_i[j])
        }
        phi_f[i] <- sum(SURV_a[i,]*waa_i*mat_i)
    }  
    
    ypr[,k] <- rowSums(YPR_a)
    eq_rec_f[,k] <- (1/b*(a-1/phi_f))
    eq_rec_f[eq_rec_f<0] <- 0
    yield[,k] <- ypr[,k]*eq_rec_f[,k]
    
    Fmsy[k] <- f[which(yield[,k]==max(yield[,k]))]
    msy[k] <- yield[,k][which(yield[,k]==max(yield[,k]))]
    SSBmsy[k] <- eq_rec_f[,k][which(yield[,k]==max(yield[,k]))]/(a-eq_rec_f[which(yield[,k]==max(yield[,k]))]*b)
}

#static
static_bmy_lrp = mean(SSBmsy[1:6])*0.4


#############----------B0
#dynamic
#get r0 from BH

r0 = (1 / b_lognorm) * (a_lognorm - (1/phi0)) 
r0[r0<0] = 1

#get b0
b0 = r0*phi0*1000 #tonnes
b0

#static
#use proxy for eq recruitment (1981 to 1984 recruitment)
df$year = 1950:2016
r = df %>% filter(year %in% 1981:1984)
eq_rec = mean(r$rec)
static_b0 = 0.25*mean(phi0[1:6])*eq_rec/1000

#plot dynamic b0 and msy

plot = data.frame(year = 1950:2018, 
                  bmsy = SSBmsy*0.4, 
                  b0 = 0.25*b0, 
                  ssb = data$ssb, 
                  static_bmsy = static_bmy_lrp, 
                  static_b0 = static_b0)

plot = melt(plot, id.vars = "year")

str(plot)

ggplot(plot, aes(year, value,
                 colour = variable,
                 linetype = variable,
                 alpha = variable)) +          # <-- add alpha mapping
    geom_line() +
    theme_bw() +
    xlab("") +
    ylab("SSB (tonnes)") +
    
    # --- COLORS: shared within pairs ---
    scale_colour_manual(
        name = "",
        values = c(
            "bmsy"        = "lightcoral",
            "static_bmsy" = "lightcoral",
            "b0"          = "lightblue",
            "static_b0"   = "lightblue",
            "ssb"         = "black"
        ),
        labels = c(
            expression(paste("Dynamic ", 0.4 * B[MSY])),
            expression(paste("Dynamic ", 0.25 * B[0])),
            "SSB",
            expression(paste("Static ", 0.4 * B[MSY])),
            expression(paste("Static ", 0.25 * B[0]))
        )
    ) +
    
    # --- LINETYPES ---
    scale_linetype_manual(
        name = "",
        values = c(
            "bmsy"        = "solid",
            "static_bmsy" = "dashed",
            "b0"          = "solid",
            "static_b0"   = "dashed",
            "ssb"         = "solid"
        ),
        labels = c(
            expression(paste("Dynamic ", 0.4 * B[MSY])),
            expression(paste("Dynamic ", 0.25 * B[0])),
            "SSB",
            expression(paste("Static ", 0.4 * B[MSY])),
            expression(paste("Static ", 0.25 * B[0]))
        )
    ) +
    
    # --- ALPHA: make all non-SSB lines half-transparent ---
    scale_alpha_manual(
        values = c(
            "bmsy"        = 0.75,
            "static_bmsy" = 0.75,
            "b0"          = 0.75,
            "static_b0"   = 0.75,
            "ssb"         = 1     # fully opaque
        ),
        guide = "none"           # hides alpha from legend
    )




# Proxy for BMSY
# (1) Bmsy proxy = The biomass corresponding to the biomass per recruit at F0.1 multiplied 
# by the average number of recruits.

ages = 2:12
#means for 1950 to 1955
m = colMeans(M[1:6,])
w = colMeans(waa[1:6,])
mataa = colMeans(mat[1:6,])
s = colMeans(sel[1:6,])

#find F0.1
y = ypr(age=ages, wgt = w, partial = s, M = m, plus=TRUE, oldest = 20, maxF = 2, incrF=0.1, graph=TRUE)
y$Reference_Points
f0.1 = y$Reference_Points[1]
f0.1

#find ssb per recruit at F0.1
l_age0.1 <- survivorship_F(f=f0.1,M=m,n_ages=length(ages),sel=s)
l_age0.1
phi0.1 <- sum(l_age0.1 *w*mataa) 

#use equilibrium recruitment
msyproxy1 = phi0.1*eq_rec
msyproxy1
0.4*msyproxy1 

# (2) The average biomass (or index of biomass) over a productive period.

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
dat = data.frame(year = year, biomass = b, production = p)
dat = dat %>% filter(year<2018)

#plot
plot(scale(dat$b), type = "l", ylim = c(-1.5, 3.25), ylab = "Scaled value", xlab = "")
lines(scale(dat$p), col = "red")

# (3) The biomass corresponding to 50% of the maximum historical biomass
max(ssb)
0.5*max(ssb)
0.5*max(ssb)*0.4
