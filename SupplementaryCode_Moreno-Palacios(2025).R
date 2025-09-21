#'                                [Supplementary code]
#==========================================================================================
#                       TIME SERIES ANALYSIS (of avian breeding data)
#     TS DECOMPOSITION AND CONSINOR MODELS OF BREEDING SEASONALITY IN COLOMBIAN BIRDS
#==========================================================================================
#'                      *By Miguel Moreno-Palacios et al. (2025)*
#                            mc.morenop[at]uniandes.edu.co
#==========================================================================================

#==========================================================================================
# 1. IMPUTATION OF MISSING VALUES IN UNIVARIATE TIME SERIES
#==========================================================================================

# Load required libraries
library(imputeTS)   # For univariate time series imputation
library(zoo)        # For converting time series objects to zoo and using zoo functions
library(TSstudio)   # For plotting time series and converting between zoo and ts objects
library(forecast)   # For forecasting methods and additional imputation functions
library(dplyr)      # For data manipulation
library(tidyr)  # For reshaping data
library(lubridate) #dates
library(xts)    # For converting ts to data.frame via xts
library(fpp2)    # For forecasting and time series functions
library(TTR)     # For time series decomposition
library(season)  # For seasonal analysis tools

# Read data
# Please change the Data file for "COL_GBIF_TS.csv" to replicate GBIF analysis.
breed <- read.csv("COL_Banding_TS.csv", header = TRUE, sep = ";")

# Confirm total records per column (NA counts excluded)
colSums(breed, na.rm = TRUE)

# Convert to monthly time series starting January 2009
breeding.ts    <- ts(breed$breeding, frequency = 12, start = c(2009, 1))
nonbreeding.ts <- ts(breed$nonbreeding, frequency = 12, start = c(2009, 1))

# Plot missing-data patterns
ggplot_na_distribution(breeding.ts)
ggplot_na_distribution(nonbreeding.ts)

# Summary statistics for missing values
statsNA(breeding.ts)

#------------------------------------------------------------------------------------------
# Imputation methods
#------------------------------------------------------------------------------------------

# 1. Kalman smoothing and structural time series (imputeTS)
breeding.imp1    <- na_kalman(breeding.ts)
nonbreeding.imp1 <- na_kalman(nonbreeding.ts)

# 2. Seasonal decomposition imputation (imputeTS)
breeding.imp2    <- na_seadec(breeding.ts)
nonbreeding.imp2 <- na_seadec(nonbreeding.ts)

# 3. Seasonal split imputation (imputeTS)
breeding.imp3    <- na_seasplit(breeding.ts)
nonbreeding.imp3 <- na_seasplit(nonbreeding.ts)

# 4. Imputation by monthly averages (zoo)
#    Convert to zoo objects and replace NAs with the monthly mean
breeding.zoo      <- as.zoo(breeding.ts)
nonbreeding.zoo   <- as.zoo(nonbreeding.ts)
breeding.imp4.zoo <- na.aggregate(breeding.zoo, months)    # Replace NA with monthly mean
nonbreeding.imp4.zoo <- na.aggregate(nonbreeding.zoo, months)

# Convert back to ts objects
breeding.imp4    <- zoo_to_ts(breeding.imp4.zoo)
nonbreeding.imp4 <- zoo_to_ts(nonbreeding.imp4.zoo)

# 5. Last observation carried forward (na.locf, zoo)
#    Ensure first value is not NA by replacing with overall mean if necessary
if (is.na(breeding.zoo[1])) {
  breeding.zoo[1] <- mean(breeding.zoo, na.rm = TRUE)
}
breeding.imp5.zoo <- na_locf(breeding.zoo)

if (is.na(nonbreeding.zoo[1])) {
  nonbreeding.zoo[1] <- mean(nonbreeding.zoo, na.rm = TRUE)
}
nonbreeding.imp5.zoo <- na_locf(nonbreeding.zoo)

breeding.imp5    <- zoo_to_ts(breeding.imp5.zoo)
nonbreeding.imp5 <- zoo_to_ts(nonbreeding.imp5.zoo)

# 6. Structural time series imputation (na.StructTS, zoo)
#    Ensure first value is not NA
if (is.na(breeding.zoo[1])) {
  breeding.zoo[1] <- mean(breeding.zoo, na.rm = TRUE)
}
breeding.imp6.zoo <- na.StructTS(breeding.zoo)

if (is.na(nonbreeding.zoo[1])) {
  nonbreeding.zoo[1] <- mean(nonbreeding.zoo, na.rm = TRUE)
}
nonbreeding.imp6.zoo <- na.StructTS(nonbreeding.zoo)

breeding.imp6    <- zoo_to_ts(breeding.imp6.zoo)
nonbreeding.imp6 <- zoo_to_ts(nonbreeding.imp6.zoo)

# 7. Interpolation (na.interp, forecast)
breeding.imp7    <- na.interp(breeding.ts)
nonbreeding.imp7 <- na.interp(nonbreeding.ts)

# Plot original series along with one imputation method (you may change for comparison)
ggplot_na_imputations(breeding.ts, breeding.imp2, title = "Breeding (Imputed: Seadec)")
ggplot_na_imputations(nonbreeding.ts, nonbreeding.imp2, title = "Non-Breeding (Imputed: Seadec)")

# Combine all imputed versions into matrices
breed.imp <- cbind(
  breeding.imp1, breeding.imp2, breeding.imp3,
  breeding.imp4, breeding.imp5, breeding.imp6, breeding.imp7
)
colnames(breed.imp) <- c(
  "kalman", "seadec", "seasplit",
  "na.aggregate", "na.locf", "na.StructTS", "na.interp"
)

nonbreed.imp <- cbind(
  nonbreeding.imp1, nonbreeding.imp2, nonbreeding.imp3,
  nonbreeding.imp4, nonbreeding.imp5, nonbreeding.imp6, nonbreeding.imp7
)
colnames(nonbreed.imp) <- c(
  "kalman", "seadec", "seasplit",
  "na.aggregate", "na.locf", "na.StructTS", "na.interp"
)

#==========================================================================================
# 2. REORGANIZATION OF TABLES FOR TIME SERIES ANALYSIS
#==========================================================================================

# Adjust each imputed series by sample size to obtain a rate per 1,000 individuals
dynamic1.ts <- ((breeding.imp1    / (breeding.imp1    + nonbreeding.imp1))    * 1000)
dynamic2.ts <- ((breeding.imp2    / (breeding.imp2    + nonbreeding.imp2))    * 1000)
dynamic3.ts <- ((breeding.imp3    / (breeding.imp3    + nonbreeding.imp3))    * 1000)
dynamic4.ts <- ((breeding.imp4    / (breeding.imp4    + nonbreeding.imp4))    * 1000)
dynamic5.ts <- ((breeding.imp5    / (breeding.imp5    + nonbreeding.imp5))    * 1000)
dynamic6.ts <- ((breeding.imp6    / (breeding.imp6    + nonbreeding.imp6))    * 1000)
dynamic7.ts <- ((breeding.imp7    / (breeding.imp7    + nonbreeding.imp7))    * 1000)

# Combine all dynamics into one matrix
dynamics <- cbind(
  dynamic1.ts, dynamic2.ts, dynamic3.ts,
  dynamic4.ts, dynamic5.ts, dynamic6.ts, dynamic7.ts
)
colnames(dynamics) <- c(
  "kalman", "seadec", "seasplit",
  "na.aggregate", "na.locf", "na.StructTS", "na.interp"
)

#------------------------------------------------------------------------------------------
# 2.1 VISUAL INSPECTION OF ALL IMPUTATION-BASED DYNAMICS
#------------------------------------------------------------------------------------------

# Plot all series together to identify any negative values or anomalies
ts.plot(
  dynamic1.ts, dynamic2.ts, dynamic3.ts,
  dynamic4.ts, dynamic5.ts, dynamic6.ts, dynamic7.ts,
  gpars = list(
    col = c("black", "red", "blue", "orange", "purple", "green", "yellow"),
    xlab = "Year",
    ylab = "Breeding Rate (%)",
    lty = 1:7
  )
)

# Combine into a single ts.union object for multiple-panel plotting
dynamic_all.TS <- ts.union(
  dynamic1.ts, dynamic2.ts, dynamic3.ts,
  dynamic4.ts, dynamic5.ts, dynamic6.ts, dynamic7.ts
)
plot.ts(dynamic_all.TS, main = "All Time Series by Imputation Method")

#------------------------------------------------------------------------------------------
# 2.2 PREPARE TABLES AND ARRAYS FOR ANALYSIS
#------------------------------------------------------------------------------------------

# Convert dynamics matrix to xts and then to data.frame
dynamics.xts <- as.xts(dynamics)
dynamics.df  <- data.frame(
  date = index(dynamics.xts),
  coredata(dynamics.xts)
)

# Array 1: Add month index (1–12) for each row
dynamics.array1 <- dynamics.df %>%
  mutate(month = (row_number() - 1) %% 12 + 1)

# Calculate average of the valid (non-negative) imputation methods
# (here assuming methods in columns 4:6 are chosen as non-negative)
dynamics.array1$averageCorrected <- rowMeans(dynamics.array1[, c(4:6)], na.rm = TRUE)

# End of data reorganization

#==========================================================================================
# 3. TRADITIONAL TIME SERIES ANALYSIS
#==========================================================================================

# Extract the corrected average series (excluding negative-valued imputation methods) and convert to ts
dynamics.array1_av.ts <- ts(
  dynamics.array1$averageCorrected,
  frequency = 12,
  start = c(2009, 1)
)

# 3.1 Time plot of the corrected series
autoplot(dynamics.array1_av.ts) +
  ggtitle("Time Plot of Corrected Breeding Rate") +
  xlab("Year") +
  ylab("Monthly Rate of Breeding Birds") +
  theme_classic()

# 3.2 Ljung-Box test for overall randomness (white noise test)
Box.test(dynamics.array1_av.ts, lag = 24, fitdf = 0, type = "Ljung-Box")

# 3.3 Autocorrelation Function (ACF) plot
ggAcf(dynamics.array1_av.ts, lag = 120, type = "correlation") +
  ggtitle("Autocorrelation Function (ACF)") +
  theme_classic()

# 3.4 Periodogram (base R and with season package)
spec <- spec.pgram(dynamics.array1_av.ts, demean = TRUE, plot = TRUE)
plot(
  spec,
  ci = -1,
  main = "",
  xlab = "Frequency",
  ylab = "Periodogram",
  log = "no"
)
peri(dynamics.array1_av.ts)  # from season package

# 3.5 Cumulative periodogram to confirm periodicity signal
cpgram(dynamics.array1_av.ts, main = "")

# 3.6 Decomposition plot (multiplicative)
dynamics.array1_av.ts %>%
  decompose(type = "multiplicative") %>%
  autoplot() +
  xlab("Year") +
  ylab("Monthly Rate of Breeding Birds")

# Extract seasonal component from decomposition
dinamic.components <- decompose(dynamics.array1_av.ts, type = "multiplicative")
seasonal_data <- dinamic.components$seasonal

# Convert seasonal component into a wide data frame (years × months)
seasonal_df <- data.frame(
  Year     = floor(time(seasonal_data)),
  Month    = cycle(seasonal_data),
  Seasonal = as.numeric(seasonal_data)
)
seasonal_wide <- reshape(
  seasonal_df,
  timevar  = "Month",
  idvar    = "Year",
  direction = "wide"
)
colnames(seasonal_wide) <- c("Year", month.abb)

# 3.7 Seasonal plot (one line per year)
ggseasonplot(dynamics.array1_av.ts, year.labels = TRUE) +
  ggtitle("Seasonal Plot") +
  xlab("Month") +
  ylab("Monthly Rate of Breeding Birds") +
  theme_classic()

# 3.8 Seasonal subseries plot (show monthly variation and overall mean)
ggsubseriesplot(dynamics.array1_av.ts) +
  ggtitle("Seasonal Subseries Plot") +
  xlab("Time") +
  ylab("Monthly Rate of Breeding Birds") +
  theme_classic()

# 3.9 Box plot of month-wise distribution
boxplot(
  dynamics.array1_av.ts ~ cycle(dynamics.array1_av.ts),
  xlab = "Months",
  ylab = "Monthly Rate of Breeding Birds",
  main = ""
)

# 3.10 Smoothed trend with jittered points (monthly averages)
ggplot(dynamics.array1, aes(x = month, y = averageCorrected)) +
  geom_smooth(method = "loess", formula = y ~ x) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c(
      "Jan", "Feb", "Mar", "Apr", "May", "Jun",
      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    )
  ) +
  labs(
    x = "Month",
    y = "Corrected Avg. Breeding Birds",
    title = "Breeding Bird Dynamics in Colombia (Annual Variation)"
  ) +
  geom_jitter()

# End of traditional time series analysis

#==========================================================================================
# 4. SEASONAL ANALYSIS (season package)
#==========================================================================================

# Prepare data for seasonal analysis (corrected values and sample size)

# Convert breed.imp matrix to xts and then to data.frame
breed.imp.xts <- as.xts(breed.imp)
breed.imp.df  <- data.frame(date = index(breed.imp.xts), coredata(breed.imp.xts))

# Add month index (1–12)
breedimp2 <- breed.imp.df %>%
  mutate(month = (row_number() - 1) %% 12 + 1)

# Compute average of imputation methods (excluding those with negative values)
breedimp2$average <- rowMeans(breedimp2[, c(4:6)], na.rm = TRUE)

# Separate date column into month and year (date format was YYYY.MM)
breedS <- breedimp2 %>%
  separate(date, into = c("month2", "year"))

# Reconstruct a proper Date column (day set to 1st of month)
breedS$fecha <- as.Date(
  paste(as.numeric(breedS$month), "01", breedS$year, sep = "-"),
  format = "%m-%d-%Y"
)

# Number of days in each month
breedS$ndaysmonth <- lubridate::days_in_month(breedS$fecha)

# Convert columns to correct types
breedS$month <- as.integer(breedS$month)
breedS$year  <- as.numeric(breedS$year)
breedS$average <- as.integer(breedS$average)

# Compute year-month decimal variable
breedS$yrmon <- breedS$year + (breedS$month - 1) / 12

# Sample size (mean of selected imputation methods for breeding and nonbreeding)
breedS$captures <- rowMeans(breed.imp[, c(3:5)]) +
  rowMeans(nonbreed.imp[, c(3:5)])
breedS$captures <- as.integer(breedS$captures)

# Adjusted breeding count (daily adjustment based on days in month)
breedS$adjb <- breedS$average * ((365.25 / 12) / breedS$ndaysmonth)

# Average corrected series from dynamics.array1
breedS$averageCorrected   <- dynamics.array1$averageCorrected
breedS$averageCorrected2  <- as.integer(round(dynamics.array1$averageCorrected))

# Corrected adjusted breeding count (using rounded corrected values)
breedS$adjbC <- breedS$averageCorrected2 * ((365.25 / 12) / breedS$ndaysmonth)

breedS <- as.data.frame(breedS)

#------------------------------------------------------------------------------------------
# 4.1 Prepare matrix for seasonal analysis (monthmean function)
#------------------------------------------------------------------------------------------

# Compute monthly means for seasonal analysis (using 'averageCorrected2' as response)
mmean <- monthmean(
  data       = breedS,
  resp       = 'averageCorrected2',
  offsetpop  = NULL,
  adjmonth   = 'average'
)
print(mmean)
plot(mmean)

#------------------------------------------------------------------------------------------
# 4.2 Linear regression model (GLM) with effect of month
#------------------------------------------------------------------------------------------

lmmodel <- monthglm(
  formula     = averageCorrected2 ~ 1,
  data        = breedS,
  family      = poisson(),
  offsetpop   = NULL,
  offsetmonth = TRUE,
  refmonth    = 9
)
plot(lmmodel)
summary(lmmodel)

# Check residuals for seasonality
seasrescheck(resid(lmmodel))
shapiro.test(lmmodel$residuals)

#------------------------------------------------------------------------------------------
# 4.3 Cosinor models (stationary)
#------------------------------------------------------------------------------------------

# Single-cycle cosinor model
res <- cosinor(
  averageCorrected2 ~ 1,
  date     = "month",
  data     = breedS,
  type     = 'monthly',
  family   = poisson(),
  cycles   = 1,
  offsetmonth = TRUE,
  offsetpop   = NULL
)
summary(res)
plot(res)

# Two-cycle cosinor model
res2 <- cosinor(
  averageCorrected2 ~ 1,
  date     = "month",
  data     = breedS,
  type     = 'monthly',
  family   = poisson(),
  cycles   = 2,
  offsetmonth = TRUE,
  offsetpop   = NULL
)
summary(res2)
plot(res2)

# Check GLM summaries, AIC, and residual diagnostics
summary(res$glm)
summary(res2$glm)
seasrescheck(resid(res))
seasrescheck(resid(res2))
shapiro.test(res$residuals)

#------------------------------------------------------------------------------------------
# 4.4 Non-stationary cosinor models
#------------------------------------------------------------------------------------------

nsmodel <- nscosinor(
  data      = breedS,
  response  = "adjbC",
  cycles    = c(12),
  niters    = 5000,
  burnin    = 1000,
  tau       = c(5, 100)
)
summary(nsmodel)
plot(nsmodel, main = "Non-Stationary Cosinor (cycles = 12)")
seasrescheck(resid(nsmodel))
shapiro.test(resid(nsmodel))
qqnorm(resid(nsmodel))

nsmodel2 <- nscosinor(
  data      = breedS,
  response  = "adjbC",
  cycles    = c(6),
  niters    = 5000,
  burnin    = 1000,
  tau       = c(5, 100)
)
summary(nsmodel2)
plot(nsmodel2, main = "Non-Stationary Cosinor (cycles = 6)")
seasrescheck(resid(nsmodel2))
shapiro.test(resid(nsmodel2))

#==========================================================================================
# END OF SCRIPT
#==========================================================================================