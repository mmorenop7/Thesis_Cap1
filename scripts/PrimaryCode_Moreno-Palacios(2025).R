#'                                [Primary code]
#==========================================================================================
#                       TIME SERIES ANALYSIS (of avian breeding data)
#         GAM MODELS AND CIRCULAR STATISTICS OF BREEDING SEASONALITY IN COLOMBIAN BIRDS
#==========================================================================================
#'                      *By Miguel Moreno-Palacios et al. (2025)*
#                            mc.morenop[at]uniandes.edu.co
#==========================================================================================
# libraries
library(lubridate)
library(mgcv) 
library(ggplot2) 
library(gratia) 
library(dplyr)
library(visreg)
library(circular)
library(imputeTS)

# Breeding records matrix (this case from excel)
data <- read.csv("Breeding_database.csv", head = TRUE, sep = ",")

#==========================================================================================
# [1] ADJUST VARIABLES
#==========================================================================================
data$elevationCat = if_else(data$Elevation <= 999, "Low elevations",
                     if_else(data$Elevation > 1000 & data$Elevation <= 1999, "Mid elevations",
                     if_else(data$Elevation > 2000, "High elevations", false=NA)))

data$breedingCor <- as.numeric(data$breedingCor)
data$region <- as.factor(data$region)
data$species <- as.factor(data$species)
data$guild <- as.factor(data$guild)
data$latitude <- as.numeric(data$latitude)
data$Elevation <- as.numeric(data$Elevation)
data$elevationCat <- as.factor(data$elevationCat)
data$database <- as.factor(data$database)

# FIND DE DATE AND DOY
data$date <- dmy(paste(data$day,data$month,data$year))
data$date <- as.Date(as.character(data$date), format = "%Y-%m-%d")
data$doy <- yday(data$date)


#==========================================================================================
# [2] GENERALIZED ADDITIVE MODELS
#==========================================================================================
# [2.1] Basic GAM
GAM_basic <- gam(breedingCor ~ s(doy, bs="cc", k=-1)
                        , family="binomial", data=data[data$status=="Native",], method = "REML")
summary(GAM_basic)

# [2.2] GAM with covariates
data_natives <- data[data$status == "Native", ]
GAM_covariates <- gam(breedingCor ~ s(doy, bs="cc", k=-1) + guild + s(year, bs="re") + Elevation + latitude
                    , family = "binomial", data=data_natives, method = "REML")
summary(GAM_covariates)

# [2.3] Elevation effect 
GAM_Elevation <- gam(breedingCor ~ s(doy, by=elevationCat, bs="cc", k=-1) + guild + s(year, bs="re") + Elevation + latitude
                     , family = "binomial", data=data[data$status=="Native",], method = "REML")
summary(GAM_Elevation)
draw(GAM_Elevation)
# [2.4] Regional effect
GAM_region <- gam(breedingCor ~ s(doy, by=region, bs="cc", k=-1) + guild + s(year, bs="re") + Elevation + latitude
                 , family = "binomial", data=data[data$status=="Native",], method = "REML")
summary(GAM_region)

# [2.5] Guild effect
GAM_guilds <- gam(breedingCor ~ s(doy, by=guild, bs="cc", k=-1) + s(year, bs="re") + Elevation + latitude
                 , family = "binomial", data=data[data$status=="Native",], method = "REML")
summary(GAM_guilds)

#==========================================================================================
# [2.6] DATABASE difference
#==========================================================================================

data_banding <- data %>% filter(database == "banding")
data_gbif <- data %>% filter(database == "GBIF")



#GAM
GAM_banding <- gam(breedingCor ~ s(doy, bs="cc", k=-1)
                   , family="binomial", data=data_banding[data_banding$status=="Native",], method = "REML")
GAM_gbif <- gam(breedingCor ~ s(doy, bs="cc", k=-1)
                , family="binomial", data=data_gbif[data_gbif$status=="Native",], method = "REML")

#==========================================================================================
# [3] CIRCULAR ANALYSIS
#==========================================================================================

# df with doy and angles
IndexAngular <- data.frame(
  numday = 1:365,
  angularday = (1:365) * (360 / 365)
)

# [3.1] identify NA values and use only presence (breeding = 1)
sum(is.na(data$breedingCor))  # Ver cuántos NA hay
data$breedingCor <- na_replace(data$breedingCor, fill = 0)
# assign angles only to "presence values"
index <- IndexAngular$numday
values <- IndexAngular$angularday
# only active breeding records
data_breeding <- data %>% filter(breedingCor == 1)  
data_breeding$doyAngular <- values[match(data_breeding$doy, index)]
# Remove dates with no angles
data_breeding <- data_breeding[!is.na(data_breeding$doyAngular), ]

# [3.2] filter database
banding_dataset <- data_breeding %>% filter(database == "banding")
gbif_dataset <- data_breeding %>% filter(database == "GBIF")

#[3.3] Circular analysis
circular_banding <- circular(banding_dataset$doyAngular, 
                             type="angles", units="degrees", zero = pi/2, rotation ="clock")
circular_gbif <- circular(gbif_dataset$doyAngular, 
                          type="angles", units="degrees", zero = pi/2, rotation ="clock")
median.circular(circular_banding , na.rm=F)
mean.circular(circular_banding , na.rm=F)
# basic circular distribution tests
watson.test(circular_banding, alpha = 0.05, dist ="uniform") # null H0
rao.spacing.test(circular_banding) # null H0
watson.test(circular_banding, alpha = 0, dist ="vonmises") # unimodal (normal circular)

#==========================================================================================
#' *FIGURES*
#==========================================================================================

# [1] Figure 2A
# smooths
F2A <- smooth_estimates(GAM_covariates) %>% 
  filter(.smooth == "s(doy)") %>%
  mutate(
    prob = plogis(.estimate), 
    lower = plogis(.estimate - 2 * .se),
    upper = plogis(.estimate + 2 * .se),
    model = "GAM_covariates"
  )
#plot
fig2A <- ggplot(F2A, aes(x = doy, y = prob)) +
  geom_line(linewidth = 0.5, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightgray") +
  labs(x = "Day of year", y = "Probability of Breeding activity") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# [2] Figures 2B (Latitude) and 2C (elevation) 
# elevation
plot_elev <- visreg(GAM_covariates, "Elevation", scale = "response", gg = TRUE)
fig2B <- plot_elev + labs(y = "Breeding activity", x = "Elevation") +
  ylim(0.10, 0.25) +
  theme_classic(base_size = 14) +  
  theme(panel.grid.minor = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        )
# latitude
plot_lat <- visreg(GAM_covariates, "latitude", scale = "response", gg = TRUE)
fig2C <- plot_lat + labs(y = "Breeding activity", x = "Latitude") +
  ylim(0.10, 0.25) +
  theme_classic(base_size = 14) +  
  theme(panel.grid.minor = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
        )

# [3] Figure 2D
#smooths
s1 <- smooth_estimates(GAM_banding) %>% 
  filter(.smooth == "s(doy)") %>%
  mutate(
    prob = plogis(.estimate), 
    lower = plogis(.estimate - 2 * .se),
    upper = plogis(.estimate + 2 * .se),
    model = "Banding"
  )

s2 <- smooth_estimates(GAM_gbif) %>% 
  filter(.smooth == "s(doy)") %>%
  mutate(
    prob = plogis(.estimate), 
    lower = plogis(.estimate - 2 * .se),
    upper = plogis(.estimate + 2 * .se),
    model = "GBIF"
  )
all_data <- bind_rows(s1, s2)

#plot
fig2D <- ggplot(all_data, aes(x = doy, y = prob, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(x = "Day of year", y = "Probability of Breeding activity", color = "Data source", fill = "Data source") +
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# [4] Figure 2E
png("/Users/miguelmoreno/Documents/R/phd_thesis/Breeding/CircPlot.png", width = 2000, height = 1650, res = 300)
#--
plot.circular(circular_banding, cex=0.3, sep=0.025, stack=F, col = "dark gray", tol = 0.01, axes = FALSE, 
              ticks = F, shrink=1.1, main = "")
rose.diag(circular_gbif, axes = F, pch = 2, ticks=F, border="lightblue", cex = 0.3, bins = 365, prop = 9.0, shrink=0.5,
          add = T, col = "darkblue", zero = pi/2.0, rotation ="clock")
rose.diag(circular_banding, axes = F, pch = 2, ticks=F, border="lightpink", cex = 0.7, bins = 365, prop = 6.0, shrink=0.5,
          add = T, col = "darkred", zero = pi/2.0, rotation ="clock")
axis.circular(at=circular(seq(0, 11*pi/6,pi/6)), 
              labels=c("A","M","F","J","D","N","O","S","A", "J", "J", "M"), cex=1.0)
# kernel density 
lines(density.circular(circular_banding, bw=50), lwd=2, lty=9, col="red")
lines(density.circular(circular_gbif, bw=50), lwd=2, lty=9, col="blue")
# vectors
arrows.circular(mean(circular_banding),lwd=1, col = "red")
arrows.circular(mean(circular_gbif),lwd=1, col = "blue")
# Legend
legend("topright", legend = c("Banding", "GBIF"), col = c("red", "blue"), lty = c(9, 9), lwd = 2, bty = "n")
#--
dev.off()

# [5] Figure 3A
data_filtered_elev <- data %>% filter(elevationCat != "NA")
fig3A <- ggplot(data = data_filtered_elev, aes(x = doy, y = breedingCor, color = factor(elevationCat))) +
  geom_smooth(method ="gam", 
              formula= y ~ s(x, bs = "cc", k=-1),
              method.args=list(family="binomial")) +
  coord_cartesian(ylim = c(0, 0.4)) +
  xlab("Day of year") + ylab("Breeding activity") +
  labs(color = "Elevation") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.8), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()
        )

# [6] Figure 3B
data_filtered_region <- data %>% filter(region != "Unknown")
fig3B <- ggplot(data = data_filtered_region, aes(x = doy, y = breedingCor, color = factor(region))) +
  geom_smooth(method ="gam", 
              formula= y ~ s(x, bs = "cc", k=-1),
              method.args=list(family="binomial"), se = TRUE) +
  coord_cartesian(ylim = c(0, 0.4)) + 
  xlab("Day of year")+ ylab("Breeding activity")+
  labs(color = "Region") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.75), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()
        )

# [7] Figure 3C
data_filtered_guild <- data %>% filter(guild != "Unknown")
fig3C <- ggplot(data = data_filtered_guild, aes(x = doy, y = breedingCor, color = factor(guild))) +
  geom_smooth(method = "gam", 
              formula = y ~ s(x, bs = "cc", k = -1),
              method.args = list(family = "binomial")) +  
  coord_cartesian(ylim = c(0, 0.4)) +
  xlab("Day of year") + ylab("Breeding activity") +
  labs(color = "Guild") + 
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.8),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank())

#==========================================================================================
# END OF SCRIPT
#==========================================================================================
