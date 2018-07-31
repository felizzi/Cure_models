################################
# Plot external mortality data #
################################

#################################################################
# Functions to plot demographic characteristics for the trial and
#  background mortality (using HMD life tables)
#################################################################

# Load packages
library(graphics)

# Plot trial data

## Extract data
### For sex
t_Sex <- table(trial_data$SEX)
t_Sex <- t_Sex / sum(t_Sex)
### For countries
t_country <- table(trial_data$COUNTRY)

## Define colors
jet.colors <- colorRampPalette(c("lightcyan", "navyblue"))

## Open graphics device and define layout
jpeg(paste0("reports/", "countries_years.jpg", sep = ""),
     width = 6, height = 5, units = "in", res = 500)
layout(matrix(c(1, 2, 2, 3, 3, 3), 2, 3, byrow = TRUE))

### Barplot for gender distribution
par(mar = c(4, 5, 2, 2))
barplot(t_Sex, main = "Gender", col = c("red", "blue"),
        cex.lab = 1.5, cex.axis = 1.2, cex = 1.2, ylab = "Percent",
        axes = FALSE, ylim = c(0, .7))
axis(2, at = c(0, 0.33, 0.6), label = c("0", "33", "60"),
     cex.lab = 1.5, cex = 1.2, cex.axis = 1.2)

### Histogram for age distribution
par(mar = c(5, 6, 2, 2))
rh <- hist(as.numeric(trial_data$AGE), main = "Age",
           cex.axis = 1.2, cex.lab = 1.5, cex = 1.2, freq = TRUE, axes = TRUE,
           xlab = "Age (years) at start of trial")

### Barplot for geographical distribution
barplot(t_country, las = 2, main = "Countries",
        ylab = "Number of patients", cex.names = 1.5, cex.lab = 1.5)

## Close graphics device
dev.off()

# Plot background survival data, for Italy, Russia and the US (as examples)

## Extract data
### For Italy
table_in_italy <- (subset(mortality_all_latest, Country_code == "ITA"))
table_in_italy_M <- subset(table_in_italy, Gender == "MALE")
table_in_italy_F <- subset(table_in_italy, Gender == "FEMALE")
### For Russia
table_in_russia <- (subset(mortality_all_latest, Country_code == "RUS"))
table_in_russia_M <- subset(table_in_russia, Gender == "MALE")
table_in_russia_F <- subset(table_in_russia, Gender == "FEMALE")

## Open graphics device
jpeg(paste(wd, "plots_countries.jpg", sep = "/"), width = 10, height = 7.5, units = "in", res = 500)
par(mar = c(4.5, 4.5, 2, 2))
layout(matrix(c(1, 2), 2, 1, byrow = TRUE))

## Plot survival curves
### For men
plot(table_in_italy_M$lx / 1e5, lty = 1, lwd = 2, type = "l", xlab = "Age",
     ylab = "Survival probability", cex.lab = 1.5, cex.axis = 1.2,
     main = "Survival Male")
lines(table_in_russia_M$lx / 1e5, col = "red", lwd = 2)
lines(USA_table_male_2013$lx / 1e5, col = "blue", lwd = 2)
legend("bottomleft",c("Italy 2009", "Russia 2010", "USA 2013"), cex = 1.2,
       col = c("black", "red", "blue"), lty = c(1, 1, 1), lwd = c(2, 2, 2))

### For women
par(mar = c(4.5, 4.5, 2, 2))
plot(table_in_italy_F$lx / 1e5, lty = 1, lwd = 2, type = "l", xlab = "Age",
     ylab = "Survival probability", cex.lab = 1.5, cex.axis = 1.2,
     main = "Survival Female")
lines(table_in_russia_F$lx / 1e5, col = "red", lwd = 2)
lines(USA_table_female_2013$lx / 1e5, col = "blue", lwd = 2)

## Close graphics device
dev.off()
