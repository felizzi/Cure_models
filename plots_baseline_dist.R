## plot for external mortality data
trial_data_1 <- read.csv(paste0("data/example_1.csv"))
trial_data <- read.csv(paste0("data/example.csv"))

require(graphics)
jet.colors <-
  colorRampPalette(c("lightcyan", "navyblue"))
jpeg(paste0("reports/","countries_years.jpg",sep = ""), width=6, height=5, units="in", res=500)
layout(matrix(c(1,2,2,3,3,3), 2, 3, byrow = TRUE))
t_Sex <- table(trial_data$SEX)
t_Sex <- t_Sex/sum(t_Sex)
par(mar=c(4,5,2,2))  
barplot(t_Sex, main = 'Gender', col = c( 'red', 'blue'), cex.lab = 1.5, cex.axis = 1.2, cex = 1.2, ylab = "Percent", axes = FALSE, ylim = c(0,.7))
axis(2, at = c(0,0.33, 0.6), label = c("0%","33%", "60%"), cex.lab = 1.5, cex = 1.2, cex.axis = 1.2)
#hist(as.numeric(seer1$Age), main = "Age", cex.axis = 1.2, cex.lab = 1.5, xlab = "Age", cex = 1.2, freq = TRUE)
par(mar = c(5,6,2,2) )
rh <- hist(as.numeric(trial_data$AGE), main = "Age", cex.axis = 1.2, cex.lab = 1.5,  cex = 1.2, freq = TRUE, axes = TRUE, xlab = "Age at start of trial")
#axis(1, at = c(20,30,50,60,70,80,90), cex.lab = 1.5, cex.axis = 1.9)
#axis(2,)
t_country <- table(trial_data$COUNTRY)
barplot(t_country, las = 2, col = jet.colors(19), main = "Countries", ylab = "Frequency", cex.names = 1.5, cex.lab = 1.5)
dev.off()


axis(2, at = NULL, ylab = "Frequency", cex.lab = 1.5, cex = 1.2, las = 1)

#barplot(tCoBRIM, las = 2, col = jet.colors(19), main = "Countries",  cex.lab = 1.5, cex = 1.2)

#axis(2, at = c(0, 250, 500), cex.lab = 1.5, cex.axis = 1.2)
#plot(0:110, table_USA_Male$lx[id1988m]/1e5, xlab = "Age", ylab = "Surv. prob.", cex.lab = 1.5, cex.axis = 1.2, col = 'blue', lwd = 2, lty  = 1, type = "l", main = "Background Mortality")
#lines(0:110, table_USA_Male$lx[id2008m]/1e5, col = 'blue', lwd = 2, lty = 2)
#lines(0:110, table_USA_Female$lx[id1988f]/1e5, col = 'red', lwd = 2, lty = 1, type = "l")
#lines(0:110, table_USA_Female$lx[id2008f]/1e5, col = 'red', lty = 2, cex.axis = 1.2, type = "l", lwd = 2)
#legend("bottomleft", c("Males 1988", "Males 2008","Females 1988", "Females 2008"), cex=1.2, col=c("blue", "blue", "red", "red"), lty = c(1,2,1,2), lwd = c(2,2,2,2));


jpeg(paste(wd,"plots_countries.jpg",sep = "/"), width=10, height=7.5, units="in", res=500)
table_in_italy <- (subset(mortality_all_latest, Country_code == "ITA"))
table_in_italy_M <- subset(table_in_italy, Gender == "MALE")
table_in_italy_F <- subset(table_in_italy, Gender == "FEMALE")
table_in_russia <- (subset(mortality_all_latest, Country_code == "RUS"))
table_in_russia_M <- subset(table_in_russia, Gender == "MALE")
table_in_russia_F <- subset(table_in_russia, Gender == "FEMALE")
par(mar = c(4.5,4.5,2,2))
layout(matrix(c(1,2),2,1, byrow = TRUE))
plot(table_in_italy_M$lx/1e5, lty = 1, lwd =2, type = "l", xlab = "Age", ylab = "Survival probability", cex.lab = 1.5, cex.axis = 1.2, main = "Survival Male")
lines(table_in_russia_M$lx/1e5, col = "red", lwd = 2)
lines(USA_table_male_2013$lx/1e5, col = "blue", lwd = 2)
legend("bottomleft", c("Italy 2009", "Russia 2010","USA 2013"), cex=1.2, col=c("black", "red", "blue"), lty = c(1,1,1), lwd = c(2,2,2));
par(mar = c(4.5,4.5,2,2) )
plot(table_in_italy_F$lx/1e5, lty = 1, lwd = 2,type = "l", xlab = "Age", ylab = "Survival probability", cex.lab = 1.5, cex.axis = 1.2, main = "Survival Female")
lines(table_in_russia_F$lx/1e5, col = "red", lwd = 2)
lines(USA_table_female_2013$lx/1e5, col = "blue", lwd =2)
dev.off()
