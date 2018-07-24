#####################################
# Clean and format mortality tables #
#####################################

##########################################################################
# Functions to clean and format the mortality tables obtained from the HMD
#  for use in the analysis
##########################################################################

# Read in list of countries
country_List <- read.csv("libraries/CountryList.csv", sep = ",")
mainDir <- paste0(getwd(), "/libraries", sep = "")

source("mortality_table_loader.R")

# Create mortality files
for (i in 1:length(country_List$Code)) {
  subDir <- country_List$Country[i]

  if (!file.exists(paste(mainDir, subDir, sep = "/"))) {
    dir.create(file.path(mainDir, subDir))
  }
  # check if the Female/Male subdirs exist
  if (!file.exists(paste(mainDir, subDir, "Female", sep = "/"))) {
    dir.create(file.path(paste(mainDir, subDir, sep = "/"), "Female"))
  }
  if (!file.exists(paste(mainDir, subDir, "Male", sep = "/"))) {
    dir.create(file.path(paste(mainDir, subDir, sep = "/"), "Male"))
  }
  table_gen(str_Country_Code = country_List$Code[i],
            str_Country = country_List$Country[i], sex = "F",
            parent_dir = "libraries")
  table_gen(str_Country_Code = country_List$Code[i],
            str_Country = country_List$Country[i], sex = "M",
            parent_dir = "libraries")
}

# Generate a common file
## Currently, only the year 2013 is considered but future updates will include
##  additional years
frame_all_countries <- data.frame()
year_in_use <- 2013

for (i in 1:length(country_List$Code)) {
  table_m <- read.csv(file = paste(mainDir, country_List$Country[i],
                                   "Male", "Mortality.csv", sep = "/"))
  table_f <- read.csv(file = paste(mainDir, country_List$Country[i],
                                   "Female", "Mortality.csv", sep = "/"))
  
  # Extract the specific year (maximum year )
  table_my <- subset(table_m, Year == max(table_m$Year))
  table_fy <- subset(table_f, Year == max(table_f$Year))
  
  table_my$Gender <- "MALE"
  table_fy$Gender <- "FEMALE"

  table_my$Country <- country_List$Country[i]
  table_fy$Country <- country_List$Country[i]
  
  table_my$Country_code <- country_List$Code3[i]
  table_fy$Country_code <- country_List$Code3[i]
  
  frame_all_countries <- rbind(frame_all_countries, table_fy)
  frame_all_countries <- rbind(frame_all_countries, table_my)
}

frame_all_countries[, "Age"] <- as.numeric(frame_all_countries[, "Age"])

# Write data to csv file
write.csv(frame_all_countries,
          file = paste(mainDir, "Mortality_all_latest.csv", sep = "/"),
          sep = ",")
