##########################################################################
# Download data from Human Mortality Database and prepare for analysis   #
##########################################################################

#############################################################################
# Obtain mortality data required for specific trial/analysis.
#  For the BRIM3 and coBRIM analyses, relevant countries are specified in
#  'country_list.csv'. The analysis is currently restricted to the year 2013.
#############################################################################

# Load custom functions ---------------------------------------------------
source("functions/funs_load_mort_table.R")


# Specify directory from which data are read ------------------------------
mainDir <- "data/mortality"


# Specify countries of interest -------------------------------------------
## Three column list: HMD code, country name, ISO3 code
country_list <- read.csv(paste0(mainDir, "/country_list.csv"))


# Bring mortality data in required form  ----------------------------------
for (i in 1:length(country_list$Code)) {
  subDir <- country_list$Country[i]

  if (!file.exists(paste(mainDir, subDir, sep = "/"))) {
    dir.create(file.path(mainDir, subDir))
  }
  
  # Check if the Male/Female subdirectories exist
  if (!file.exists(paste(mainDir, subDir, "Male", sep = "/"))) {
    dir.create(file.path(paste(mainDir, subDir, sep = "/"), "Male"))
  }
  if (!file.exists(paste(mainDir, subDir, "Female", sep = "/"))) {
    dir.create(file.path(paste(mainDir, subDir, sep = "/"), "Female"))
  }
  
  # Generate mortality tables
  table_gen(str_Country_Code = country_list$Code[i],
            str_Country = country_list$Country[i],
            sex = "M",
            parent_dir = "data")
  table_gen(str_Country_Code = country_list$Code[i],
            str_Country = country_list$Country[i],
            sex = "F",
            parent_dir = "data")
}


# Generate a file containing all mortality data ---------------------------

# Initialize data frame and choose year as 2013 ====
frame_all_countries <- data.frame()
year_in_use <- 2013

# Read, extract, combine, and write data ====
for (i in 1:length(country_list$Code)) {
  table_m <- read.csv(file = paste(mainDir, country_list$Country[i],
                                   "Male", "Mortality.csv", sep = "/"))
  table_f <- read.csv(file = paste(mainDir, country_list$Country[i],
                                   "Female", "Mortality.csv", sep = "/"))
  
  # Extract the specific year (maximum year )
  table_my <- subset(table_m, Year == max(table_m$Year))
  table_fy <- subset(table_f, Year == max(table_f$Year))
  table_my$Gender <- "MALE"
  table_fy$Gender <- "FEMALE"

  table_my$Country <- country_list$Country[i]
  table_fy$Country <- country_list$Country[i]
  table_my$Country_code <- country_list$Code3[i]
  table_fy$Country_code <- country_list$Code3[i]
  frame_all_countries <- rbind(frame_all_countries, table_fy)
  frame_all_countries <- rbind(frame_all_countries, table_my)
}

frame_all_countries[, "Age"] <- as.numeric(frame_all_countries[, "Age"])

write.csv(frame_all_countries,
          file = paste(mainDir, "mortality_all_latest.csv", sep = "/"))
