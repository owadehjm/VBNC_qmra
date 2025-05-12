# Computation of all the R codes
# This code runs the risk assessment models for scenarios both for the processed and bagged lettuce
# Get the R source code

source("processed.code.R")

# Processed lettuce
# Create a data frame for all scenarios
scenarios <- data.frame(
  Temp=c(0,9,9,9,9,9,9,17,17,17,17,17,17),
  Time=c(0, 0, 24, 48, 72, 96, 120,0, 24, 48, 72, 96, 120))

# Apply your function to each row and store results
processed_list <- lapply(1:nrow(scenarios), function(i) {
  result <- processed_scenario(scenarios$Temp[i], scenarios$Time[i])
  data.frame(Temp = scenarios$Temp[i],Time = scenarios$Time[i],
             result)
})

processed_df <- do.call(rbind, lapply(processed_list, as.data.frame))

processed_df$Time<-ifelse(processed_df$Temp==0 & processed_df$Time==0,
                          "Baseline",processed_df$Time)

write.csv(processed_df,"processed.dat.csv")

# Bagged lettuce
source("bagged.code.R")

# Apply your function to each row and store results
bagged_list <- lapply(1:nrow(scenarios), function(i) {
  result <- bagged_scenario(scenarios$Temp[i], scenarios$Time[i])
  # Combine with Temp and Time
  cbind(Temp = scenarios$Temp[i], Time = scenarios$Time[i], result)
})

bagged_df <- do.call(rbind, lapply(bagged_list, as.data.frame))

bagged_df$Time<-ifelse(bagged_df$Temp==0 & bagged_df$Time==0,
                       "Baseline",bagged_df$Time)

write.csv(bagged_df,"bagged.dat.csv")

