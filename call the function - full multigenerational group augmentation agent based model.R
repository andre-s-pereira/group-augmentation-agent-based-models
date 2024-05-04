setwd("C:/Users/... ... ...") # where you stored the function for the agent based model

library(tictoc) # just to keep track of how fast the model is running
source("agent_based_model.r")

tic()
model<-agent_based_model(5000,20) #n_groups=5000, n_years=20
toc()