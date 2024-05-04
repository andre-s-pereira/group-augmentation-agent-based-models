agent_based_model <- function(n_groups, n_years) { # create function - then we run this function from another  R script
  
  final_dataset <- data.frame() # Create the final dataset, see below
  
  for (curgroup in 1:n_groups) {
    
    ##################################################
    #### Create the environment and the variables ####
    ##################################################
    
    group_id<- curgroup
    p_nonkin<- runif(1,0,1) # Probability that potential immigrants are non-kin
    year_harsh<- sort(sample(0:20,sample(1:21,1,replace=FALSE),replace=FALSE))
    year_opt<- sort(sample(0:20,sample(1:21,1,replace=FALSE),replace=FALSE))  
    harshness<- runif(1,0,1) #in the Thesis this is "how costly it is to diverge from optimal group size" - here i called it "harshness for simplicity";;;; needs to come back here for year 0, just in case year 0 is not sorted in the years_harsh
    optimal_gs<- 10 #for year 0 (no chance in optimal_gs) just in case there's no change in optimal_gs as per year_gs
    
    ###########################################
    #### create the original group members ####
    ########################################### 
    
    # this is our initial group
    # all individuals are related
    # individuals have random reproductive value
    
    group1 <- data.frame(
      age = sample(0:9, 10, replace=T), # 0 to 9 because the true age is +1 (line 34 of the code)
      gene = 0, #0= core of related group members or related to the core of related group members, 1= non-kin
      prob_acceptance = 0,
      prob_reproduction = 0,
      total_reproduction = 0,
      prob_survival = 0,
      survived = 0
    )
    
    for(curyear in 0:n_years){
      
      group1$age<- group1$age+1
      
      # no individual aged 10 should survive past the loop but this is just a "sanity check"
      group1<- subset(group1, age<11) 
      if(nrow(group1) == 0) {break}
      
      group1_start<-nrow(group1) #save group 1 as it was before group augmentation and deaths
      
      for(a in 1:length(year_harsh)){
        if(curyear==year_harsh[a]){harshness<-runif(1,0,1)}
          }
      for(a in 1:length(year_opt)){
        if(curyear==year_opt[a] & optimal_gs*2.5<1000){optimal_gs<-sample((optimal_gs+1):(optimal_gs*2.5),1)}
      }
      
      deviation<- ifelse(optimal_gs-nrow(group1)<=0, 0, optimal_gs-nrow(group1))
      immigrant_pool <- rpois(1,deviation)
      
      #####################
      #### Recruitment ####
      #####################
      
      n_recruits<-0
      
      if(nrow(group1)<2*optimal_gs){
        
      potential_immigrants<-data.frame()
      if (immigrant_pool>0){    
        potential_immigrants<- data.frame(
          age = sample(1:10,immigrant_pool,replace=T), # Recruits are already adults (otherwise they wouldn't migrate)
          gene = 0, #To prepare for the next step
          prob_acceptance = 0, 
          prob_reproduction = 0,
          total_reproduction = 0,
          prob_survival = 0,
          survived = 0,
          accepted = 0
        )
        for(d in 1:immigrant_pool){
          potential_immigrants$gene[d] <- rbinom(n=1,size=1,prob=p_nonkin)
          if(potential_immigrants$gene[d]>0){potential_immigrants$gene[d]<-sample(1:5,1)}
        }
      }
      
      if (immigrant_pool>0){    
        for (b in 1:immigrant_pool) { 
          n_votes <- 0
          for (c in 1:nrow(group1)){
            group1$prob_acceptance[c]<-runif(1,0,1)
            if(group1$prob_acceptance[c]>0.5){
              n_votes= n_votes+1
            }
          }
          if (n_votes >= (nrow(group1)/2)) {
            n_recruits = n_recruits+1
            potential_immigrants$accepted[b]<-1
          }
        }
      }
      
      }
      
      if(n_recruits>0){
        recruits<- subset(potential_immigrants, accepted==1, select = c(age, gene, prob_acceptance, prob_reproduction, total_reproduction, prob_survival, survived))
      } else {
        recruits<-data.frame()
      }
      
      ######################
      #### Reproduction ####
      ######################
      
      n_reproduced <- 0
      breeders<-c() # to keep track of what individuals reproduced
      if(curyear>0){group.size.fitness.rep=1-((final_group$distance_to_optimal_group_size/final_group$optimal_group_size)*(harshness_year_before))     
      } else {group.size.fitness.rep=1}
      
      for (e in 1:nrow(group1)) {
        
        #probability that each individual reproduces
        group1$prob_reproduction[e]= group.size.fitness.rep * (1-((group1$age[e]-1)/11)) 
        
        #did the individual reproduce?
        reproduced=0
        reproduced=rbinom(1,1,prob=group1$prob_reproduction[e])
        if (reproduced == 1) {
          n_reproduced = n_reproduced+1
          breeders[n_reproduced]=group1$gene[e]
          group1$total_reproduction[e]= group1$total_reproduction[e]+1
        }
      }
      
      if(length(breeders)>0){
        offspring <- data.frame(
          age = 0,
          gene = breeders,
          prob_acceptance = 0,
          prob_reproduction = 0, 
          total_reproduction = 0,
          prob_survival = 0,
          survived = 0
        )
      } else {
        offspring<-data.frame()
      }
      
      n.offspring.kin<-sum(breeders==0)
      n.offspring.nonkin<-sum(breeders>0)
      pro_rep<-mean(group1$prob_reproduction)
      
      ###################################
      #### Now that all individuals who were not recruited this round have reproduced, we add the recruits and offspring (easy fix around the recruits not voting/ reproducing)
      ###################################
      
      # First we find the difference between the additions and the maximum group size
      
      if(nrow(group1)<2*optimal_gs){
        
        n_additions=n_recruits+n_reproduced
        dev=(2*optimal_gs-nrow(group1))-n_additions
        
        if(dev<0){
          offs<-rep(0,nrow(offspring))
          imm<-rep(1,nrow(recruits))
          added<-c(offs,imm)
          final<-sample(added,(n_additions+dev)) # it's "+" because dev is negative, so this is what we have to calculate to know how many recruits need to be added to reach maximal group size
          offspring<-offspring[sample(nrow(offspring), sum(final==0)), ]
          recruits<-recruits[sample(nrow(recruits), sum(final==1)), ]
        }
        
        # Now we add the immigrants and offspring that decided not to (re)emmigrate 
        additions<-rbind(recruits,offspring)
        group1<-rbind(group1,additions)
        
      }
      
      #################################
      #### Disappearance of individuals
      #################################
      
      pre_disappearances<-nrow(group1)
      distance_to_optimal_group_size_neg <- nrow(group1)-optimal_gs
      distance_to_optimal_group_size <-abs(nrow(group1)-optimal_gs)
      group.size.fitness.surv<-  1-((distance_to_optimal_group_size/optimal_gs)*(harshness))
      
      for (f in 1:nrow(group1)) {
        #probability that each individual survives
        group1$prob_survival[f]= group.size.fitness.surv * (1-(group1$age[f]/10))
        group1$survived[f] = rbinom(1,1,prob=group1$prob_survival[f])
      }
      
      pro_survival<-mean(group1$prob_survival)
      n_disappearances<-nrow(subset(group1, survived==0))
      group1<-subset(group1, survived==1)
      harshness_year_before<-harshness
      
      ####################################
      #### Final dataset of the group ####
      ####################################
      
      if (nrow(group1)==0){
          final_group <- data.frame(  
          group = group_id,
          year = curyear,
          type_of_model = "recruit kin and non-kin",
          optimal_group_size = optimal_gs,
          deviation = deviation,
          harshness = harshness,
          n_recruits = n_recruits,
          n_offspring_kin = n.offspring.kin,
          n_offspring_nonkin = n.offspring.nonkin,
          group_size_fitness.rep = group.size.fitness.rep,
          group_size_fitness.surv = group.size.fitness.surv,
          prob_rep = pro_rep,
          n_reproduced = n_reproduced,
          p_reproduced = n_reproduced/group1_start,
          final_group_size_pre_deaths = pre_disappearances,
          final_group_size_post_deaths = nrow(group1),
          n_kin = 0,
          n_nonkin = 0,
          p_nonkin = p_nonkin,
          distance_to_optimal_group_size_neg= distance_to_optimal_group_size_neg,
          distance_to_optimal_group_size= distance_to_optimal_group_size,
          proportion_of_kin_after_deaths= 0,
          prob_survival = pro_survival,
          n_disappearances=n_disappearances,
          mix_related=NA
        ) 
      } else {
          final_group <- data.frame(   
          group = group_id,
          year = curyear,
          type_of_model = "recruit kin and non-kin",
          optimal_group_size = optimal_gs,
          deviation = deviation,
          harshness=harshness,
          n_recruits = n_recruits,
          n_offspring_kin = n.offspring.kin,
          n_offspring_nonkin = n.offspring.nonkin,
          group_size_fitness.rep = group.size.fitness.rep,
          group_size_fitness.surv = group.size.fitness.surv,
          prob_rep = pro_rep,
          n_reproduced = n_reproduced,
          p_reproduced = n_reproduced/group1_start,
          final_group_size_pre_deaths = pre_disappearances,
          final_group_size_post_deaths = nrow(group1),
          n_kin = sum(group1$gene==0),
          n_nonkin = sum(group1$gene>0),
          p_nonkin = p_nonkin,
          distance_to_optimal_group_size_neg= distance_to_optimal_group_size_neg,
          distance_to_optimal_group_size= distance_to_optimal_group_size,
          proportion_of_kin_after_deaths= sum(group1$gene==0)/nrow(group1),
          n_disappearances=n_disappearances,
          prob_survival = pro_survival,
          mix_related=ifelse(sum(group1$gene>0)>0,1,0)
          )
      }
      
      final_dataset <- rbind(final_dataset, final_group) 
      
      if(nrow(group1) == 0) {break}
      
      table<-c(curgroup,curyear)
      print(table)
    }
    
    save(final_dataset, file="agent_based_model_dataset.Rdata")
  }
  
}

