#This R script is conducted using toydata on the Spike gene as example.

#Install packages and set library package
install.packages("prettyR")

library(prettyR)
library(scales)

#Set working dictionary
setwd("./S-EPS/data")

#The analyzed regions include Australia, Brazil, California, Germany, India, Israel, Japan, Mexico, New York, Russia, Singapore, South Africa, and United Kingdom
region <- c("AUS", "BRA", "CA", "GER", "IND",
            "ISR", "JPN", "MEX", "NY", "RUS",
            "SA", "SGP", "UK")

#Step 1: Identify Delta-characteristic key mutations
#Use three substeps to filter out the substitutions that meet the criteria for Uniqueness, Substitutes for Fitness and Coverage.
#Step 1.1: Fitness criteria.
#Identify the mutations that ever reached dominance in at least one of the selected regions during the period (the Fitness criteria)
#Run time: ~1 min

fit_sub <- function(selected_region) {
  list_fit_sub <- vector()
  sample_size <- read.csv(file = "input/sample_size.csv", 
                          header = TRUE,
                          row.names = 1)
  month <- rownames(sample_size)
  reference_aa_seq <- read.csv(file = "input/wuhan_reference_spike_aa.csv",
                               header = TRUE,
                               row.names = 1,
                               colClasses = "character")
  for (l in 1:length(selected_region)) {
    for (t in 1:length(month)) {
      if (file.exists(paste0("input/first stage sampling/", 
                             selected_region[l], "/", month[t], ".csv"))) {
        monthly_aa <- read.csv(file = paste0("input/first stage sampling/", 
                                             selected_region[l], "/", month[t], ".csv"),
                               header = TRUE,
                               row.names = 1,
                               colClasses = "character")
        for (j in 1:ncol(monthly_aa)) {
          if (length(which(monthly_aa[,j] != reference_aa_seq[1,j] &
                           monthly_aa[,j] != "---")) >= nrow(monthly_aa)/2) {
            list_fit_sub <- append(list_fit_sub,
                                       paste0(reference_aa_seq[1,j],
                                              j,
                                              Mode(monthly_aa[,j])))
          }
        }
      }
    }
    list_fit_sub <- unique(list_fit_sub)
  }
  return(list_fit_sub)
}

sub_fitness <- fit_sub(selected_region = region)
#sub_fitness: list of substitutions that meet the Fitness criteria of key mutations as described in the Methods section

#write.table(sub_fitness, file = "output/sub_fitness.csv", row.names = FALSE, col.names = "sub_fitness")

#Step 1.2 Uniqueness criteria
#Identify the mutations that categorized to Delta-characteristic mutations (the Uniqueness criteria) 

uniq_sub <- function(list_sub) {
  delta_sub <- read.table(file = "input/delta mutations.csv",
                          colClasses = "character")
  list_unique_sub <- intersect(as.matrix(delta_sub), list_sub)
  return(list_unique_sub)
}

sub_fitness_uniqueness <- uniq_sub(list_sub = sub_fitness)
#sub_fitness_uniqueness: list of substitutions that meet the Uniqueness and Fitnes criteria of key mutations as described in the Methods section
#write.table(sub_fitness_uniqueness, file = "output/sub_fitness_uniqueness.csv", row.names = FALSE, col.names = "sub_fitness_uniqueness")

#Step 1.3 Coverage criteria
#Identify the mutations that have been observed with a prevalence exceeding 20% in more than half of the regions (the Coverage criteria)
#Run time: ~1 min

cov_sub <- function(list_sub, selected_region) {
  sample_size <- read.csv(file = "input/sample_size.csv", 
                          header = TRUE,
                          row.names = 1)
  month <- rownames(sample_size)
  sub_prev <- as.data.frame(matrix(nrow = length(month),
                                   ncol = length(list_sub)*length(selected_region)))
  
  for (l in 1:length(selected_region)) {
    for (t in 1:length(month)) {
      if (file.exists(paste0("input/first stage sampling/", 
                             selected_region[l], "/", month[t], ".csv"))) {
        monthly_aa <- read.csv(file = paste0("input/first stage sampling/", 
                                             selected_region[l], "/", month[t], ".csv"),
                               header = TRUE,
                               row.names = 1,
                               colClasses = "character")
        for (k in 1:length(list_sub)) {
          site <- as.numeric(substr(list_sub[k],
                                    2,
                                    nchar(list_sub[k])-1))
          sub_aa <- as.character(substr(list_sub[k],
                                        nchar(list_sub[k]),
                                        nchar(list_sub[k])))
          sub_prev[t,length(selected_region)*k- length(selected_region) + l] <- length(which(monthly_aa[,site] == sub_aa)) /
            nrow(monthly_aa)
        }
      }
    }
  }
  non_cov_sub <- vector()
  for (k in 1:length(list_sub)) {
    max_prev <- vector()
    for (l in 1:length(selected_region)) {
      max_prev[l] <- max(sub_prev[,length(selected_region)*k-length(selected_region)+l], na.rm = TRUE)
    }
    if (length(which(max_prev >= 0.2)) >= length(selected_region)/2) {
      non_cov_sub <- append(non_cov_sub, selected_region[k])
    }
  }
  list_cov_sub <- setdiff(list_sub, non_cov_sub)
  return(list_cov_sub)
}

key_sub <- cov_sub(sub_fitness_uniqueness, selected_region = region)

#key_sub: list of key substitutions that meet the Uniqueness, Fitness, and Coverage criteria as described in the Methods section
#write.table(key_sub, file = "output/key_sub.csv", row.names = FALSE, col.names = "key_sub")




#Step 2: Determine occurrence ranking using Equal Power Sampling (EPS) Framework
#Run time: ~1 min

sub_detection_rank <- function(selected_region, list_sub) {
  sample_size <- read.csv(file = "input/sample_size.csv", 
                          header = TRUE,
                          row.names = 1)
  month <- rownames(sample_size)
  sub_prev1 <- as.data.frame(matrix(nrow = length(month),
                                    ncol = length(list_sub)*length(selected_region)))
  
  for (l in 1:length(selected_region)) {
    for (t in 1:length(month)) {
      if (file.exists(paste0("input/first stage sampling/", 
                             selected_region[l], "/", month[t], ".csv"))) {
        monthly_aa <- read.csv(file = paste0("input/first stage sampling/", 
                                             selected_region[l], "/", month[t], ".csv"),
                               header = TRUE,
                               row.names = 1,
                               colClasses = "character")
        for (k in 1:length(list_sub)) {
          site <- as.numeric(substr(list_sub[k],
                                    2,
                                    nchar(list_sub[k])-1))
          sub_aa <- as.character(substr(list_sub[k],
                                        nchar(list_sub[k]),
                                        nchar(list_sub[k])))
          sub_prev1[t,length(selected_region)*k- length(selected_region) + l] <- length(which(monthly_aa[,site] == sub_aa)) /
            nrow(monthly_aa)
        }
      }
    }
  }
  sub_prev2 <- as.data.frame(matrix(nrow = length(month),
                                    ncol = length(list_sub)*length(selected_region)))
  
  for (l in 1:length(selected_region)) {
    for (t in 1:length(month)) {
      if (file.exists(paste0("input/first stage sampling/", 
                             selected_region[l], "/", month[t], ".csv"))) {
        monthly_aa <- read.csv(file = paste0("input/second stage sampling/", 
                                             selected_region[l], "/", month[t], ".csv"),
                               header = TRUE,
                               row.names = 1,
                               colClasses = "character")
        for (k in 1:length(list_sub)) {
          site <- as.numeric(substr(list_sub[k],
                                    2,
                                    nchar(list_sub[k])-1))
          sub_aa <- as.character(substr(list_sub[k],
                                        nchar(list_sub[k]),
                                        nchar(list_sub[k])))
          sub_prev2[t,length(selected_region)*k- length(selected_region) + l] <- length(which(monthly_aa[,site] == sub_aa)) /
            nrow(monthly_aa)
        }
      }
    }
  }
  list_sub_summary <- as.data.frame(matrix(nrow = length(list_sub),
                                           ncol = length(selected_region)*4))
  rownames(list_sub_summary) <- list_sub
  colnames(list_sub_summary) <- c(paste0(selected_region, "_detection_t1"),
                                  paste0(selected_region, "_detection_p1"),
                                  paste0(selected_region, "_detection_t2"),
                                  paste0(selected_region, "_detection_rank"))
  for (k in 1:length(list_sub)) {
    for (l in 1:length(selected_region)) {
      monthly_prev1 <- sub_prev1[,length(selected_region)*k - length(selected_region) + l]
      monthly_prev2 <- sub_prev2[,length(selected_region)*k - length(selected_region) + l]
      if (max(monthly_prev1, na.rm = TRUE) >= 0.2) {
        prevailed_t1 <- min(which(monthly_prev1 >= 0.2))
        list_sub_summary[k,l] <- max(which(monthly_prev1[1:prevailed_t1] == 0)) + 1
        list_sub_summary[k,l + length(selected_region)] <- monthly_prev1[list_sub_summary[k,l]]
      }
      if (max(monthly_prev2, na.rm = TRUE) >= 0.2) {
        prevailed_t2 <- min(which(monthly_prev2 >= 0.2))
        list_sub_summary[k,l + length(selected_region)*2] <- max(which(monthly_prev2[1:prevailed_t2] == 0)) + 1
      }
    }
    #Step 2.1: Identification of source region
    #Determine the source region using the first stage sampling data
    detection_rank1 <- rank(list_sub_summary[k,1:(length(selected_region))],
                            na.last = "keep")
    #There may be two or more regions that discover novel key mutations earliest in the same month,
    #and we recognize the one with the highest prevalence in that month as the source
    #pot_source_k: potential source of mutation k (region(s) with the earliest detection of key mutation k)
    pot_source_k <- which(detection_rank1 == min(detection_rank1, na.rm = TRUE))
    source_k <- pot_source_k[which.max(detection_rank1[pot_source_k])]
    list_sub_summary[k,length(selected_region)*3 + source_k] <- 1
    
    #Step 2.2: Occurrence ranking of the remaining regions
    #Regions with insufficient sample size to detect mutation occurrence are not included in the ranking
    #Determine the occurrence rank using the second stage sampling data
    detection_rank2 <- rank(list_sub_summary[k,(length(selected_region)*2 + 1):(length(selected_region)*3)], 
                            na.last = "keep")
    for (l in setdiff(1:length(selected_region), source_k)) {
      if ((!is.na(detection_rank2[l]))&
          (detection_rank2[l] != min(detection_rank2,na.rm = TRUE)) &
          (sample_size[list_sub_summary[k,l + length(selected_region)*2]-1,l] < 163)) {
        detection_rank2[l] <- NA
      }
    }
    detection_rank2[source_k] <- 1
    
    #Step 2.3: Normalization
    #Normalize the ranking to a scale of 1-13 to facilitate joint analysis of all mutations
    orig_ranking_max <- max(rank(detection_rank2[-source_k], na.last = "keep"), na.rm = TRUE)
    orig_ranking_min <- min(rank(detection_rank2[-source_k], na.last = "keep"), na.rm = TRUE)
    detection_rank2[-source_k] <- round(((length(selected_region)-2)/(orig_ranking_max-orig_ranking_min))*rank(detection_rank2[-source_k], na.last = "keep") + 
                                          (2*orig_ranking_max-length(selected_region)*orig_ranking_min)/(orig_ranking_max-orig_ranking_min))
    list_sub_summary[k,(length(selected_region)*3 + 1):(length(selected_region)*4)] <- detection_rank2
  }
  sub_rank <- list_sub_summary[,(length(selected_region)*3 + 1):(length(selected_region)*4)]
  return(sub_rank)
}

key_sub_rank <- sub_detection_rank(selected_region = region, list_sub = key_sub)

#key_sub_rank: Detection rank of each region along the transmission pathway of each key mutation
#write.csv(key_sub_rank, file = "output/key_sub_rank.csv")




#Step 3: Calculate rank statistics probability
rank_prob <- function(sub_rank) {
  rank_probability <- as.data.frame(matrix(nrow = length(region),
                                           ncol = length(region)))
  colnames(rank_probability) <- region
  rownames(rank_probability) <- paste("Probability of rank =",1:length(region))
  for (l in 1:length(region)) {
    for (r in 1:length(region)) {
      rank_probability[r,l] <- percent(length(which(sub_rank[,l] == r)) / nrow(sub_rank), accuracy = 0.1)
    }
  }
  cols_order <- order(as.matrix(rank_probability[1,]), 
                      as.matrix(rank_probability[2,]), 
                      decreasing = TRUE)
  rank_probability_sorted <- rank_probability[, cols_order]
  return(rank_probability_sorted)
}
prob_key_sub_rank <- rank_prob(key_sub_rank)

#prob_key_sub_rank: matrix of probability of each region being the rth order along the key mutation transmission among the selected regions
#write.csv(prob_key_sub_rank, file = "output/prob_key_sub_rank.csv")