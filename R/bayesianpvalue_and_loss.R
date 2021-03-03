out_dir <- "chessB_pn/"

## posterior predictive checking - tempest model
load(paste0(out_dir, "surv_sim.RData")) # load surv_sim.RData (same data as csv files in postpredcheck folder)
library(magrittr)

mt = readr::read_csv(file="chessB_pn/input/mt.csv", col_names=F) %>% as.matrix() # observed response time data (40 x 251)
mi = readr::read_csv(file="chessB_pn/input/mi.csv", col_names=F) %>% as.matrix() # observed response data (40 x 251)

## time = sim_data[[item]][,1:251]
## pp = sim_data[[item]][,252:502]

ppp_matrix_rt <- matrix(0, nrow = 40, ncol = 251) # setting an empty matrix for bayesianpvalue of each cell
for (i in 1:nrow(ppp_matrix_rt)) {
  time <- sim_data[[i]][, 1:251] # simulated response time data of an item i for every examinee (3000 x 251)
  for (j in 1:ncol(ppp_matrix_rt)) {
    ppp_matrix_rt[i, j] <- sum(mt[i, j] < time[, j]) / 3000 # calculate bayesianpvalue for each item&examinee pair
  }
}

## bayesianpvalue for the correct responses
ppp_matrix_rt_correct <- matrix(0, nrow = 40, ncol = 251) # setting an empty matrix for bayesianpvalue of each cell
for (i in 1:nrow(ppp_matrix_rt_correct)) {
  for (j in 1:ncol(ppp_matrix_rt_correct)) {
    if (mi[i, j] == 1) { # if an examinee got an item correct,
      ppp_matrix_rt_correct[i, j] <- ppp_matrix_rt[i, j] # extract bayesianpvalue for the correct responses
    } else {
      ppp_matrix_rt_correct[i, j] <- NA
    }
  }
}

## incorrect
ppp_matrix_rt_incorrect <- matrix(0, nrow = 40, ncol = 251) # setting an empty matrix for bayesianpvalue of each cell
for (i in 1:nrow(ppp_matrix_rt_incorrect)) {
  for (j in 1:ncol(ppp_matrix_rt_incorrect)) {
    if (mi[i, j] == 0) { # if an examinee got an item incorrect,
      ppp_matrix_rt_incorrect[i, j] <- ppp_matrix_rt[i, j] # extract bayesianpvalue for the incorrect responses
    } else {
      ppp_matrix_rt_incorrect[i, j] <- NA
    }
  }
}

sum(is.na(ppp_matrix_rt_correct)) # 4950 -> number of correct responses: 5090, number of incorrect responses: 4950

## 0.05 < bayesianpvalue <0.95
sum((ppp_matrix_rt > 0.05) & (ppp_matrix_rt < 0.95)) / (40 * 251) # 0.8705179
sum((ppp_matrix_rt_correct > 0.05) & (ppp_matrix_rt_correct < 0.95), na.rm = TRUE) / 5090 # 0.802947
sum((ppp_matrix_rt_incorrect > 0.05) & (ppp_matrix_rt_incorrect < 0.95), na.rm = TRUE) / 4950 # 0.94
## 0.1 < proportion <0.9
sum((ppp_matrix_rt > 0.1) & (ppp_matrix_rt < 0.9)) / (40 * 251) # 0.7534861
sum((ppp_matrix_rt_correct > 0.1) & (ppp_matrix_rt_correct < 0.9), na.rm = TRUE) / 5090 # 0.6638507
sum((ppp_matrix_rt_incorrect > 0.1) & (ppp_matrix_rt_incorrect < 0.9), na.rm = TRUE) / 4950 # 0.8456566
## 0.2 < proportion <0.8
sum((ppp_matrix_rt > 0.2) & (ppp_matrix_rt < 0.8)) / (40 * 251) # 0.5660359
sum((ppp_matrix_rt_correct > 0.2) & (ppp_matrix_rt_correct < 0.8), na.rm = TRUE) / 5090 # 0.4960707
sum((ppp_matrix_rt_incorrect > 0.2) & (ppp_matrix_rt_incorrect < 0.8), na.rm = TRUE) / 4950 # 0.6379798

## bayesianpvalue histograms
hist(ppp_matrix_rt, nclass = 100)
abline(v = 0.5)

par(mfrow = c(1, 2))
hist(ppp_matrix_rt_correct, nclass = 100, ylim = c(0, 250))
abline(v = 0.5)
hist(ppp_matrix_rt_incorrect, nclass = 100, ylim = c(0, 250))
abline(v = 0.5)

## Loss function(L1)_response time
l1_matrix_rt <- matrix(0, nrow = 40, ncol = 251) # setting an empty matrix for bayesianpvalue of each cell
for (i in 1:40) {
  time <- sim_data[[i]][, 1:251] # simulated response time data of an item i for every examinee (3000 x 251)
  for (j in 1:251) {
    l1_matrix_rt[i, j] <- sum(abs(mt[i, j] - time[, j])) # calculating loss1 value for each item&examinee pair
  }
}

## loss1 for the correct responses
l1_matrix_rt_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_rt_correct)) {
  for (j in 1:ncol(l1_matrix_rt_correct)) {
    if (mi[i, j] == 1) {
      l1_matrix_rt_correct[i, j] <- l1_matrix_rt[i, j]
    } else {
      l1_matrix_rt_correct[i, j] <- NA
    }
  }
}

## loss1 for the incorrect responses
l1_matrix_rt_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_rt_incorrect)) {
  for (j in 1:ncol(l1_matrix_rt_incorrect)) {
    if (mi[i, j] == 0) {
      l1_matrix_rt_incorrect[i, j] <- l1_matrix_rt[i, j]
    } else {
      l1_matrix_rt_incorrect[i, j] <- NA
    }
  }
}

## Loss function(L2)_response time
l2_matrix_rt <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:40) {
  time <- sim_data[[i]][, 1:251]
  for (j in 1:251) {
    l2_matrix_rt[i, j] <- sum((mt[i, j] - time[, j])^2)
  }
}

## Loss2 for the correct responses
l2_matrix_rt_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_rt_correct)) {
  for (j in 1:ncol(l2_matrix_rt_correct)) {
    if (mi[i, j] == 1) {
      l2_matrix_rt_correct[i, j] <- l2_matrix_rt[i, j]
    } else {
      l2_matrix_rt_correct[i, j] <- NA
    }
  }
}

## loss2 for the incorrect responses
l2_matrix_rt_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_rt_incorrect)) {
  for (j in 1:ncol(l2_matrix_rt_incorrect)) {
    if (mi[i, j] == 0) {
      l2_matrix_rt_incorrect[i, j] <- l2_matrix_rt[i, j]
    } else {
      l2_matrix_rt_incorrect[i, j] <- NA
    }
  }
}

## histograms of loss values
par(mfrow = c(1, 2))
hist(l1_matrix_rt, nclass = 100)
hist(l2_matrix_rt, nclass = 100)

mean(l1_matrix_rt)
mean(l2_matrix_rt)

par(mfrow = c(2, 2))
hist(l1_matrix_rt_correct, nclass = 100)
hist(l1_matrix_rt_incorrect, nclass = 100)
hist(l2_matrix_rt_correct, nclass = 100)
hist(l2_matrix_rt_incorrect, nclass = 100)

mean(l1_matrix_rt_correct, na.rm = TRUE) # 31180.79
mean(l1_matrix_rt_incorrect, na.rm = TRUE) # 28380.64
mean(l2_matrix_rt_correct, na.rm = TRUE) # 568381.4
mean(l2_matrix_rt_incorrect, na.rm = TRUE) # 424044.6: incorrect에서의 loss가 더 작다.


## loss function(L1)_probability of getting an answer correct
l1_matrix <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:40) {
  prob <- sim_data[[i]][, 252:502] # simulated probability of getting an item i correct for every examinee (3000 x 251)
  for (j in 1:251) {
    l1_matrix[i, j] <- sum(abs(mi[i, j] - prob[, j]))
  }
}

## loss1 for the correct responses
l1_matrix_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_correct)) {
  for (j in 1:ncol(l1_matrix_correct)) {
    if (mi[i, j] == 1) {
      l1_matrix_correct[i, j] <- l1_matrix[i, j]
    } else {
      l1_matrix_correct[i, j] <- NA
    }
  }
}

## loss1 for the incorrect responses
l1_matrix_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_incorrect)) {
  for (j in 1:ncol(l1_matrix_incorrect)) {
    if (mi[i, j] == 0) {
      l1_matrix_incorrect[i, j] <- l1_matrix[i, j]
    } else {
      l1_matrix_incorrect[i, j] <- NA
    }
  }
}

## loss function(L2)_probability of getting an answer correct
l2_matrix <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:40) {
  prob <- sim_data[[i]][, 252:502]
  for (j in 1:251) {
    l2_matrix[i, j] <- sum((mi[i, j] - prob[, j])^2)
  }
}

l2_matrix_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_correct)) {
  for (j in 1:ncol(l2_matrix_correct)) {
    if (mi[i, j] == 1) {
      l2_matrix_correct[i, j] <- l2_matrix[i, j]
    } else {
      l2_matrix_correct[i, j] <- NA
    }
  }
}

l2_matrix_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_incorrect)) {
  for (j in 1:ncol(l2_matrix_incorrect)) {
    if (mi[i, j] == 0) {
      l2_matrix_incorrect[i, j] <- l2_matrix[i, j]
    } else {
      l2_matrix_incorrect[i, j] <- NA
    }
  }
}

par(mfrow = c(1, 2))
hist(l1_matrix, nclass = 100, ylim = c(0, 500))
hist(l2_matrix, nclass = 100, ylim = c(0, 500))

par(mfrow = c(2, 2))
hist(l1_matrix_correct, nclass = 100)
hist(l1_matrix_incorrect, nclass = 100)
hist(l2_matrix_correct, nclass = 100)
hist(l2_matrix_incorrect, nclass = 100)

mean(l1_matrix_correct, na.rm = TRUE) # 1610.3
mean(l1_matrix_incorrect, na.rm = TRUE) # 1356.137
mean(l2_matrix_correct, na.rm = TRUE) # 1177.837
mean(l2_matrix_incorrect, na.rm = TRUE) # 909.3655


### ---------------- THIS IS THE END OF THE CODES FOR THE TEMPEST MODEL ----------------- ##
### ---------------- Same codes were applied to the no latent model ------------------- ##


out_dir <- "chessB_no_latent/"

## posterior predictive checking - tempest model
load(paste0(out_dir, "surv_sim.RData")) # load surv_sim.RData (same data as csv files in postpredcheck folder)
library(magrittr)

mt = readr::read_csv(file="chessB_pn/input/mt.csv", col_names=F) %>% as.matrix() # observed response time data (40 x 251)
mi = readr::read_csv(file="chessB_pn/input/mi.csv", col_names=F) %>% as.matrix() # observed response data (40 x 251)
## time = sim_data[[item]][,1:251]
## pp = sim_data[[item]][,252:502]

ppp_matrix_rt_no_latent <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(ppp_matrix_rt_no_latent)) {
  time <- sim_data[[i]][, 1:251]
  for (j in 1:ncol(ppp_matrix_rt_no_latent)) {
    ppp_matrix_rt_no_latent[i, j] <- sum(mt[i, j] < time[, j]) / 3000
  }
}

## correct
ppp_matrix_rt_no_latent_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(ppp_matrix_rt_no_latent_correct)) {
  for (j in 1:ncol(ppp_matrix_rt_no_latent)) {
    if (mi[i, j] == 1) {
      ppp_matrix_rt_no_latent_correct[i, j] <- ppp_matrix_rt_no_latent[i, j]
    } else {
      ppp_matrix_rt_no_latent_correct[i, j] <- NA
    }
  }
}

## incorrect
ppp_matrix_rt_no_latent_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(ppp_matrix_rt_no_latent_incorrect)) {
  for (j in 1:ncol(ppp_matrix_rt_no_latent)) {
    if (mi[i, j] == 0) {
      ppp_matrix_rt_no_latent_incorrect[i, j] <- ppp_matrix_rt_no_latent[i, j]
    } else {
      ppp_matrix_rt_no_latent_incorrect[i, j] <- NA
    }
  }
}

## 0.05 < proportion <0.95 # correct/incorrect 구분해보기
sum((ppp_matrix_rt_no_latent > 0.05) & (ppp_matrix_rt_no_latent < 0.95)) / (40 * 251) # 0.947012
sum((ppp_matrix_rt_no_latent_correct > 0.05) & (ppp_matrix_rt_no_latent_correct < 0.95), na.rm = TRUE) / 5090 # 0.9495088
sum((ppp_matrix_rt_no_latent_incorrect > 0.05) & (ppp_matrix_rt_no_latent_incorrect < 0.95), na.rm = TRUE) / 4950 # 0.9444444
## 0.1 < proportion <0.9
sum((ppp_matrix_rt_no_latent > 0.1) & (ppp_matrix_rt_no_latent < 0.9)) / (40 * 251) # 0.8318725
sum((ppp_matrix_rt_no_latent_correct > 0.1) & (ppp_matrix_rt_no_latent_correct < 0.9), na.rm = TRUE) / 5090 # 0.8227898
sum((ppp_matrix_rt_no_latent_incorrect > 0.1) & (ppp_matrix_rt_no_latent_incorrect < 0.9), na.rm = TRUE) / 4950 # 0.8412121

## 0.2 < proportion <0.8
sum((ppp_matrix_rt_no_latent > 0.2) & (ppp_matrix_rt_no_latent < 0.8)) / (40 * 251) # 0.5477092
sum((ppp_matrix_rt_no_latent_correct > 0.2) & (ppp_matrix_rt_no_latent_correct < 0.8), na.rm = TRUE) / 5090 # 0.5485265
sum((ppp_matrix_rt_no_latent_incorrect > 0.2) & (ppp_matrix_rt_no_latent_incorrect < 0.8), na.rm = TRUE) / 4950 # 0.5468687

hist(ppp_matrix_rt_no_latent, nclass = 100, ylim = c(0, 250))
abline(v = 0.5)
par(mfrow = c(1, 2))
hist(ppp_matrix_rt_no_latent_correct, nclass = 100, ylim = c(0, 250))
abline(v = 0.5)
hist(ppp_matrix_rt_no_latent_incorrect, nclass = 100, ylim = c(0, 250))
abline(v = 0.5)


## loss function(L1, L2)_rt correct/incorrect 구분해보기
l1_matrix_rt_no_latent <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:40) {
  time <- sim_data[[i]][, 1:251]
  for (j in 1:251) {
    l1_matrix_rt_no_latent[i, j] <- sum(abs(mt[i, j] - time[, j]))
  }
}

l1_matrix_rt_no_latent_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_rt_no_latent_correct)) {
  for (j in 1:ncol(l1_matrix_rt_no_latent_correct)) {
    if (mi[i, j] == 1) {
      l1_matrix_rt_no_latent_correct[i, j] <- l1_matrix_rt_no_latent[i, j]
    } else {
      l1_matrix_rt_no_latent_correct[i, j] <- NA
    }
  }
}

l1_matrix_rt_no_latent_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_rt_no_latent_incorrect)) {
  for (j in 1:ncol(l1_matrix_rt_no_latent_incorrect)) {
    if (mi[i, j] == 0) {
      l1_matrix_rt_no_latent_incorrect[i, j] <- l1_matrix_rt_no_latent[i, j]
    } else {
      l1_matrix_rt_no_latent_incorrect[i, j] <- NA
    }
  }
}

l2_matrix_rt_no_latent <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:40) {
  time <- sim_data[[i]][, 1:251]
  for (j in 1:251) {
    l2_matrix_rt_no_latent[i, j] <- sum((mt[i, j] - time[, j])^2)
  }
}


l2_matrix_rt_no_latent_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_rt_no_latent_correct)) {
  for (j in 1:ncol(l2_matrix_rt_no_latent_correct)) {
    if (mi[i, j] == 1) {
      l2_matrix_rt_no_latent_correct[i, j] <- l2_matrix_rt_no_latent[i, j]
    } else {
      l2_matrix_rt_no_latent_correct[i, j] <- NA
    }
  }
}

l2_matrix_rt_no_latent_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_rt_no_latent_incorrect)) {
  for (j in 1:ncol(l2_matrix_rt_no_latent_incorrect)) {
    if (mi[i, j] == 0) {
      l2_matrix_rt_no_latent_incorrect[i, j] <- l2_matrix_rt_no_latent[i, j]
    } else {
      l2_matrix_rt_no_latent_incorrect[i, j] <- NA
    }
  }
}

par(mfrow = c(1, 2))
hist(l1_matrix_rt_no_latent, nclass = 100)
hist(l2_matrix_rt_no_latent, nclass = 100)

mean(l1_matrix_rt_no_latent)
mean(l2_matrix_rt_no_latent)

par(mfrow = c(2, 2))
hist(l1_matrix_rt_no_latent_correct, nclass = 100)
hist(l1_matrix_rt_no_latent_incorrect, nclass = 100)
hist(l2_matrix_rt_no_latent_correct, nclass = 100)
hist(l2_matrix_rt_no_latent_incorrect, nclass = 100)

mean(l1_matrix_rt_no_latent_correct, na.rm = TRUE) # 16237.5
mean(l1_matrix_rt_no_latent_incorrect, na.rm = TRUE) # 26631.3
mean(l2_matrix_rt_no_latent_correct, na.rm = TRUE) # 165749.7
mean(l2_matrix_rt_no_latent_incorrect, na.rm = TRUE) # 379816.2

## loss function(L1, L2)_response
l1_matrix_no_latent <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:40) {
  prob <- sim_data[[i]][, 252:502]
  for (j in 1:251) {
    l1_matrix_no_latent[i, j] <- sum(abs(mi[i, j] - prob[, j]))
  }
}

l1_matrix_no_latent_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_no_latent_correct)) {
  for (j in 1:ncol(l1_matrix_no_latent_correct)) {
    if (mi[i, j] == 1) {
      l1_matrix_no_latent_correct[i, j] <- l1_matrix_no_latent[i, j]
    } else {
      l1_matrix_no_latent_correct[i, j] <- NA
    }
  }
}

l1_matrix_no_latent_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l1_matrix_no_latent_incorrect)) {
  for (j in 1:ncol(l1_matrix_no_latent_incorrect)) {
    if (mi[i, j] == 0) {
      l1_matrix_no_latent_incorrect[i, j] <- l1_matrix_no_latent[i, j]
    } else {
      l1_matrix_no_latent_incorrect[i, j] <- NA
    }
  }
}


l2_matrix_no_latent <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:40) {
  prob <- sim_data[[i]][, 252:502]
  for (j in 1:251) {
    l2_matrix_no_latent[i, j] <- sum((mi[i, j] - prob[, j])^2)
  }
}

l2_matrix_no_latent_correct <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_no_latent_correct)) {
  for (j in 1:ncol(l2_matrix_no_latent_correct)) {
    if (mi[i, j] == 1) {
      l2_matrix_no_latent_correct[i, j] <- l2_matrix_no_latent[i, j]
    } else {
      l2_matrix_no_latent_correct[i, j] <- NA
    }
  }
}

l2_matrix_no_latent_incorrect <- matrix(0, nrow = 40, ncol = 251)
for (i in 1:nrow(l2_matrix_no_latent_incorrect)) {
  for (j in 1:ncol(l2_matrix_no_latent_incorrect)) {
    if (mi[i, j] == 0) {
      l2_matrix_no_latent_incorrect[i, j] <- l2_matrix_no_latent[i, j]
    } else {
      l2_matrix_no_latent_incorrect[i, j] <- NA
    }
  }
}

par(mfrow = c(1, 2))
hist(l1_matrix_no_latent, nclass = 100)
hist(l2_matrix_no_latent, nclass = 100)

mean(l1_matrix_no_latent)
mean(l2_matrix_no_latent)

par(mfrow = c(2, 2))
hist(l1_matrix_no_latent_correct, nclass = 100)
hist(l1_matrix_no_latent_incorrect, nclass = 100)
hist(l2_matrix_no_latent_correct, nclass = 100)
hist(l2_matrix_no_latent_incorrect, nclass = 100)

mean(l1_matrix_no_latent_correct, na.rm = TRUE) # 1314.027
mean(l1_matrix_no_latent_incorrect, na.rm = TRUE) # 1282.058
mean(l2_matrix_no_latent_correct, na.rm = TRUE) # 862.7586
mean(l2_matrix_no_latent_incorrect, na.rm = TRUE) # 829.8821
