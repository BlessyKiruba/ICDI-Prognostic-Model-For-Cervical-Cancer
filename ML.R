mode = NULL
single_ml = NULL
alpha_for_Enet = NULL
direction_for_stepcox = NULL
double_ml1 = NULL
double_ml2 = NULL
nodesize = NULL
seed = NULL
cores_for_parallel = NULL

alpha_for_Enet <- 0.1
cores_for_parallel <- 6
direction_for_stepcox <- "both"

library(Matrix)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(mixOmics)
library(data.table)

result <- data.frame()
ml.res = list()
riskscore = list()
est_data <- as.data.frame(df_TCGA)
val_data_list <- list(
  TCGA_CESC = df_TCGA,
  CGCI_CC = df_CGI,
  GSE52903 = df_GSE)
pre_var <- colnames(est_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
est_dd$OS.time <- as.numeric(est_dd$OS.time)    


val_dd_list <- lapply(
  val_data_list,
  function(df) {
    df$OS.time <- as.numeric(df$OS.time)      # make OS.time numeric
    df[ , c("OS.time", "OS", pre_var), drop = FALSE]   # keep wanted cols
  }
)


rf_nodesize <- 5
seed <- 1

###################################
############# 1.1 RSF #############
###################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- "RSF"
result <- rbind(result, cc)
ml.res[["RSF"]] = fit
riskscore[["RSF"]] = rs

###################################
#######1.2 RSF + CoxBoost #########
###################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], 
                              est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                      2)]), trace = TRUE, start.penalty = 500, 
                              parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], 
                        est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                2)]), maxstepno = 500, K = 10, type = "verweij", 
                        penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, 
                                                "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                  penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            newdata = x[, -c(1, 2)], newtime = x[, 1], 
                                            newstatus = x[, 2], type = "lp")))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "CoxBoost")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "CoxBoost")]] = fit
  riskscore[[paste0("RSF + ", "CoxBoost")]] = rs
}

###################################
######### 1.3 RSF + Enet ##########
###################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, 
                    nfolds = 10)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "link", newx = as.matrix(x[, -c(1, 
                                                                                     2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("RSF + ", "Enet", "[α=", 
                       alpha, "]")
    result <- rbind(result, cc)
    ml.res[[paste0("RSF + ", "Enet", "[α=", alpha, 
                   "]")]] = fit
    riskscore[[paste0("RSF + ", "Enet", "[α=", 
                      alpha, "]")]] = rs
  }
}
###################################
############ 1.4 RSF + GBM ########
###################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
             data = est_dd2, distribution = "coxph", n.trees = 10000, 
             interaction.depth = 3, n.minobsinnode = 10, 
             shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
             data = est_dd2, distribution = "coxph", n.trees = best, 
             interaction.depth = 3, n.minobsinnode = 10, 
             shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            x, n.trees = best, type = "link")))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "GBM")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "GBM")]] = list(fit = fit, 
                                           best = best)
  riskscore[[paste0("RSF + ", "GBM")]] = rs
}

###################################
######## 1.1 RSF + Lasso ##########
###################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                  alpha = 1, type.measure = "C")
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "response", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "Lasso")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "Lasso")]] = fit
  riskscore[[paste0("RSF + ", "Lasso")]] = rs
}
######################################
####### 1.6 RSF  + plsRcox ###########
######################################
message("---1-6.RSF + plsRcox ---")
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, 
                                               rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                              nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, 
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "lp", newdata = x[, -c(1, 2)])))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "plsRcox")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "plsRcox")]] = fit
  riskscore[[paste0("RSF + ", "plsRcox")]] = rs
}

######################################
####### 1.7 RSF  + Ridge   ###########
######################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                  alpha = 0, type.measure = "C")
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "response", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "Ridge")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "Ridge")]] = fit
  riskscore[[paste0("RSF + ", "Ridge")]] = rs
}

######################################
####### 1.8 RSF  + StepCox ###########
######################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  for (direction in c("both", "backward", "forward")) {
    fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), 
                direction = direction)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict(fit, type = "risk", 
                                   newdata = x))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("RSF + ", "StepCox", "[", 
                       direction, "]")
    result <- rbind(result, cc)
    ml.res[[paste0("RSF + ", "StepCox", "[", direction, 
                   "]")]] = fit
    riskscore[[paste0("RSF + ", "StepCox", "[", 
                      direction, "]")]] = rs
  }
}
######################################
####### 1.9 RSF  + SuperPC ###########
######################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, 
               censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                                                                                  2)])
  set.seed(seed)
  fit <- superpc.train(data = data, type = "survival", 
                       s0.perc = 0.5)
  repeat {
    tryCatch({
      cv.fit <- superpc.cv(fit, data, n.threshold = 20, 
                           n.fold = 10, n.components = 3, min.features = 2, 
                           max.features = nrow(data$x), compute.fullcv = TRUE, 
                           compute.preval = TRUE)
      break
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      cat("Retrying...\n")
      Sys.sleep(1)
    })
  }
  rs <- lapply(val_dd_list2, function(w) {
    test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, 
                 censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                                                                        2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
    ])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[, 1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "SuperPC")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "SuperPC")]] = list(fit, 
                                               cv.fit)
  riskscore[[paste0("RSF + ", "SuperPC")]] = rs
}

##############################
#### 10 RSF + survival-SVM ###
##############################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, 
                    gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            x)$predicted))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "survival-SVM")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "survival-SVM")]] = fit
  riskscore[[paste0("RSF + ", "survival-SVM")]] = rs
}

##############################
########## 2. Enet ###########
##############################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
for (alpha in seq(0.1,0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, 
                  nfolds = 10)
  rs <- lapply(val_dd_list, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "link", newx = as.matrix(x[, -c(1, 
                                                                                   2)]), s = fit$lambda.min)))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("Enet", "[α=", alpha, "]")
  result <- rbind(result, cc)
  ml.res[[paste0("Enet", "[α=", alpha, "]")]] = fit
  riskscore[[paste0("Enet", "[α=", alpha, "]")]] = rs
}
##############################
######### 3.StepCox ##########
##############################
for (direction in c("both", "backward", "forward")) {
  message(paste0("---3.StepCox", "[", direction, 
                 "]---"))
  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), 
              direction = direction)
  rs <- lapply(val_dd_list, function(x) {
    cbind(x[, 1:2], RS = predict(fit, type = "risk", 
                                 newdata = x))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("StepCox", "[", direction, 
                     "]")
  result <- rbind(result, cc)
  ml.res[[paste0("StepCox", "[", direction, "]")]] = fit
  riskscore[[paste0("StepCox", "[", direction, 
                    "]")]] = rs
}
direction = "both"
if (T) {
  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), 
              direction = direction)
  rid <- names(coef(fit))
  if (length(rid) > 1) {
    est_dd2 <- est_dd[, c("OS.time", "OS", 
                          rid)]
    val_dd_list2 <- lapply(val_dd_list, 
                           function(x) {
                             x[, c("OS.time", "OS", rid)]
                           })
    set.seed(seed)
    message(paste0("---3.StepCox", "[", direction, 
                   "]", " + CoxBoost ---"))
    pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], 
                                est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                        2)]), trace = TRUE, start.penalty = 500, 
                                parallel = T)
    cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], 
                          est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                  2)]), maxstepno = 500, K = 10, type = "verweij", 
                          penalty = pen$penalty)
    fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, 
                                                  "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              newdata = x[, -c(1, 2)], newtime = x[, 
                                                                                   1], newstatus = x[, 2], type = "lp")))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + CoxBoost")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + CoxBoost")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + CoxBoost")]] = rs
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    for (alpha in seq(0.1, 0.9, 0.1)) {
      set.seed(seed)
      message(paste0("--- 3.StepCox", "[", direction, 
                     "]", " + Enet", "[α=", alpha, "] ---"))
      fit = cv.glmnet(x1, x2, family = "cox", 
                      alpha = alpha, nfolds = 10)
      rs <- lapply(val_dd_list2, function(x) {
        cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                                type = "link", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
      })
      cc <- data.frame(Cindex = sapply(rs, function(x) {
        as.numeric(summary(coxph(Surv(OS.time, 
                                      OS) ~ RS, x))$concordance[1])
      })) %>% rownames_to_column("ID")
      cc$Model <- paste0("StepCox", "[", direction, 
                         "]", " + Enet", "[α=", alpha, "]")
      result <- rbind(result, cc)
      ml.res[[paste0("StepCox", "[", direction, 
                     "]", " + Enet", "[α=", alpha, "]")]] = fit
      riskscore[[paste0("StepCox", "[", direction, 
                        "]", " + Enet", "[α=", alpha, "]")]] = rs
    }
    set.seed(seed)
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + GBM ---"))
    fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
               data = est_dd2, distribution = "coxph", 
               n.trees = 10000, interaction.depth = 3, 
               n.minobsinnode = 10, shrinkage = 0.001, 
               cv.folds = 10, n.cores = cores_for_parallel)
    best <- which.min(fit$cv.error)
    set.seed(seed)
    fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
               data = est_dd2, distribution = "coxph", 
               n.trees = best, interaction.depth = 3, n.minobsinnode = 10, 
               shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              x, n.trees = best, type = "link")))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + GBM")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + GBM")]] = list(fit = fit, best = best)
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + GBM")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + Lasso ---"))
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 1, type.measure = "C")
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "response", newx = as.matrix(x[, 
                                                                                    -c(1, 2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + Lasso")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + Lasso")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + Lasso")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + plsRcox ---"))
    set.seed(seed)
    cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, 
                                                 rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                                nt = 10, verbose = FALSE)
    fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, 
                   event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "lp", newdata = x[, -c(1, 2)])))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + plsRcox")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + plsRcox")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + plsRcox")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + Ridge ---"))
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 0, type.measure = "C")
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "response", newx = as.matrix(x[, 
                                                                                    -c(1, 2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + Ridge")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + Ridge")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + Ridge")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + RSF ---"))
    set.seed(seed)
    fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, 
                 ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                 importance = T, proximity = T, forest = T, 
                 seed = seed)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + RSF")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + RSF")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + RSF")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + SuperPC ---"))
    data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, 
                 censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                                                                                    2)])
    set.seed(seed)
    fit <- superpc.train(data = data, type = "survival", 
                         s0.perc = 0.5)
    repeat {
      tryCatch({
        cv.fit <- superpc.cv(fit, data, n.threshold = 20, 
                             n.fold = 10, n.components = 3, min.features = 2, 
                             max.features = nrow(data$x), compute.fullcv = TRUE, 
                             compute.preval = TRUE)
        break
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        cat("Retrying...\n")
        Sys.sleep(1)
      })
    }
    rs <- lapply(val_dd_list2, function(w) {
      test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, 
                   censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                                                                          2)])
      ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
      ])], n.components = 1)
      rr <- as.numeric(ff$v.pred)
      rr2 <- cbind(w[, 1:2], RS = rr)
      return(rr2)
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + SuperPC")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + SuperPC")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + survival-SVM ---"))
    fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, 
                      gamma.mu = 1)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              x)$predicted))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + survival-SVM")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + survival-SVM")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + survival-SVM")]] = rs
  }
  else {
    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
  }
}
direction = "backward"
if (T) {
  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), 
              direction = direction)
  rid <- names(coef(fit))
  if (length(rid) > 1) {
    est_dd2 <- est_dd[, c("OS.time", "OS", 
                          rid)]
    val_dd_list2 <- lapply(val_dd_list, 
                           function(x) {
                             x[, c("OS.time", "OS", rid)]
                           })
    set.seed(seed)
    message(paste0("---3.StepCox", "[", direction, 
                   "]", " + CoxBoost ---"))
    pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], 
                                est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                        2)]), trace = TRUE, start.penalty = 500, 
                                parallel = T)
    cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], 
                          est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                  2)]), maxstepno = 500, K = 10, type = "verweij", 
                          penalty = pen$penalty)
    fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, 
                                                  "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              newdata = x[, -c(1, 2)], newtime = x[, 
                                                                                   1], newstatus = x[, 2], type = "lp")))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + CoxBoost")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + CoxBoost")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + CoxBoost")]] = rs
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    for (alpha in seq(0.1, 0.9, 0.1)) {
      set.seed(seed)
      message(paste0("--- 3.StepCox", "[", direction, 
                     "]", " + Enet", "[α=", alpha, "] ---"))
      fit = cv.glmnet(x1, x2, family = "cox", 
                      alpha = alpha, nfolds = 10)
      rs <- lapply(val_dd_list2, function(x) {
        cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                                type = "link", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
      })
      cc <- data.frame(Cindex = sapply(rs, function(x) {
        as.numeric(summary(coxph(Surv(OS.time, 
                                      OS) ~ RS, x))$concordance[1])
      })) %>% rownames_to_column("ID")
      cc$Model <- paste0("StepCox", "[", direction, 
                         "]", " + Enet", "[α=", alpha, "]")
      result <- rbind(result, cc)
      ml.res[[paste0("StepCox", "[", direction, 
                     "]", " + Enet", "[α=", alpha, "]")]] = fit
      riskscore[[paste0("StepCox", "[", direction, 
                        "]", " + Enet", "[α=", alpha, "]")]] = rs
    }
    set.seed(seed)
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + GBM ---"))
    fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
               data = est_dd2, distribution = "coxph", 
               n.trees = 10000, interaction.depth = 3, 
               n.minobsinnode = 10, shrinkage = 0.001, 
               cv.folds = 10, n.cores = cores_for_parallel)
    best <- which.min(fit$cv.error)
    set.seed(seed)
    fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
               data = est_dd2, distribution = "coxph", 
               n.trees = best, interaction.depth = 3, n.minobsinnode = 10, 
               shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              x, n.trees = best, type = "link")))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + GBM")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + GBM")]] = list(fit = fit, best = best)
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + GBM")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + Lasso ---"))
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 1, type.measure = "C")
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "response", newx = as.matrix(x[, 
                                                                                    -c(1, 2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + Lasso")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + Lasso")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + Lasso")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + plsRcox ---"))
    set.seed(seed)
    cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, 
                                                 rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                                nt = 10, verbose = FALSE)
    fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, 
                   event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "lp", newdata = x[, -c(1, 2)])))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + plsRcox")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + plsRcox")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + plsRcox")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + Ridge ---"))
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 0, type.measure = "C")
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "response", newx = as.matrix(x[, 
                                                                                    -c(1, 2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + Ridge")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + Ridge")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + Ridge")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + RSF ---"))
    set.seed(seed)
    fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, 
                 ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                 importance = T, proximity = T, forest = T, 
                 seed = seed)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + RSF")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + RSF")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + RSF")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + SuperPC ---"))
    data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, 
                 censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                                                                                    2)])
    set.seed(seed)
    fit <- superpc.train(data = data, type = "survival", 
                         s0.perc = 0.5)
    repeat {
      tryCatch({
        cv.fit <- superpc.cv(fit, data, n.threshold = 20, 
                             n.fold = 10, n.components = 3, min.features = 2, 
                             max.features = nrow(data$x), compute.fullcv = TRUE, 
                             compute.preval = TRUE)
        break
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        cat("Retrying...\n")
        Sys.sleep(1)
      })
    }
    rs <- lapply(val_dd_list2, function(w) {
      test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, 
                   censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                                                                          2)])
      ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
      ])], n.components = 1)
      rr <- as.numeric(ff$v.pred)
      rr2 <- cbind(w[, 1:2], RS = rr)
      return(rr2)
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + SuperPC")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + SuperPC")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + survival-SVM ---"))
    fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, 
                      gamma.mu = 1)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              x)$predicted))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + survival-SVM")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + survival-SVM")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + survival-SVM")]] = rs
  }
  else {
    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
  }
}
direction = "forward"
if (T) {
  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), 
              direction = direction)
  rid <- names(coef(fit))
  if (length(rid) > 1) {
    est_dd2 <- est_dd[, c("OS.time", "OS", 
                          rid)]
    val_dd_list2 <- lapply(val_dd_list, 
                           function(x) {
                             x[, c("OS.time", "OS", rid)]
                           })
    set.seed(seed)
    message(paste0("---3.StepCox", "[", direction, 
                   "]", " + CoxBoost ---"))
    pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], 
                                est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                        2)]), trace = TRUE, start.penalty = 500, 
                                parallel = T)
    cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], 
                          est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                  2)]), maxstepno = 500, K = 10, type = "verweij", 
                          penalty = pen$penalty)
    fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, 
                                                  "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              newdata = x[, -c(1, 2)], newtime = x[, 
                                                                                   1], newstatus = x[, 2], type = "lp")))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + CoxBoost")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + CoxBoost")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + CoxBoost")]] = rs
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    for (alpha in seq(0.1, 0.9, 0.1)) {
      set.seed(seed)
      message(paste0("--- 3.StepCox", "[", direction, 
                     "]", " + Enet", "[α=", alpha, "] ---"))
      fit = cv.glmnet(x1, x2, family = "cox", 
                      alpha = alpha, nfolds = 10)
      rs <- lapply(val_dd_list2, function(x) {
        cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                                type = "link", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
      })
      cc <- data.frame(Cindex = sapply(rs, function(x) {
        as.numeric(summary(coxph(Surv(OS.time, 
                                      OS) ~ RS, x))$concordance[1])
      })) %>% rownames_to_column("ID")
      cc$Model <- paste0("StepCox", "[", direction, 
                         "]", " + Enet", "[α=", alpha, "]")
      result <- rbind(result, cc)
      ml.res[[paste0("StepCox", "[", direction, 
                     "]", " + Enet", "[α=", alpha, "]")]] = fit
      riskscore[[paste0("StepCox", "[", direction, 
                        "]", " + Enet", "[α=", alpha, "]")]] = rs
    }
    set.seed(seed)
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + GBM ---"))
    fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
               data = est_dd2, distribution = "coxph", 
               n.trees = 10000, interaction.depth = 3, 
               n.minobsinnode = 10, shrinkage = 0.001, 
               cv.folds = 10, n.cores = cores_for_parallel)
    best <- which.min(fit$cv.error)
    set.seed(seed)
    fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
               data = est_dd2, distribution = "coxph", 
               n.trees = best, interaction.depth = 3, n.minobsinnode = 10, 
               shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              x, n.trees = best, type = "link")))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + GBM")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + GBM")]] = list(fit = fit, best = best)
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + GBM")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + Lasso ---"))
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 1, type.measure = "C")
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "response", newx = as.matrix(x[, 
                                                                                    -c(1, 2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + Lasso")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + Lasso")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + Lasso")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + plsRcox ---"))
    set.seed(seed)
    cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, 
                                                 rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                                nt = 10, verbose = FALSE)
    fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, 
                   event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "lp", newdata = x[, -c(1, 2)])))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + plsRcox")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + plsRcox")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + plsRcox")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + Ridge ---"))
    x1 <- as.matrix(est_dd2[, rid])
    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 0, type.measure = "C")
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "response", newx = as.matrix(x[, 
                                                                                    -c(1, 2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + Ridge")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + Ridge")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + Ridge")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + RSF ---"))
    set.seed(seed)
    fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, 
                 ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                 importance = T, proximity = T, forest = T, 
                 seed = seed)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + RSF")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + RSF")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + RSF")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + SuperPC ---"))
    data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, 
                 censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                                                                                    2)])
    set.seed(seed)
    fit <- superpc.train(data = data, type = "survival", 
                         s0.perc = 0.5)
    repeat {
      tryCatch({
        cv.fit <- superpc.cv(fit, data, n.threshold = 20, 
                             n.fold = 10, n.components = 3, min.features = 2, 
                             max.features = nrow(data$x), compute.fullcv = TRUE, 
                             compute.preval = TRUE)
        break
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        cat("Retrying...\n")
        Sys.sleep(1)
      })
    }
    rs <- lapply(val_dd_list2, function(w) {
      test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, 
                   censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                                                                          2)])
      ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
      ])], n.components = 1)
      rr <- as.numeric(ff$v.pred)
      rr2 <- cbind(w[, 1:2], RS = rr)
      return(rr2)
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + SuperPC")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + SuperPC")]] = rs
    message(paste0("--- 3.StepCox", "[", direction, 
                   "]", " + survival-SVM ---"))
    fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, 
                      gamma.mu = 1)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              x)$predicted))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("StepCox", "[", direction, 
                       "]", " + survival-SVM")
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox", "[", direction, 
                   "]", " + survival-SVM")]] = fit
    riskscore[[paste0("StepCox", "[", direction, 
                      "]", " + survival-SVM")]] = rs
  }
  else {
    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
  }
}


###################################
######### 4.1 CoxBoost ############
###################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                          newdata = x[, -c(1, 2)], newtime = x[, 1], 
                                          newstatus = x[, 2], type = "lp")))
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- paste0("CoxBoost")
result <- rbind(result, cc)
ml.res[[paste0("CoxBoost")]] = fit
riskscore[[paste0("CoxBoost")]] = rs

###################################
######### 4. 2 CoxBoost + ENET#####
###################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    message(paste0("--- 4.CoxBoost", " + Enet", 
                   "[α=", alpha, "] ---"))
    fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, 
                    nfolds = 10)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                              type = "link", newx = as.matrix(x[, -c(1, 
                                                                                     2)]), s = fit$lambda.min)))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("CoxBoost", " + Enet", 
                       "[α=", alpha, "]")
    result <- rbind(result, cc)
    ml.res[[paste0("CoxBoost", " + Enet", "[α=", 
                   alpha, "]")]] = fit
    riskscore[[paste0("CoxBoost", " + Enet", "[α=", 
                      alpha, "]")]] = rs
  }
}
###################################
#########4.3 CoxBoost + GBM########
###################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
             data = est_dd2, distribution = "coxph", n.trees = 10000, 
             interaction.depth = 3, n.minobsinnode = 10, 
             shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
             data = est_dd2, distribution = "coxph", n.trees = best, 
             interaction.depth = 3, n.minobsinnode = 10, 
             shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            x, n.trees = best, type = "link")))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("CoxBoost + ", "GBM")
  result <- rbind(result, cc)
  ml.res[[paste0("CoxBoost + ", "GBM")]] = list(fit = fit, 
                                                best = best)
  riskscore[[paste0("CoxBoost + ", "GBM")]] = rs
}
###################################
#########4.4 CoxBoost + Lasso######
###################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                  alpha = 1, type.measure = "C")
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "response", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("CoxBoost + ", "Lasso")
  result <- rbind(result, cc)
  ml.res[[paste0("CoxBoost + ", "Lasso")]] = fit
  riskscore[[paste0("CoxBoost + ", "Lasso")]] = rs
}
###################################
######### 4.5 CoxBoost + plsRcox###
###################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, 
                                               rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                              nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, 
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "lp", newdata = x[, -c(1, 2)])))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("CoxBoost + ", "plsRcox")
  result <- rbind(result, cc)
  ml.res[[paste0("CoxBoost + ", "plsRcox")]] = fit
  riskscore[[paste0("CoxBoost + ", "plsRcox")]] = rs
}

###################################
#########4.6 CoxBoost + Ridge######
###################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                  alpha = 0, type.measure = "C")
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "response", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("CoxBoost + ", "Ridge")
  result <- rbind(result, cc)
  ml.res[[paste0("CoxBoost + ", "Ridge")]] = fit
  riskscore[[paste0("CoxBoost + ", "Ridge")]] = rs
}

###################################
######### 4.7 CoxBoost + StepCox###
###################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  for (direction in c("both", "backward", "forward")) {
    message(paste0("--- 4.CoxBoost + ", "StepCox", 
                   "[", direction, "] ---"))
    fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), 
                direction = direction)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict(fit, type = "risk", 
                                   newdata = x))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("CoxBoost + ", "StepCox", 
                       "[", direction, "]")
    result <- rbind(result, cc)
    ml.res[[paste0("CoxBoost + ", "StepCox", "[", 
                   direction, "]")]] = fit
    riskscore[[paste0("CoxBoost + ", "StepCox", 
                      "[", direction, "]")]] = rs
  }
}

###################################
####### 4.8 CoxBoost + SuperPC#####
###################################

set.seed(seed)
message(paste0("--- 4.CoxBoost + ", "SuperPC ---"))
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, 
               censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                                                                                  2)])
  set.seed(seed)
  fit <- superpc.train(data = data, type = "survival", 
                       s0.perc = 0.5)
  repeat {
    tryCatch({
      cv.fit <- superpc.cv(fit, data, n.threshold = 20, 
                           n.fold = 10, n.components = 3, min.features = 2, 
                           max.features = nrow(data$x), compute.fullcv = TRUE, 
                           compute.preval = TRUE)
      break
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      cat("Retrying...\n")
      Sys.sleep(1)
    })
  }
  rs <- lapply(val_dd_list2, function(w) {
    test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, 
                 censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                                                                        2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
    ])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[, 1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("CoxBoost + ", "SuperPC")
  result <- rbind(result, cc)
  ml.res[[paste0("CoxBoost + ", "SuperPC")]] = list(fit = fit, 
                                                    cv.fit = cv.fit)
  riskscore[[paste0("CoxBoost + ", "SuperPC")]] = rs
}
###################################
#### 4.9 CoxBoost + survival-SVM###
###################################
message(paste0("--- 4.CoxBoost + ", "survival-SVM ---"))
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], 
                            est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, 
                                                  "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                      K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], 
                as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, 
                    gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            x)$predicted))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("CoxBoost + ", "survival-SVM")
  result <- rbind(result, cc)
  ml.res[[paste0("CoxBoost + ", "survival-SVM")]] = fit
  riskscore[[paste0("CoxBoost + ", "survival-SVM")]] = rs
}

#################################
############# 5.plsRcox #########
#################################
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd[, pre_var], 
                                 time = est_dd$OS.time, status = est_dd$OS), 
                            nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd[, pre_var], time = est_dd$OS.time, 
               event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                          type = "lp", newdata = x[, -c(1, 2)])))
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- paste0("plsRcox")
result <- rbind(result, cc)
ml.res[[paste0("plsRcox")]] = fit
riskscore[[paste0("plsRcox")]] = rs
##########################
####### 6.superpc ########
##########################

data <- list(x = t(est_dd[, -c(1, 2)]), y = est_dd$OS.time, 
             censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 
                                                                              2)])
set.seed(seed)
fit <- superpc.train(data = data, type = "survival", 
                     s0.perc = 0.5)
repeat {
  tryCatch({
    cv.fit <- superpc.cv(fit, data, n.threshold = 20, 
                         n.fold = 10, n.components = 3, min.features = 2, 
                         max.features = nrow(data$x), compute.fullcv = TRUE, 
                         compute.preval = TRUE)
    break
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    cat("Retrying...\n")
    Sys.sleep(1)
  })
}
rs <- lapply(val_dd_list, function(w) {
  test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, 
               censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                                                                      2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
  ])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- paste0("SuperPC")
result <- rbind(result, cc)
ml.res[[paste0("SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
riskscore[[paste0("SuperPC")]] = rs
############################
####### 7.GBM ##############
###########################
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd, 
           distribution = "coxph", n.trees = 10000, interaction.depth = 3, 
           n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, 
           n.cores = cores_for_parallel)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd, 
           distribution = "coxph", n.trees = best, interaction.depth = 3, 
           n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, 
           n.cores = cores_for_parallel)
rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                          x, n.trees = best, type = "link")))
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- paste0("GBM")
result <- rbind(result, cc)
ml.res[[paste0("GBM")]] = list(fit = fit, best = best)
riskscore[[paste0("GBM")]] = rs

###########################
#### 8.survivalsvm ########
###########################
fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd, 
                  gamma.mu = 1)
rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                          x)$predicted))
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- paste0("survival - SVM")
result <- rbind(result, cc)
ml.res[[paste0("survival - SVM")]] = fit
riskscore[[paste0("survival - SVM")]] = rs

##############################
########### 9.Ridge ##########
##############################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = glmnet(x1, x2, family = "cox", alpha = 0, 
             lambda = NULL)
cv.fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                   type.measure = "C")
rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                          type = "response", newx = as.matrix(x[, -c(1, 
                                                                                     2)]), s = cv.fit$lambda.min)))
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- paste0("Ridge")
result <- rbind(result, cc)
ml.res[[paste0("Ridge")]] = list(fit = fit, cv.fit = cv.fit)
riskscore[[paste0("Ridge")]] = rs



##############################
#########10.Lasso#############
##############################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                alpha = 1, type.measure = "C")
rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                          type = "response", newx = as.matrix(x[, -c(1, 
                                                                                     2)]), s = fit$lambda.min)))
})
cc <- data.frame(Cindex = sapply(rs, function(x) {
  as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                             RS, x))$concordance[1])
})) %>% rownames_to_column("ID")
cc$Model <- paste0("Lasso")
result <- rbind(result, cc)
ml.res[[paste0("Lasso")]] = fit
riskscore[[paste0("Lasso")]] = rs

##############################
###10.Lasso + CoxBoost #######
##############################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                alpha = 1, type.measure = "C")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min")
rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 
                                             0)]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], 
                              est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                      2)]), trace = TRUE, start.penalty = 500, 
                              parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], 
                        est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                                                                2)]), maxstepno = 500, K = 10, type = "verweij", 
                        penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, 
                                                "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                  penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            newdata = x[, -c(1, 2)], newtime = x[, 1], 
                                            newstatus = x[, 2], type = "lp")))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("Lasso + ", "CoxBoost")
  result <- rbind(result, cc)
  ml.res[[paste0("Lasso + ", "CoxBoost")]] = fit
  riskscore[[paste0("Lasso + ", "CoxBoost")]] = rs
}

#####################################
######### 10.Lasso + GBM#############
#####################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                alpha = 1, type.measure = "C")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min")
rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 
                                             0)]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
             data = est_dd2, distribution = "coxph", n.trees = 10000, 
             interaction.depth = 3, n.minobsinnode = 10, 
             shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS) ~ ., 
             data = est_dd2, distribution = "coxph", n.trees = best, 
             interaction.depth = 3, n.minobsinnode = 10, 
             shrinkage = 0.001, cv.folds = 10, n.cores = cores_for_parallel)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            x, n.trees = best, type = "link")))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("Lasso + ", "GBM")
  result <- rbind(result, cc)
  riskscore[[paste0("Lasso + ", "GBM")]] = rs
  ml.res[[paste0("Lasso + ", "GBM")]] = list(fit = fit, 
                                             best = best)
}

#################################################
############# 10.Lasso + plsRcox ################
#################################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                alpha = 1, type.measure = "C")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min")
rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 
                                             0)]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, 
                                               rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                              nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, 
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "lp", newdata = x[, -c(1, 2)])))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("Lasso + ", "plsRcox")
  result <- rbind(result, cc)
  ml.res[[paste0("Lasso + ", "plsRcox")]] = fit
  riskscore[[paste0("Lasso + ", "plsRcox")]] = rs
}}

###########################################
###############10.Lasso + ", "RSF #########
###########################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                alpha = 1, type.measure = "C")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min")
rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 
                                             0)]
rid <- rid[-1]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, 
               ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
               importance = T, proximity = T, forest = T, 
               seed = seed)
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("Lasso", " + RSF")
  result <- rbind(result, cc)
  riskscore[[paste0("Lasso + ", "RSF")]] = rs
  ml.res[[paste0("Lasso", " + RSF")]] = fit
}
###################################
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ####################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                alpha = 1, type.measure = "C")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min")
rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 
                                             0)]
if (length(rid) > 1) {
  est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(val_dd_list, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  for (direction in c("both", "backward", "forward")) {
    message(paste0("--- 10.Lasso + ", "StepCox", 
                   "[", direction, "] ---"))
    fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), 
                direction = direction)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict(fit, type = "risk", 
                                   newdata = x))
    })
    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                 RS, x))$concordance[1])
    })) %>% rownames_to_column("ID")
    cc$Model <- paste0("Lasso + ", "StepCox", 
                       "[", direction, "]")
    result <- rbind(result, cc)
    ml.res[[paste0("Lasso + ", "StepCox", "[", 
                   direction, "]")]] = fit
    riskscore[[paste0("Lasso + ", "StepCox", "[", 
                      direction, "]")]] = rs
    
    ###############10.Lasso + ", "superPC############
    x1 <- as.matrix(est_dd[, pre_var])
    x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 1, type.measure = "C")
    fit$lambda.min
    myCoefs <- coef(fit, s = "lambda.min")
    rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 
                                                 0)]
    if (length(rid) > 1) {
      est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
      val_dd_list2 <- lapply(val_dd_list, 
                             function(x) {
                               x[, c("OS.time", "OS", rid)]
                             })
      data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, 
                   censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                                                                                      2)])
      set.seed(seed)
      fit <- superpc.train(data = data, type = "survival", 
                           s0.perc = 0.5)
      repeat {
        tryCatch({
          cv.fit <- superpc.cv(fit, data, n.threshold = 20, 
                               n.fold = 10, n.components = 3, min.features = 2, 
                               max.features = nrow(data$x), compute.fullcv = TRUE, 
                               compute.preval = TRUE)
          break
        }, error = function(e) {
          cat("Error:", conditionMessage(e), "\n")
          cat("Retrying...\n")
          Sys.sleep(1)
        })
      }
      rs <- lapply(val_dd_list2, function(w) {
        test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, 
                     censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                                                                            2)])
        ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
        ])], n.components = 1)
        rr <- as.numeric(ff$v.pred)
        rr2 <- cbind(w[, 1:2], RS = rr)
        return(rr2)
      })
      cc <- data.frame(Cindex = sapply(rs, function(x) {
        as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                   RS, x))$concordance[1])
      })) %>% rownames_to_column("ID")
      cc$Model <- paste0("Lasso + ", "SuperPC")
      result <- rbind(result, cc)
      ml.res[[paste0("Lasso + ", "SuperPC")]] = list(fit = fit, 
                                                     cv.fit = cv.fit)
      rs = returnIDtoRS(rs.table.list = rs, rawtableID = val_dd_list)
      riskscore[[paste0("Lasso + ", "SuperPC")]] = rs
    }
    ###########10.Lasso +survival-SVM#########
    #########################################
    x1 <- as.matrix(est_dd[, pre_var])
    x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
    set.seed(seed)
    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                    alpha = 1, type.measure = "C")
    fit$lambda.min
    myCoefs <- coef(fit, s = "lambda.min")
    rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 
                                                 0)]
    if (length(rid) > 1) {
      est_dd2 <- est_dd[, c("OS.time", "OS", rid)]
      val_dd_list2 <- lapply(val_dd_list, 
                             function(x) {
                               x[, c("OS.time", "OS", rid)]
                             })
      fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, 
                        gamma.mu = 1)
      rs <- lapply(val_dd_list2, function(x) {
        cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                                x)$predicted))
      })
      cc <- data.frame(Cindex = sapply(rs, function(x) {
        as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                                   RS, x))$concordance[1])
      })) %>% rownames_to_column("ID")
      cc$Model <- paste0("Lasso + ", "survival-SVM")
      result <- rbind(result, cc)
      ml.res[[paste0("Lasso + ", "survival-SVM")]] = fit
      riskscore[[paste0("Lasso + ", "survival-SVM")]] = rs
    }
    
    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, 
                             riskscore = riskscore, Sig.genes = pre_var)
    
    write_xlsx(result, "ML_prognostic.xlsx")