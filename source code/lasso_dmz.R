# set working directory
setwd("Z:/Projects/Martina_Berenil/Baseline model/raw data/")

# load data
data <- read.csv('RNasedata_refine.csv', header = TRUE)


# Creat the evaluation function: eval_results, which contained RMSE and Rsquare
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square)
}

x <- as.matrix(data[,-c(1,2)])
y <- as.matrix(data[2])

# lasso regression
library(glmnet)
set.seed(1)
lambdas <- 10^seq(2, -10, length = 1000)

# use sv.glmnet to find the best lambda for lasso from 5-fold cv
lasso_reg <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
plot(lasso_reg)

# plot the shrinkage graph with multiple lambda values
lasso_model <- glmnet(x, y, alpha = 1, nlambda =100,standardize = TRUE)
print(lasso_model)
p1 <- plot(lasso_model,xvar="lambda",label = T, lwd=4,cex.lab= 2,cex.axis=2,xlim = c(-4.5,0.5), ylim=c(-20,20))
p1.lty=2
box(lwd=4)

# find the lambda with lowest mean-squrared error from cv
lambda_best_lasso <- lasso_reg$lambda.min
lambda_best_lasso

# build the lasso regression model using selected descriptors, using refined lambda
lasso_model <- glmnet(x, y, alpha = 1, lambda =exp(-8),standardize = TRUE)
summary(lasso_model)

# find the non-zero coefficients and their names
lasso.coef  <-  predict(lasso_model,type="coefficients")
lasso.coef
lasso.coef[lasso.coef!=0] 
lasso_nonzerocoef <- predict(lasso_model,type="nonzero")
lasso_nonzerocoef
colnames(data[,lasso_nonzerocoef$s0+1])



# use stepwise further select the model
library(car)
data_step <- data[,append(lasso_nonzerocoef$s0+1,2,0)]
mdl_null = lm(RNase~1, data=data_step)
summary(mdl_null)

mdl_full = lm(RNase ~ ., data=data_step)
summary(mdl_full)

mdl_stepwise = step(mdl_null, scope = formula(c(mdl_full, mdl_null)), direction="both", trace=1, criteria="BIC", k=log(nrow(data)))
summary(mdl_stepwise)


# exhaustively search for all combinations
# m = number of features in the model, data_step contains all non-zero descriptor candidates from above lasso, "results" summarizes all results
data_step <- data[,append(lasso_nonzerocoef$s0+1,2,0)]

# exhaustively search for all combinations
data_step <- data[-1]
m <- 4
idx <- combn(rep(1:(length(data_step)-1)),m)
results <- NULL
idrows <- NULL
start_time <- Sys.time()
for (i in 1:ncol(idx)) {

  data_exhau <- data_step[,append(idx[,i]+1,1,0)]
  mdl_exhau <- lm(RNase~.,data=data_exhau)
  
  fitted <- mdl_exhau$fitted.values
  a <- eval_results(data$RNase,fitted,data)
    
    if (a[2]>=0.6) {
      result <- data.frame(train=a)
      results<- rbind(results,result)
      idrows <- rbind(idrows, i)
    }
    
    if (i%%10000 == 0){
      print(paste("This is the end of the run: ", i))
      time_passed <- Sys.time()-start_time
      print(paste("Time passed: ", time_passed))
      start_time <- Sys.time()
    }
}

# idrows find all candidates with top performance, and print out the model summary for statistical significance check as well as loocv
id_top <- idrows[which(results$train.Rsquare>=.75)]
results_wloocv <- NULL

for (val in id_top) {
  data_exhau <- data_step[,append(idx[,val]+1,1,0)]
  mdl_exhau <- lm(RNase~.,data=data_exhau)
  fitted <- mdl_exhau$fitted.values
  a <- eval_results(data$RNase,fitted,data)[2]
  s <- summary(mdl_exhau)
  print(s)
  yhat_loocv <- NULL
  for (i in 1:nrow(data)){
    data_loocv_test <- data_exhau[i,]
    data_loocv_train <- data_exhau[-i,]
    mdl_loocv <- lm(RNase~.,data=data_loocv_train)
    predict_loocv <- predict(mdl_loocv,newdata = data_loocv_test)
    yhat_loocv <- rbind(yhat_loocv, predict_loocv)
  }
  loocv_r2 <- eval_results(data$RNase,yhat_loocv,data)$Rsquare
  result <- data.frame(train=a, loocv=loocv_r2)
  results_wloocv <- rbind(results_wloocv,result)
  cat("loocv:", loocv_r2)
 }
# plot the curve for the top model

library(ggplot2)
# load the model

mdl <- lm(formula = "RNase~1+FP_460+FP_694+FP_1297+FP_1508", 
          data = trainingset)
  mdl <- lm(formula = "RNase~1+FP_141+FP_360+FP_380+FP_829+FP_1386", 
          data = trainingset)
mdl <- lm(formula = "RNase~1+AM1_dipole+BCUT_SMR_2+std_dim2+FP_1172+FP_1601", 
          data = trainingset)
summary(mdl)
predict <- predict(mdl,newdata = data)
id  <-  numeric(59)
id[-trainid] <- 1
data_plot <- cbind(predict,data$RNase,id)
colnames(data_plot) <- c("predict", "obs","id")

ggplot(as.data.frame(data_plot), aes(x=obs,y=predict))+
  ggtitle(expression("Baseline model of lnk"[off]*"")) +
  xlab(expression("Observed lnk"[off]*"")) + ylab(expression("Predicted lnk"[off]*""))+
  # THE DATA POINT
  geom_point(aes(color = factor(id)),size = 5,alpha =1) +
  xlim(min(data$RNase)-2,max(data$RNase)+2)+
  ylim(min(data$RNase)-2,max(data$RNase)+2)+
  scale_color_manual(labels = c("Training set", "Test set"), values = c("dodgerblue", "red2"))+
  
  # title
  theme_bw()+
  theme(axis.ticks.length=unit(.4,"lines"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size=20), 
        axis.title = element_text(size = 25,face = 'bold'),title =element_text(size = 25,face = 'bold') )+
# legend  
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size=20, face="bold"))+
  theme(legend.position = c(0.80, 0.1))+
# rec
  theme(panel.background = element_rect(colour = "black", size = 3.5))+
# ref line  
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=1.5)

ggsave("RNase_FP+MOE.tiff", units="in", width=8, height=8, dpi=600)




# loop search for opt lambda
results <- NULL
lambda <- seq(0.001,5, length=10000)
for (val in lambda) {
  lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda = val,standardize = TRUE)
  summary(lasso_model)
  lasso.coef = predict(lasso_model,type="coefficients")
  
  lasso_fittings <- predict(lasso_model, s = val, newx = x_train)
  lasso_predictions <- predict(lasso_model, s = val, newx = x_test)
  R2_test <- as.matrix(eval_results(y_test, lasso_predictions, testset)[2])
  R2_train <- as.matrix(eval_results(y_train, lasso_fittings, trainingset)[2])
  idx_nonzero <- which(lasso.coef!=0)
  len <- length(idx_nonzero)
  feature <- colnames(trainingset[,idx_nonzero[-1]])
  feature <- paste0(feature,collapse = ",")
  result <- data.frame(R1=R2_test, 
                        R2=R2_train,
                        len=len-1,
                       lambda=val,
                       Variable=feature)
  results<- rbind(results,result)
}

plot(x = 1,
     type = "n",
     xlim = c(0, 5), 
     ylim = c(0, 1),
     pch = 16,
     xlab = "lambda", 
     ylab = "R2",
     main = "Training process")
# add training and test distribution
points(x=results[,4],y=results[,1],pch = 16,
       col = transparent("red", trans.val = .2),
       cex = 0.5)
points(x = results[,4],
       y = results[,2],
       pch = 16,
       col = transparent("steelblue3", trans.val = .5),
       cex = 0.5)


legend("bottomright",
       legend = c("training", "test"),
       col = transparent(c('steelblue3', 'red'), .2),
       pch = c(16, 16),cex = 2,
       bg = "white")

lasso_pred <- predict(lasso,as.data.frame(x_test))
lasso_fitted <- predict(lasso,trainingset[-1])
a <- eval_results(y_test,lasso_pred, testset)
a
b <- eval_results(y,lasso_fitted, trainingset)
b

plot(x = 1,
     type = "n",
     xlim=c(0,16),
     ylim=c(0,16),
     xlab = "Pred. RNase",
     ylab = "Obs. RNase",
     pch = 16,
     main = expression("lnk"[on]*"-GBM prediction"),
     cex.lab=2, 
     cex.axis=2,
     cex.main =2)
box(lwd=2)



#text(x = c(-17, -14,-17,-14),
#     y = c(-6,-6,-7,-7),
#     labels = c("R2_train = ", round(b[2],2), "R2_test =", round(a[2],2)),cex = 1.5)

# Add training points
points(x = fitted,
       y = as.matrix(trainingset[,1]),
       pch = 16,
       col = 'blue',
       cex = 2)

# Add test points
points(x = pred_koff,
       y = y_test,
       pch = 16,
       col = 'red',
       cex = 2)

# Add ref line
abline(coef = c(0,1), 
       lty = 2,lwd = 4)

legend("bottomright",
       legend = c("training", "test"),
       col = c('blue', 'red'),
       pch = c(16, 16),cex = 2,
       bg = "white")
