# Prediction
# Linear regression
ukb = ukb[,c(1,5)]
true = merge(true, ukb, "IID")
true = na.omit(true)
predicted = merge(predicted, ukb, "IID")
x_train = true[,c(3:19,21,22,48)]
# x_train = true[1:5000,c(3:19,21,22,48)]
y_train = true[,20]
# y_train = true[1:5000,20]
x = cbind(x_train, y_train)


x_test = predicted[,c(3:19,21,22,36)]
# x_test = true[20000:25000,c(3:19,21,22,48)]

linear = lm(y_train ~ ., data = x)
summary(linear)
test = true[20000:25000,c(3:19,20,21,22,48)]
test$linear = predict(linear,x_test)
cor.test(test$avMSE, test$linear)
plot(density(predicted$predicted_avMSE), col = "blue"); lines(density(true$avMSE)); lines(density(predicteds), col = "red")

comp = cbind(predicted, predicteds)
cor.test(comp$predicted_avMSE, comp$predicteds)
summary(lm(predicteds ~ Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Age+UniEdu+Sex+AgeSpexWear,comp))
summary(lm(predicted_avMSE ~ Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Age+UniEdu+Sex+AgeSpexWear,predicted))


# Decission tree
library(rpart)
fit = rpart(y_train ~ ., data = x,method="class")
summary(fit)

test$tree = predict(fit,x_test)


# SVM
library(e1071)
fit = svm(y_train ~ ., data = x)
summary(fit)

test$svm = predict(fit,x_test)
cor.test(test$avMSE, test$svm)
comp = cbind(predicted, predicteds)
lines(density(predicteds), col = "gold3")
cor.test(comp$predicted_avMSE, comp$predicteds)
summary(lm(predicteds ~ Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Age+UniEdu+Sex+AgeSpexWear,comp))


# Naive Bayes
fit = naiveBayes(y_train ~ ., data = x)
summary(fit)

test$Bayes = predict(fit,x_test)
cor.test(test$avMSE, test$Bayes)
comp = cbind(predicted, predicteds)
lines(density(predicteds), col = "violet")
cor.test(comp$predicted_avMSE, comp$predicteds)
summary(lm(predicteds ~ Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Age+UniEdu+Sex+AgeSpexWear,comp))


# KNN method
library(class)
fit = knn(y_train ~ ., data = x, k=5)
summary(fit)

predicteds = predict(fit,x_test)


# K Means
library(cluster)
fit = kmeans(x, 3) # 5 cluster solution
summary(fit)
x_test$kmeans = predict(fit,x_test)

# Random Forest
library(randomForest)
fit = randomForest(y_train ~ ., x, ntree=500)
summary(fit)
test$rf = predict(fit,x_test)
cor.test(test$avMSE, test$rf)

# Dimensionality Reduction
library(stats)
pca = princomp(x_train, cor = TRUE)
train_reduced  = predict(pca,x_train)
test_reduced  = predict(pca,x_test)


# Gradient Boosting
library(caret)
fitControl = trainControl( method = "repeatedcv", number = 4, repeats = 4)
fit = train(y_test ~ ., data = x, method = "gbm", trControl = fitControl,verbose = FALSE)
predicted = predict(fit,x_test,type= "prob")[,2]


# XGBoost
require(caret)
TrainControl = trainControl( method = "repeatedcv", number = 10, repeats = 4)
model = train(y ~ ., data = x, method = "xgbLinear", trControl = TrainControl,verbose = FALSE)
# OR 
model = train(y ~ ., data = x, method = "xgbTree", trControl = TrainControl,verbose = FALSE)
predicted = predict(model, x_test)


# LightGBM
library(RLightGBM)
data(example.binary)
#Parameters

num_iterations <- 100
config <- list(objective = "binary",  metric="binary_logloss,auc", learning_rate = 0.1, num_leaves = 63, tree_learner = "serial", feature_fraction = 0.8, bagging_freq = 5, bagging_fraction = 0.8, min_data_in_leaf = 50, min_sum_hessian_in_leaf = 5.0)

#Create data handle and booster
handle.data <- lgbm.data.create(x)

lgbm.data.setField(handle.data, "label", y)

handle.booster <- lgbm.booster.create(handle.data, lapply(config, as.character))

#Train for num_iterations iterations and eval every 5 steps

lgbm.booster.train(handle.booster, num_iterations, 5)

#Predict
pred <- lgbm.booster.predict(handle.booster, x.test)

#Test accuracy
sum(y.test == (y.pred > 0.5)) / length(y.test)

#Save model (can be loaded again via lgbm.booster.load(filename))
#lgbm.booster.save(handle.booster, filename = "/tmp/model.txt")
#If you're familiar with the Caret package in R, this is another way of implementing the LightGBM.

require(caret)
require(RLightGBM)
data(iris)

model <-caretModel.LGBM()

fit <- train(Species ~ ., data = iris, method=model, verbosity = 0)
print(fit)
y.pred <- predict(fit, iris[,1:4])

library(Matrix)
model.sparse <- caretModel.LGBM.sparse()

#Generate a sparse matrix
mat <- Matrix(as.matrix(iris[,1:4]), sparse = T)
fit <- train(data.frame(idx = 1:nrow(iris)), iris$Species, method = model.sparse, matrix = mat, verbosity = 0)
print(fit)


# CatBoost
set.seed(1)

require(titanic)

require(caret)

require(catboost)

tt <- titanic::titanic_train[complete.cases(titanic::titanic_train),]

data <- as.data.frame(as.matrix(tt), stringsAsFactors = TRUE)

drop_columns = c("PassengerId", "Survived", "Name", "Ticket", "Cabin")

x <- data[,!(names(data) %in% drop_columns)]y <- data[,c("Survived")]

fit_control <- trainControl(method = "cv", number = 4,classProbs = TRUE)

grid <- expand.grid(depth = c(4, 6, 8),learning_rate = 0.1,iterations = 100, l2_leaf_reg = 1e-3,            rsm = 0.95, border_count = 64)

report <- train(x, as.factor(make.names(y)),method = catboost.caret,verbose = TRUE, preProc = NULL,tuneGrid = grid, trControl = fit_control)

print(report)

importance <- varImp(report, scale = FALSE)

print(importance)
