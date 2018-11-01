explain <- function(x, 
                    y, 
                    data, 
                    method, 
                    na.action = na.omit, 
                    fallen.leaves = TRUE, 
                    split.border.col = 1, 
                    box.palette = "auto", 
                    n.trees=10000, 
                    matching.method = "nearest", 
                    trControl=NULL, 
                    split.fun=NULL, ...) {
  
  # Arguments:
  # x = dataframe or matrix
  # y = outcome vector
  # method = propensity score calculaiton method - one of "glm", "rpart", "rf", "cbps", "svm", "twang"
  # na.action = both for generating propensity scores and matchit
  # fallen.leaves = decision tree aesthetic 
  # split.border.col = color of the split box borders
  # box.palette = decision tree color scheme (default = greens)
  # n.trees = number of gbm iterations passed on to gbm (TWANG)
  # matching.method = specifies a matching method in MatchIt
  # trControl = defines trControl() function for caret's train()
  
  library(caret)
  library(MatchIt)
  library(rpart)
  library(partykit)
  
  # combining names and creating dataframe
  dat.new <- as.data.frame(cbind(y, x))
  names(dat.new)[1] <- "y.var"
  
  if (is.null(trControl)) {
    trControl = trainControl(method = "repeatedcv", repeats = 3)
  }
  
  if (method == "glm") {
    # logistic regression for propensity scores
    lr1 <- glm(y.var ~ ., data = dat.new, family = "binomial", na.action = na.action)
    # matching
    m.out <- matchit(y.var ~ fitted(lr1, type = "prob"), data = dat.new, na.action = na.action, method = matching.method, 
                     replace = TRUE)
    # indicator for matched/unmatched
    newy <- ifelse(m.out$weight == 0, 0, 1)
    newy <- factor(newy, levels = c(0, 1), labels = c("out", "in"))
    dat.new$newy <- newy
    # grow decision tree for determining matched/unmatched
    tree1 <- train(x = x, y = newy, method = "rpart", trControl = trControl)
    # plot tree
    library(rpart.plot)
    rpart.plot(tree1$finalModel, fallen.leaves = fallen.leaves, split.border.col = split.border.col, box.palette = box.palette, split.fun = split.fun)
  } 
  
  else if (method == "rpart") {
    y2 <- ifelse(dat.new$y.var==0, "no", "yes")
    ps.tree1 <- train(x = x, y = y2, method = "rpart", trControl = trControl)
    ps1 <- predict(ps.tree1$finalModel)[, 2]
    dat.new$ps1 <- ps1
    m.out2 <- matchit(y.var ~ ps1, data = dat.new, na.action = na.action, method = matching.method, replace = TRUE)
    newy2 <- ifelse(m.out2$weight == 0, 0, 1)
    dat.new$newy2 <- as.factor(newy2)
    tree2 <- train(x = x, y = as.factor(newy2), method = "rpart", trControl = trControl)
    rpart.plot(tree2$finalModel, fallen.leaves = fallen.leaves, split.border.col = split.border.col, box.palette = box.palette, split.fun = split.fun)
  } 
  
  else if (method == "rf") {
    library(randomForest)
    y2 <- factor(dat.new$y.var, levels=c(0, 1), labels=c("no", "yes"))
    ps.tree2 <- randomForest(x = x, y = y2, trControl = trControl, data=dat.new)
    ps2 <- ps.tree2$votes[,2]
    dat.new$ps2 <- ps2
    m.out3 <- matchit(y.var ~ ps2,data=dat.new,na.action=na.action, method = matching.method, replace=TRUE)
    newy3 <- ifelse(m.out3$weight==0, 0, 1)
    dat.new$newy3 <- factor(newy3, levels=c(0,1), labels=c("out", "in"))
    tree3 <- train(x = x, y = dat.new$newy3, method = "rpart", trControl = trControl)
    rpart.plot(tree3$finalModel, fallen.leaves=fallen.leaves, split.border.col=split.border.col, box.palette = box.palette, split.fun = split.fun)
  } 
  
  else if (method == "cbps") {
    library(CBPS)
    fit <- CBPS(y.var ~ ., data = dat.new, ATT = FALSE)
    ps3 <- fitted(fit)
    dat.new$ps3 <- ps3
    m.out4 <- matchit(y.var~ps3,data=dat.new,na.action=na.action,method = matching.method, replace=TRUE)
    newy4 <- ifelse(m.out4$weight==0, 0, 1)
    dat.new$newy4 <- factor(newy4, levels=c(0,1), labels=c("out", "in"))
    tree4 <- train(x = x, y = as.factor(dat.new$newy4), method = "rpart", trControl = trControl)
    rpart.plot(tree4$finalModel, fallen.leaves=TRUE, split.border.col=split.border.col, box.palette = box.palette, split.fun=split.fun)
  } 
  
  else if (method=="svm") {
    dat.new$y.var <- as.factor(dat.new$y.var)
    library(e1071)
    attach(dat.new)
    x2 <- data.matrix(x)
    svm1 <- svm(x2,dat.new$y.var, probability=T)
    ps4 <- predict(svm1, x2, probability=T)
    dat.new$ps4 <- attr(ps4,"probabilities")[,1]
    m.out5 <- matchit(y.var~ps4,data=dat.new,na.action=na.action,method=matching.method,replace=TRUE)
    newy5 <- ifelse(m.out5$weight==0, 0, 1)
    dat.new$newy5 <- factor(newy5, levels=c(0,1), labels=c("out", "in"))
    tree5 <- train(x = x, y = dat.new$newy5, method = "rpart", trControl = trControl)
    rpart.plot(tree5$finalModel, fallen.leaves=TRUE, split.border.col=1, box.palette = "Grays", split.fun=split.fun)
  } 
  
  else if (method=="twang") {
    library(twang)
    dat.new2 <- dat.new
    ps.hn<- ps(y.var ~ .,
               dat.new2,
               n.trees=n.trees,
               interaction.depth=2,
               shrinkage=0.01,
               perm.test.iters=0,
               stop.method=c("es.mean","ks.max"),
               estimand = "ATE",
               verbose=FALSE)
    ps5 <- ps.hn[["ps"]][["es.mean.ATE"]]
    m.out6 <- matchit(y.var ~ ps5, data=dat.new2, na.action=na.action, method=matching.method, replace=T)
    newy6 <- ifelse(m.out6$weight==0, 0, 1)
    dat.new2$newy6 <- factor(newy6, levels=c(0,1), labels=c("out", "in"))
    tree6 <- train(x = x, y = dat.new2$newy6, method = "rpart", trControl = trControl)
    rpart.plot(tree6$finalModel, fallen.leaves=fallen.leaves, split.border.col=split.border.col, box.palette = box.palette, split.fun=split.fun)
  }
}
