#########################################################################################
##
## R ROUTINES FOR HOMWEWORK - STATISTICAL MODELLING AND INFERENCE.

#########################################################################################


lasso.bic.glm <- function(y, x, extended=FALSE, family='gaussian', alpha=1) {
  #Select model in LASSO path with best BIC (using LASSO regression estimates)
  #Input
  # - y: vector with response variable
  # - x: design matrix
  # - extended: whether to use EBIC (Chen and Chen 2008) instead of BIC
  # - family: exponential family distribution to fit model with (for logit use binomial)
  # - alpha: to select lasso or ridge regression
  #Output: list with the following elements
  # - coef: LASSO-estimated regression coefficient with lambda set via BIC
  # - ypred: predicted y
  # - lambda.opt: optimal value of lambda
  # - lambda: data.frame with bic and number of selected variables for each value of lambda
  require(glmnet)
  fit <- glmnet(x=x,y=y,family=family,alpha=alpha)
  pred <- cbind(1,x) %*% rbind(fit$a0,fit$beta)
  n <- length(y)
  p <- colSums(fit$beta!=0) + 1
  if (!extended){
    bic <- deviance(fit) + log(n)*p 
  } else {
    bic <- deviance(fit) + log(n)*p + 2*log(choose(ncol(x),p))
  }
  sel <- which.min(bic)
  beta <- c(fit$a0[sel],fit$beta[,sel]); names(beta)[1]= 'Intercept'
  ypred <- pred[,sel]
  ans <- list(coef=beta,ypred=ypred,lambda.opt=fit$lambda[sel],lambda=data.frame(lambda=fit$lambda,bic=bic,nvars=p))
  return(ans)
}

kfoldCV.glm <- function(y, x, K=10, seed, family = "gaussian", predict_type = "link") {
  ## Perform K-fold cross-validation for least-squares regression estimate
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - family: generalises to exponential family regressions
  ## - predict_type: type of prediction. For logit use 'response' to get probabilities.
  ## - seed: random number generator seed (optional)
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  if (ncol(x)>0) {
    for (k in 1:K) {
      sel <- subset==k
      fit <- glm(y[!sel] ~ ., data=x[!sel,,drop=FALSE], family = family)
      pred[sel] <- predict.glm(fit, newdata=x[sel,,drop=FALSE], type = predict_type)
    }
  } else {
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}

kfoldCV.lasso.glm <- function(y, x, K=10, seed, criterion='cv', family = 'gaussian', predict_type = "link") {
  ## Perform K-fold cross-validation for LASSO regression estimate (lambda set either via cross-val or BIC or EBIC)
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - seed: random number generator seed (optional)
  ## - family: generalises to exponential family regressions
  ## - predict_type: type of prediction. For logit use 'response' to get probabilities (see ?predict.glmnet for options).
  ## - criterion: the criterion to select the penalization parameter, either cross-val or BIC or EBIC
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  require(glmnet)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  cat("Starting cross-validation")
  if (ncol(x)>0) {  #if there are some covariates
    for (k in 1:K) {
      sel <- subset==k
      if (criterion=='cv') {
        fit <- cv.glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], family = family, nfolds=10)
        pred[sel] <- predict(fit, newx= x[sel,,drop=FALSE], type = predict_type) 
      } else if (criterion=='bic'){
        fit <- lasso.bic.glm(y=y[!sel],x=x[!sel,,drop=FALSE], family = family, extended = FALSE)
        model <- glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], lambda = fit$lambda.opt, family=family)
        pred[sel] <- predict(model, newx= x[sel,,drop=FALSE], type = predict_type)
      } else if (criterion=='ebic'){
        fit <- lasso.bic.glm(y=y[!sel],x=x[!sel,,drop=FALSE], family = family, extended = TRUE)
        model <- glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], lambda = fit$lambda.opt, family=family)
        pred[sel] <- predict(model, newx= x[sel,,drop=FALSE], type = predict_type)
      } else { stop("method.lambda not implemented") }
      cat(".")
    }
  } else { #if there are no covariates, just use the intercept
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  cat("\n")
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}



