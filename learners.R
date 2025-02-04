SL.library <- c(
  "SL.glm",
  "SL.firth",
  "SL.glmnet",
  "SL.gam2",
  "SL.earth",
  "SL.ranger2"
)

# custom learners

SL.ranger2 <- function(..., mtry = 2, min.node.size = 5) {
  SL.ranger(..., mtry = mtry, min.node.size = min.node.size)
}

SL.gam2 <-
  function(Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4, ...) {

  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,")", sep=""), collapse = "+"), "+", paste(colnames(X[, !cts.x, drop=FALSE]), collapse = "+")))
  } else {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""), collapse = "+")))
  }
  # fix for when all variables are binomial
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+"), sep = ""))
  }
  fit.gam <- gam::gam(gam.model, data = X, family = family, control = gam::gam.control(maxit = 50, bf.maxit = 50), weights = obsWeights)
  if(packageVersion('gam') >= "1.15") {
    pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response") # updated gam class in version 1.15
  } else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the
    gam package with 'update.packages('gam')'")
  }
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam")
  return(out)
  }


SL.firth <- function(Y, X, newX, ...) {

  if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  fit.firth <- logistf(Y ~ ., data = X)

  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }

  pred <- predict(fit.firth, newdata = newX, type = "response")
  fit <- list(object = fit.firth)
  class(fit) <- "SL.firth"
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.firth <-
  function(object, newdata, ...) {
    # newdata must be a dataframe, not a matrix.
    if (is.matrix(newdata)) {
      newdata = as.data.frame(newdata)
    }
    pred <- predict(object = object$object,
      newdata = newdata, type = "response")
    pred
  }
