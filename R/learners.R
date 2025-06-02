univariate_library <- c(
  "SL.glm",
  "SL.mgcv",
  "SL.earth",
  "SL.bayesglm"
)

SL.library <- list(
  "SL.glm",
  "SL.glmnet",
  "SL.bayesglm",
  "SL.bartMachine",
  "SL.mgcv",
  "SL.earth",
  "SL.ranger",
  "SL.ranger2",
  "SL.xgboost",
  "SL.xgboost3",
  c("SL.glm", "screen.glmnet"),
  c("SL.xgboost", "screen.glmnet"),
  c("SL.ranger", "screen.glmnet"),
  c("SL.earth", "screen.glmnet")
)


### custom learners

## GAM

SL.mgcv <- function(
  Y,
  X,
  newX,
  family,
  obsWeights,
  deg.gam = 10,
  cts.num = 5,
  ...
) {
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))

  if (sum(!cts.x) == 0) {
    gam.model <- as.formula(paste(
      "Y~",
      paste(
        paste(
          "s(",
          colnames(X[, cts.x, drop = FALSE]),
          ")",
          sep = ""
        ),
        collapse = "+"
      )
    ))
  } else {
    gam.model <- as.formula(paste(
      "Y~",
      paste(
        paste(
          "s(",
          colnames(X[, cts.x, drop = FALSE]),
          ")",
          sep = ""
        ),
        collapse = "+"
      ),
      "+",
      paste(
        colnames(X[,
          !cts.x,
          drop = FALSE
        ]),
        collapse = "+"
      )
    ))
  }

  fit.gam <- mgcv::gam(
    gam.model,
    data = X,
    family = family,
    weights = obsWeights,
    method = "REML"
  )
  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.mgcv")
  return(out)
}

predict.SL.mgcv <- function(object, newdata, ...) {
  pred <- mgcv::predict.gam(
    object = object$object,
    newdata = newdata,
    type = "response"
  )
}

## Bayesglm

SL.bayesglm <- function(Y, X, newX, family, obsWeights, ...) {
  fit.glm <- arm::bayesglm(
    Y ~ .,
    data = X,
    family = family,
    weights = obsWeights,
    prior.scale = 1,
    prior.df = Inf,
    prior.scale.for.intercept = 3,
    prior.df.for.intercept = Inf,
    scaled = FALSE
  )
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.bayesglm")
  return(out)
}

## xgboost

SL.xgboost3 <- function(..., max_depth = 3, minobspernode = 3) {
  SL.xgboost(..., max_depth = max_depth, minobspernode = minobspernode)
}

## random forest

SL.ranger2 <- function(..., mtry = 3, min.node.size = 5) {
  SL.ranger(..., mtry = mtry, min.node.size = min.node.size)
}
