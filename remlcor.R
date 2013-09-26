# Modified version of the metaSEM reml() function.  Added an argument "cor.constraints" that accepts
# a p x p matrix of correlations between the between-studies variances.  Values in the diagonal
# must be equal to 1, and values in the off-diagonal must be less than 1.  The function assumes that
# if there is a cor.constraints matrix, there is also an RE.constraints matrix
remlcor <- function (y, v, x, data, RE.constraints = NULL, RE.startvalues = 0.1, 
          RE.lbound = 1e-10, cor.constraints = NULL, intervals.type = c("z", "LB"), model.name = "Variance component with REML", 
          suppressWarnings = TRUE, ...) 
{
  # Define variables
  mf <- match.call()
  if (missing(data)) {
    data <- sys.frame(sys.parent())
  }
  else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  my.y <- mf[[match("y", names(mf))]]
  my.v <- mf[[match("v", names(mf))]]
  
  # Get the y and v
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  
  # Get the number of ys, vs, and xs
  if (is.vector(y)) {
    no.y <- 1
    no.studies <- length(y)
  }
  # If there are no additional columns to the variance-covariance matrix, the number of variances
  # is set to 1
  else {
    no.y <- ncol(y)
    no.studies <- nrow(y)
  }
  if (is.vector(v)) {
    no.v <- 1
  }
  else {
    no.v <- ncol(v)
  }
  # If there are no moderators
  if (missing(x)) {
    no.x <- 0
  }
  else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) {
      no.x <- 1
    }
    else {
      no.x <- ncol(x)
    }
  }
  # Check that there are a valid number of variances and covariances
  if (no.v != no.y * (no.y + 1)/2) 
    stop(paste("The expected no. of columns in v is ", no.y * 
      (no.y + 1)/2, " while the observed no. of columns in v is ", 
               no.v, ".", sep = ""))
  v.labels <- vech(outer(1:no.y, 1:no.y, function(x, y) paste("v", 
                                                              x, "_", y, sep = "")))
  y.labels <- paste("y", 1:no.y, sep = "")
  x.labels <- paste("x", 1:no.x, sep = "")
  if (missing(x)) {
    x.design <- matrix(1, ncol = 1, nrow = no.studies)
    p <- no.y
    no.x <- 0
    input.df <- as.matrix(cbind(y, v))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels))
  }
  else {
    x.design <- cbind(1, x)
    no.x <- ncol(x.design)
    p <- no.x * no.y
    input.df <- as.matrix(cbind(y, v, x))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels, 
                                       x.labels))
  }
  Y <- c(t(y))
  miss.vec <- is.na(Y)
  Y <- Y[!miss.vec]
  Y <- matrix(Y, ncol = 1)
  no.es <- length(Y)
  numStats <- no.es - p
  Y <- as.mxMatrix(Y)
  fn1 <- function(x, no.y) {
    temp <- lapply(x, function(x, k) {
      diag(x = x, nrow = k, ncol = k)
    }, k = no.y)
    do.call(cbind, temp)
  }
  temp <- lapply(split(x.design, 1:nrow(x.design)), fn1, no.y = no.y)
  X <- do.call(rbind, temp)
  X <- X[!miss.vec, , drop = FALSE]
  X <- as.mxMatrix(X)
  if (is.matrix(RE.lbound)) {
    if (!all(dim(RE.lbound) == c(no.y, no.y))) 
      stop("Dimensions of \"RE.lbound\" are incorrect.")
  }
  else {
    lbound <- matrix(NA, nrow = no.y, ncol = no.y)
    diag(lbound) <- RE.lbound
  }
  lbound <- bdiagRep(lbound, no.studies)
  lbound <- lbound[!miss.vec, !miss.vec]
  free <- bdiagRep(matrix(TRUE, nrow = no.y, ncol = no.y), 
                   no.studies)
  free <- free[!miss.vec, !miss.vec]
  free <- as.logical(vech(free))
  # If there are no constraints on the random effects
  if (is.null(RE.constraints)) {
    # Put the start values in proper format
    if (is.matrix(RE.startvalues)) {
      if (!all(dim(RE.startvalues) == c(no.y, no.y))) 
        stop("Dimensions of \"RE.startvalues\" are incorrect.")
    }
    else {
      values <- diag(x = RE.startvalues, nrow = no.y, ncol = no.y)
    }
    values <- bdiagRep(values, no.studies)
    values <- values[!miss.vec, !miss.vec]
    Tau.labels <- vech(outer(1:no.y, 1:no.y, function(x, 
                                                      y) {
      paste("Tau2_", x, "_", y, sep = "")
    }))
    Tau.labels <- bdiagRep(vec2symMat(Tau.labels), no.studies)
    Tau.labels <- Tau.labels[!miss.vec, !miss.vec]
    Tau.labels[Tau.labels == "0"] <- NA
    Tau <- mxMatrix("Symm", ncol = no.es, nrow = no.es, free = free, 
                    labels = vech(Tau.labels), lbound = vech(lbound), 
                    values = vech(values), name = "Tau")
  }
  # If there are constraints on the random effects
  else {
    if (!is.matrix(RE.constraints)) 
      RE.constraints <- as.matrix(RE.constraints)
    if (!all(dim(RE.constraints) == c(no.y, no.y))) 
      stop("Dimensions of \"RE.constraints\" are incorrect.")
    
    # If there are constraints on the correlations between random effects
    if (!is.null(cor.constraints))
    {
      if (!is.matrix(cor.constraints))
        cor.constraints <- as.matrix(cor.constraints)
      
      # Variances estimated according to the specified constraints
      sd <- bdiagRep(RE.constraints, no.studies)
      sd <- sd[!miss.vec, !miss.vec]
      sd <- as.mxMatrix(sd, lbound = lbound, name = "sd")
      # Correlation matrix
      cor <- bdiagRep(cor, no.studies)
      cor <- cor[!miss.vec, !miss.vec]
      cor <- mxMatrix(type = "Full", values = cor, free = F, name = "cor")
      # Quadratic product of sd and cor
      Tau <- mxAlgebra(sqrt(sd) %&% cor, name = "Tau")
    }
    else
    {
      Tau <- bdiagRep(RE.constraints, no.studies)
      Tau <- Tau[!miss.vec, !miss.vec]
      Tau <- as.mxMatrix(Tau, lbound = lbound, name = "Tau")
    }
  }
  if (no.y == 1) {
    V <- diag(x = v, nrow = no.studies, ncol = no.studies)
  }
  else {
    V <- matrix2bdiag(v)
  }
  V <- V[!miss.vec, !miss.vec]
  V <- as.mxMatrix(V)
  W <- mxAlgebra(solve(V + Tau), name = "W")
  alpha <- mxAlgebra(solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% 
    Y, name = "alpha")
  obj <- mxAlgebra((log(det(V + Tau)) + log(det(t(X) %*% W %*% 
    X)) + t(Y - X %*% alpha) %*% W %*% (Y - X %*% alpha)), 
                   name = "obj")
  if(!is.null(cor.constraints))
  {
    reml.model <- mxModel(model = model.name, X, Y, V, W, Tau, cor, sd,
                          alpha, obj, mxAlgebraObjective("obj"), mxCI("Tau"))
  }
  else
  {
    reml.model <- mxModel(model = model.name, X, Y, V, W, Tau, 
                          alpha, obj, mxAlgebraObjective("obj"), mxCI("Tau"))
  }
  intervals.type <- match.arg(intervals.type)
  switch(intervals.type, z = mx.fit <- tryCatch(mxRun(reml.model, 
                                                      intervals = FALSE, suppressWarnings = suppressWarnings, 
                                                      ...), error = function(e) e), LB = mx.fit <- tryCatch(mxRun(reml.model, 
                                                                                                                  intervals = TRUE, suppressWarnings = suppressWarnings, 
                                                                                                                  ...), error = function(e) e))
  if (inherits(mx.fit, "error")) {
    cat("Error in running the mxModel:\n")
    warning(print(mx.fit))
  }
  mx.fit@runstate$objectives[[1]]@numObs <- no.studies
  mx.fit@runstate$objectives[[1]]@numStats <- numStats
  out <- list(call = mf, data = input.df, no.y = no.y, no.x = no.x, 
              miss.vec = miss.vec, mx.fit = mx.fit, intervals.type = intervals.type)
  class(out) <- "reml"
  return(out)
}