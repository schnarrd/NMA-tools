# Modified version of the metaSEM meta() function.  Added an argument "cor.constraints" that accepts
# a p x p matrix of correlations between the between-studies variances.  Values in the diagonal
# must be equal to 1, and values in the off-diagonal must be less than 1.  The function assumes that
# if there is a cor.constraints matrix, there is also an RE.constraints matrix
metacor <- function (y, v, x, data, intercept.constraints = NULL, coef.constraints = NULL, 
                     RE.constraints = NULL, RE.startvalues = 0.1, RE.lbound = 1e-10, cor.constraints = NULL,
                     intervals.type = c("z", "LB"), I2 = "I2q", R2 = TRUE, model.name = "Meta analysis with ML", 
                     suppressWarnings = TRUE, ...) 
{
  # Define variables
  mf <- match.call()
  if (missing(data)) {
    data <- sys.frame(sys.parent())
  }
  else {
    if (!is.data.frame(data)) 
      data <- data.frame(data)
  }
  my.y <- mf[[match("y", names(mf))]]
  my.v <- mf[[match("v", names(mf))]]
  
  # Get the y and v
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  
  # Get the number of ys, vs, and xs
  if (is.vector(y)) 
    no.y <- 1
  else no.y <- ncol(y)
  # If there are no additional columns to the variance-covariance matrix, the number of variances
  # is set to 1
  if (is.vector(v)) 
    no.v <- 1
  else no.v <- ncol(v)
  # If there are no moderators
  if (missing(x)) 
    no.x <- 0
  else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) 
      no.x <- 1
    else no.x <- ncol(x)
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
  if (no.y == 1) {
    y[is.na(v)] <- NA
  }
  else {
    index <- matrix(0, nrow = no.y, ncol = no.y)
    index[lower.tri(index, diag = TRUE)] <- seq(1, no.y * 
      (no.y + 1)/2)
    index <- diag(index)
    y[is.na(v[, index])] <- NA
  }
  
  # Prevents an error
  v[is.na(v)] <- 1e+10
  # If there are no moderators, the input data consist of the ys and the variance-covariance matrix
  if (no.x == 0) {
    input.df <- as.matrix(cbind(y, v))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels))
    miss.x <- rep(FALSE, nrow(input.df))
  }
  # If there are moderators, add the xs to the input data
  else {
    input.df <- as.matrix(cbind(y, v, x))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels, 
                                       x.labels))
    if (no.x == 1) 
      miss.x <- is.na(x)
    else miss.x <- apply(is.na(x), 1, any)
  }
  # Remove missing values
  my.df <- input.df[!miss.x, ]
  # Check for intercept constraints.  If there are none, estimate
  # the intercepts freely
  if (is.null(intercept.constraints)) {
    Inter <- matrix(paste("0*Intercept", 1:no.y, sep = ""), 
                    nrow = 1, ncol = no.y)
  }
  # If there are intercept constraints, check that the intercept
  # constraints are a valid matrix
  else {
    if (!is.matrix(intercept.constraints)) 
      intercept.constraints <- t(as.matrix(intercept.constraints))
    if (!all(dim(intercept.constraints) == c(1, no.y))) 
      stop("Dimensions of \"intercept.constraints\" are incorrect.")
    Inter <- intercept.constraints
  }
  # Put the intercept matrix in the proper format for OpenMx using as.mxMatrix
  Inter <- as.mxMatrix(t(Inter), name = "Inter")
  # If there are no moderators, X is a 1x1 unit matrix, Beta1 is also a 1x1 unit matrix,
  # and Beta is a zero matrix
  if (no.x == 0) {
    X <- mxMatrix("Unit", nrow = 1, ncol = 1, name = "X")
    Beta1 <- mxAlgebra(Inter, name = "Beta1")
    Beta <- mxMatrix("Zero", nrow = 1, ncol = 1, name = "Beta")
  }
  # If there are moderators
  else {
    # Check for constraints on the coefficient values
    if (is.null(coef.constraints)) {
      # If there are none, then allow the coefficients to be freely estimated, given the
      # outcomes and moderators
      yVar <- paste("y", seq(1, no.y), sep = "", collapse = "+")
      xVar <- paste("x", seq(1, no.x), sep = "", collapse = "+")
      startValues <- tryCatch(eval(parse(text = paste("t(coefficients(lm(cbind(", 
                                                      yVar, ")~", xVar, ", data=data.frame(my.df))))", 
                                                      sep = ""))))
      # If the startValues object is an error, plug in default start values
      if (inherits(startValues, "error")) 
        startValues <- matrix(0, nrow = no.y, ncol = (no.x + 
          1))
      
      A.labels <- outer(1:no.y, 1:no.x, function(y, x) paste("*Slope", 
                                                             y, "_", x, sep = ""))
      Beta <- matrix(paste(startValues[, -1], A.labels, 
                           sep = ""), nrow = no.y, ncol = no.x)
    }
    # If there are no constraints on the coefficient values
    else {
      # Make the constraints a matrix and check that the dimensions of the matrix are valid
      if (!is.matrix(coef.constraints)) 
        coef.constraints <- as.matrix(coef.constraints)
      coef.dim <- dim(coef.constraints)
      if (!coef.dim[1] == no.y | !(coef.dim[2] %in% c(no.x, 
                                                      no.x + no.y))) 
        stop("Dimensions of \"coef.constraints\" are incorrect.")
      Beta <- coef.constraints
    }
    # Put the coefficient constraints in a valid OpenMx format
    Beta <- as.mxMatrix(Beta)
    Beta1 <- mxAlgebra(cbind(Inter, Beta), name = "Beta1")
    X <- mxMatrix("Full", nrow = 1, ncol = (1 + no.x), free = FALSE, 
                  values = c(1, rep(NA, no.x)), labels = c(NA, paste("data.x", 
                                                                     1:no.x, sep = "")), name = "X")
  }
  # Create an mxAlgebra object to hold the transformation matrix required to estimate
  # the moderator coefficient values
  expMean <- mxAlgebra(X %*% t(Beta1), name = "expMean")
  # Format the matrix containing the lower bounds on the random effects in an appropriate manner
  if (is.matrix(RE.lbound)) {
    if (!all(dim(RE.lbound) == c(no.y, no.y))) 
      warning("Dimensions of \"RE.lbound\" are incorrect.")
    lbound <- RE.lbound
  }
  else {
    lbound <- matrix(NA, nrow = no.y, ncol = no.y)
    diag(lbound) <- RE.lbound
  }
  # If there are no constraints on the random effects
  if (is.null(RE.constraints)) {
    # Format the start values appropriately
    if (is.matrix(RE.startvalues)) {
      if (!all(dim(RE.startvalues) == c(no.y, no.y))) 
        warning("Dimensions of \"RE.startvalues\" are incorrect.")
      values <- vech(RE.startvalues)
    }
    else {
      values <- vech(diag(x = RE.startvalues, nrow = no.y, 
                          ncol = no.y))
    }
    # Tau is a transformation matrix required to treat the meta-analysis as an SEM
    Tau.labels <- vech(outer(1:no.y, 1:no.y, function(x, 
                                                      y) {
      paste("Tau2_", x, "_", y, sep = "")
    }))
    Tau <- mxMatrix("Symm", ncol = no.y, nrow = no.y, free = TRUE, 
                    labels = Tau.labels, lbound = vech(lbound), values = values, 
                    name = "Tau")
  }
  # Otherwise, there are constraints on the random effects
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
      
      # Each variance is estimated as a free parameter
      sd <- as.mxMatrix(RE.constraints, name = "sd", lbound = c(lbound))
      # Correlation matrix
      cor <- mxMatrix(type = "Full", values = cor, free = F, name = "cor")
      # Quadratic product of sd and cor
      Tau <- mxAlgebra(sqrt(sd) %&% cor, name = "Tau")
    }
    # If there are, get the random effect constraints into OpenMx format
    else
    {
      Tau <- as.mxMatrix(RE.constraints, lbound = c(lbound), 
                         name = "Tau")      
    }    
  }
  V <- mxMatrix("Symm", ncol = no.y, nrow = no.y, free = FALSE, 
                labels = paste("data.", v.labels, sep = ""), name = "V")
  # The expected covariance matrix
  expCov <- mxAlgebra(V + Tau, name = "expCov")
  # Initialize the comparison model fit object
  mx0.fit <- NA
  # If there are no moderators
  if (no.x == 0) {
    # Grab the I2 values
    I2 <- match.arg(I2, c("I2q", "I2hm", "I2am"), several.ok = TRUE)
    # Grab the variances
    v_het <- input.df[, paste("v", 1:no.y, "_", 1:no.y, sep = ""), 
                      drop = FALSE]
    # Calculate the inverses abd the inverse squares
    sum.w <- apply(v_het, 2, function(x) sum(1/x))
    sum.w2 <- apply(v_het, 2, function(x) sum(1/x^2))
    # Calculate the number of studies by counting all the columns for which the 
    # heterogeneity is a number below 1 trillion (which is the heterogeneity value
    # assigned to studies with missing values)
    no.studies <- apply(v_het, 2, function(x) sum(x < 1e+09))
    # Calculate the transformation matrices required to treat the meta-analysis as an SEM
    qV <- matrix((no.studies - 1) * sum.w/(sum.w^2 - sum.w2), 
                 nrow = 1)
    hmV <- matrix(no.studies/sum.w, nrow = 1)
    amV <- apply(v_het, 2, function(x) mean(x[x < 1e+09]))
    amV <- matrix(amV, nrow = 1)
    V_het <- rbind(qV, hmV, amV)
    V_het <- matrix(t(V_het[c("I2q", "I2hm", "I2am") %in% 
      I2, ]), ncol = 1)
    V_het <- as.mxMatrix(V_het)
    # Required to estimate the intercept
    One <- mxMatrix("Unit", nrow = length(I2), ncol = 1, 
                    name = "One")
    Tau_het <- mxAlgebra(One %x% diag2vec(Tau), name = "Tau_het")
    # Fit index
    I2_values <- mxAlgebra(Tau_het/(Tau_het + V_het), name = "I2_values")

    # If there were constraints on the correlations, add the cor and sd objects to the mxModel object
    if(!is.null(cor.constraints))
    {
      mx.model <- mxModel(model = model.name, mxData(observed = my.df, 
                                                     type = "raw"), mxFIMLObjective(covariance = "expCov", 
                                                                                    means = "expMean", dimnames = y.labels), Inter, Beta, 
                          Beta1, expMean, X, expCov, Tau, cor, sd, V, One, V_het, Tau_het, 
                          I2_values, mxCI(c("Tau", "Inter", "I2_values")))
    }
    # Otherwise, fit the model as normal
    else
    {
      mx.model <- mxModel(model = model.name, mxData(observed = my.df, 
                                                     type = "raw"), mxFIMLObjective(covariance = "expCov", 
                                                                                    means = "expMean", dimnames = y.labels), Inter, Beta, 
                          Beta1, expMean, X, expCov, Tau, V, One, V_het, Tau_het, 
                          I2_values, mxCI(c("Tau", "Inter", "I2_values")))      
    }
  }
  # If there are no moderators
  else {
    # If there are constraints on the correlations, add the cor and sd objects to the mxModel
    if(!is.null(cor.constraints))
    {
      mx.model <- mxModel(model = model.name, mxData(observed = my.df, 
                                                   type = "raw"), mxFIMLObjective(covariance = "expCov", 
                                                                                  means = "expMean", dimnames = y.labels), Inter, Beta, 
                        Beta1, expMean, X, expCov, Tau, cor, sd, V, mxCI(c("Tau", 
                                                                  "Inter", "Beta")))
    }
    # Otherwise, fit the model as usual
    else
    {
      mx.model <- mxModel(model = model.name, mxData(observed = my.df, 
                                                     type = "raw"), mxFIMLObjective(covariance = "expCov", 
                                                                                    means = "expMean", dimnames = y.labels), Inter, Beta, 
                          Beta1, expMean, X, expCov, Tau, V, mxCI(c("Tau", 
                                                                             "Inter", "Beta")))
    }
    # Calculate the base model if R2 was requested
    if (R2) 
      mx0.fit <- tryCatch(meta(y = y, v = v, data = my.df, 
                               model.name = "No predictor", suppressWarnings = TRUE, 
                               silent = TRUE), error = function(e) e)
  }
  # Calculate the requested intervals
  intervals.type <- match.arg(intervals.type)
  switch(intervals.type, z = mx.fit <- tryCatch(mxRun(mx.model, 
                                                      intervals = FALSE, suppressWarnings = suppressWarnings, 
                                                      ...), error = function(e) e), LB = mx.fit <- tryCatch(mxRun(mx.model, 
                                                                                                                  intervals = TRUE, suppressWarnings = suppressWarnings, 
                                                                                                                  ...), error = function(e) e))
  # Throw an error if any errors occurred 
  if (inherits(mx.fit, "error")) {
    cat("Error in running mxModel:\n")
    warning(print(mx.fit))
  }
  
  # Organize and return the output
  out <- list(call = mf, data = input.df, no.y = no.y, no.x = no.x, 
              miss.x = miss.x, mx.model = mx.model, I2 = I2, R2 = R2, 
              mx.fit = mx.fit, mx0.fit = mx0.fit, intervals.type = intervals.type)
  class(out) <- "meta"
  return(out)
}