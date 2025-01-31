##############################################################################
# recode: shift and re-map a numeric variable so its labels start at 1 and
#         go up to the number of unique values.
##############################################################################
recode <- function(x) {
  # shift so smallest value is 1
  x_shifted <- x - min(x) + 1
  
  # tabulate frequencies
  tab <- as.data.frame(table(x_shifted))
  tab <- tab[tab$Freq > 0, ]  # keep only values that occur
  # tab looks like:   x_shifted   Freq
  
  # new labels from 1 to number of unique values
  # 'Var1' is the old label, want to match to row index
  unique_vals <- as.numeric(as.character(tab$x_shifted))
  new_labels  <- seq_len(nrow(tab))  # 1 ... K

  # map old -> new
  old_to_new <- setNames(new_labels, unique_vals)
  # e.g. old_to_new[2] = 1, old_to_new[3] = 2, etc.

  # recode x
  x_recoded <- old_to_new[as.character(x_shifted)]
  return(as.numeric(x_recoded))
}


##############################################################################
# acfcomp2: compute sample autocorrelations for lags 1..na, similarly to
#           the acfcomp2.m in Matlab
##############################################################################
acfcomp2 <- function(tseries, na_lags = 10) {
  n <- length(tseries)
  # default for na_lags might be something like ceiling(10*log10(n)) if needed
  # but we’ll just go with user input
  result <- numeric(na_lags)
  for (j in seq_len(na_lags)) {
    # correlation with lag j
    result[j] <- cor(tseries[1:(n-j+1)], tseries[j:n])
  }
  return(result)
}


##############################################################################
# cluster_se: cluster-robust standard errors by a grouping variable
#   x      = design matrix used in the regression
#   e      = vector of residuals
#   XpXinv = inverse(X'X)
#   group  = integer or factor specifying the cluster for each observation
#   k      = number of parameters (if omitted, use ncol(XpXinv))
#
# This replicates the logic from the Matlab code.
##############################################################################
cluster_se <- function(x, e, XpXinv, group, k = NULL) {
  n <- length(e)
  if (is.null(k)) {
    k <- ncol(XpXinv)
  }

  # Sum (x_i' e_i)(x_i' e_i)' over clusters
  # cluster index goes from 1 .. max(group), or use unique() if factor
  # Make sure 'group' is integer or factor.
  if (is.factor(group)) group <- as.integer(group)
  max_grp <- max(group)

  V <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  for (i in seq_len(max_grp)) {
    I <- which(group == i)
    XiEi <- t(x[I, , drop=FALSE]) %*% e[I]
    V <- V + XiEi %*% t(XiEi)
  }

  # middle = V from the Matlab code
  # cluster variance formula:
  #  ((n-1)/(n-k)) * (max_grp/(max_grp-1)) * XpXinv * V * XpXinv
  factor1 <- ((n - 1)/(n - k)) * (max_grp/(max_grp - 1))
  vcluster <- factor1 * (XpXinv %*% V %*% XpXinv)
  se <- sqrt(diag(vcluster))

  list(se = se, vcluster = vcluster, middle = V)
}


##############################################################################
# findNonCollinear: removes columns that are nearly collinear with existing
#                    ones, up to some tolerance. This is an ad-hoc approach
#                    that tries to replicate the Matlab code's logic.
##############################################################################
findNonCollinear <- function(Z, tol = 1e-10) {
  # Replicates the same logic as Matlab:
  #   for ii = 1:p
  #       use = keep ~= p-ii+1;
  #       e = Z(:,p-ii+1)-Z(:,use)*pinv(Z(:,use)'*Z(:,use))*(Z(:,use)'*Z(:,p-ii+1));
  #       if sum(e.^2)/sum(Z(:,p-ii+1).^2) < tol
  #           keep = setdiff(keep,p-ii+1);
  #       end
  #   end

  p <- ncol(Z)
  keep <- seq_len(p)

  for (ii in seq_len(p)) {
    # The column to check is p-ii+1, going backward from p down to 1
    col_check <- p - ii + 1
    # "use" are the columns we haven't dropped yet, except col_check
    use <- keep[keep != col_check]
    # If we've already dropped everything else, no reason to keep checking
    if (length(use) == 0) break

    z_col <- Z[, col_check]
    z_use <- Z[, use, drop = FALSE]

    # e = z_col - z_use * pinv( z_use' * z_use ) * (z_use' * z_col)
    # We'll replicate Matlab pinv(...) with MASS::ginv(...)
    xtx <- crossprod(z_use)            # z_use' * z_use
    pinv_xtx <- MASS::ginv(xtx)        # pseudo-inverse
    xty <- crossprod(z_use, z_col)     # z_use' * z_col
    e <- z_col - z_use %*% (pinv_xtx %*% xty)

    # ratio = sum(e^2)/sum(z_col^2)
    ratio <- sum(e^2) / sum(z_col^2)
    if (ratio < tol) {
      # remove this col_check from keep
      keep <- setdiff(keep, col_check)
    }
  }

  # sort so we return indices in ascending order, just like the Matlab version
  return(sort(keep))
}


##############################################################################
# LassoShooting2: coordinate descent (“shooting”) algorithm for LASSO,
#                 trying to match structure of the Matlab code as closely
#                 as possible.
#
# Arguments:
#   X        = design matrix
#   y        = outcome vector
#   lambda   = the penalty scaling factor (the code multiplies this further)
#   Ups      = vector of penalization multipliers (the code uses lambda*Ups[j])
#   maxIter  = maximum iterations
#   verbose  = how verbose to be
#   optTol   = tolerance for changes in coefficients
#   ...
##############################################################################
LassoShooting2 <- function(X, y, lambda, Ups,
                           maxIter = 10000,
                           verbose = 2,
                           optTol = 1e-5,
                           zeroThreshold = 1e-4,
                           beta = NULL,
                           XX = NULL,
                           Xy = NULL) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(XX)) {
    XX <- t(X) %*% X
  }
  if (is.null(Xy)) {
    Xy <- t(X) %*% y
  }
  if (is.null(beta)) {
    # Start from a ridge solution or LS with small diagonal augment
    beta <- solve(XX + diag(lambda, p, p), Xy)
  }

  # Precompute for speed
  XX2 <- XX * 2
  Xy2 <- Xy * 2

  if (verbose == 2) {
    cat(sprintf("%10s %10s %15s %15s %15s\n",
                "iter","shoots","L1norm(B)","Change(B)","ObjFn"))
  }

  iter <- 0
  while (iter < maxIter) {
    beta_old <- beta
    for (j in seq_len(p)) {
      # sum(XX2[j, ] * beta) - XX2[j, j]*beta[j] - Xy2[j]
      S0 <- (XX2[j, , drop=FALSE] %*% beta)[1] - XX2[j, j]*beta[j] - Xy2[j]
      # Check sign for thresholding:
      if (S0 > lambda * Ups[j]) {
        beta[j] <- (lambda * Ups[j] - S0) / XX2[j, j]
      } else if (S0 < -lambda * Ups[j]) {
        beta[j] <- (-lambda * Ups[j] - S0) / XX2[j, j]
      } else {
        beta[j] <- 0
      }
    }
    iter <- iter + 1

    if (verbose == 2) {
      # compute objective: sum((X*beta-y)^2) + sum(abs(lambda*Ups*beta))
      objVal <- sum((X %*% beta - y)^2) + sum(abs(lambda * Ups * beta))
      cat(sprintf("%10d %10d %15.2e %15.2e %15.2e\n",
                  iter, iter*p,
                  sum(abs(beta)),
                  sum(abs(beta - beta_old)),
                  objVal))
    }

    # termination
    if (sum(abs(beta - beta_old)) < optTol) {
      break
    }
  }

  if (verbose > 0) {
    cat(sprintf("Number of iterations: %d\nTotal Shoots: %d\n", iter, iter*p))
  }

  return(beta)
}

#--------------------------------------------------------------------
# dummy: replicate the "dummy.m" logic
#        dummy(x) => OUT1: dummy variables for x
#        dummy(x,y) => OUT1: dummy for x, OUT2: dummy for y, OUT3: interaction
#--------------------------------------------------------------------
dummy <- function(x, y=NULL) {
  
  if (is.null(y)) {
    # single-argument version => create dummies for x (except the first category)
    t1 <- table(x)
    # remove zero freq
    t1 <- t1[t1 > 0]
    
    categories <- as.numeric(names(t1))
    # If K categories, then we create K-1 dummy columns
    # (like the Matlab code that skips the very first row in tabulate)
    
    K1 <- length(categories)
    # replicate the "table(i+1,1)" logic => skip the first category for OUT1
    OUT1 <- matrix(0, nrow=length(x), ncol=K1-1)
    for (i in seq_len(K1-1)) {
      catval <- categories[i+1]
      OUT1[,i] <- as.numeric(x == catval)
    }
    OUT2 <- NULL
    OUT3 <- NULL
    return(list(OUT1=OUT1, OUT2=OUT2, OUT3=OUT3))
    
  } else {
    # two-argument version => dummies for x, dummies for y, plus all interactions
    t1 <- table(x)
    t1 <- t1[t1 > 0]
    categories_x <- as.numeric(names(t1))
    K1 <- length(categories_x)
    
    # build OUT1 (dummy for x skipping first cat)
    OUT1 <- matrix(0, nrow=length(x), ncol=K1-1)
    for (i in seq_len(K1-1)) {
      catval <- categories_x[i+1]
      OUT1[,i] <- as.numeric(x == catval)
    }
    
    t2 <- table(y)
    t2 <- t2[t2 > 0]
    categories_y <- as.numeric(names(t2))
    K2 <- length(categories_y)
    
    # build OUT2 (dummy for y skipping first cat)
    OUT2 <- matrix(0, nrow=length(y), ncol=K2-1)
    for (j in seq_len(K2-1)) {
      catval <- categories_y[j+1]
      OUT2[,j] <- as.numeric(y == catval)
    }
    
    # build OUT3 => all interactions of x and y dummies
    # (K1-1)*(K2-1) columns
    # The Matlab code does 3D array out3, then reshapes to 2D
    out3_array <- array(0, dim=c(length(x), (K2-1), (K1-1)))
    for (i in seq_len(K1-1)) {
      for (j in seq_len(K2-1)) {
        out3_array[, j, i] <- as.numeric((x == categories_x[i+1]) & 
                                         (y == categories_y[j+1]))
      }
    }
    out3_mat <- matrix(out3_array, nrow=length(x), ncol=(K1-1)*(K2-1))
    
    return(list(OUT1=OUT1, OUT2=OUT2, OUT3=out3_mat))
  }
}

#--------------------------------------------------------------------
# dummyvar: replicate "dummyvar" from Matlab. 
#   In Matlab, dummyvar(x) returns a full set of dummies for each category
#   If x is numeric, we treat each distinct x as a category.
#--------------------------------------------------------------------
dummyvar <- function(x) {
  # Easiest is model.matrix:
  # But that includes an intercept by default, so do:
  # model.matrix(~ factor(x) - 1)
  m <- model.matrix(~ factor(x) - 1)
  return(m)
}