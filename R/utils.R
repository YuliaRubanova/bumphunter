greaterOrEqual <- function(x,y) {
  precision <- sqrt(.Machine$double.eps)
  (x >= y) | (abs(x-y) <= precision)
}

closeSockets <- function() {
    allCon <- showConnections()
    socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
    sapply(socketCon, function(ii) close.connection(getConnection(ii)) )
}

foreachCleanup <- function() {
    if (exists(".revoDoParCluster", where=doParallel:::.options)) {
        if(!is.null(doParallel:::.options$.revoDoParCluster))
            stopCluster(doParallel:::.options$.revoDoParCluster)
        remove(".revoDoParCluster", envir=doParallel:::.options)
    }
}

getLengthMatrixOrVector <- function(x)
{
  if (is.null(dim(x)))
  {
    # This is probably a vector
    res <- length(x)
  } else {
    # This is probably a matrix
    res <- nrow(x)
  }
  return(res)
}

toVerticalMatrix <- function(L)
{
  if (is.vector(L))
    return(matrix(L, ncol=1))
  else 
    return(as.matrix(L))
}

toHorizontalMatrix <- function(L)
{
  if (is.vector(L))
    return(matrix(L, nrow=1))
  else 
    return(as.matrix(L))
}

LogicalAnd <- function(a)
{
  res <- TRUE
  for (i in 1:ncol(a))
  {
    res <- res & a[,i]
  }
  return(res)
}
