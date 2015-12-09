c

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
    return(as.matrix(L, ncol=1))
  else 
    return(as.matrix(L))
}

LogicalAnd <- function(a)
{
  res <- a[1]
  for (i in 2:length(a))
  {
    res <- res & a[i]
  }
  return(res)
}
