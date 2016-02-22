library(doRNG)

computation.tots <- function(tab, V, L, D, workers, nulltabs, current_covariate) {
  LDiff <- cbind(tab$L, abs(tab[,paste0("covariate.diff", current_covariate)]))
  chunksize <- ceiling(nrow(LDiff)/workers)
  subL <- NULL
  tots <- foreach(subL = iter(LDiff, by = "row", chunksize = chunksize),
                  .combine = "cbind", .packages = "bumphunter") %dorng%
  {
    greaterOrEqual <- function(x,y) {
      precision <- sqrt(.Machine$double.eps)
      (x >= y) | (abs(x-y) <= precision)
    }

    BiggerNullBumpIndicesForPermutation <- function(list_for_comparison, value, pairwise=F)
    {
      tmp <- sapply(1:length(value), function(j) {
        if (pairwise & !is.vector(list_for_comparison))
          to_compare <- list_for_comparison[,j]
        else
          to_compare <- list_for_comparison
        greaterOrEqual(to_compare, value[j])
      } )
      return(LogicalAnd(tmp))
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
    
    apply(subL, 1, function(x) {
      res <- sapply(seq(along = D), function(i) {
        sum(greaterOrEqual(L[[i]], x[1]) &
              greaterOrEqual(abs(D[[i]]),
                             x[2]))
      })
      c(mean(res > 0), sum(res))
    })
  }
  attributes(tots)[["rng"]] <- NULL
  rate1 <- tots[1, ]
  pvalues1 <- tots[2, ]/sum(sapply(nulltabs, nrow))
  return(list(rate1 = rate1, pvalues1 = pvalues1))
}





computation.tots2 <- function(tab, A, workers, nulltabs) {
  Avalue <- matrix(tab$area, ncol = 1)
  chunksize <- ceiling(nrow(Avalue)/workers)
  subA <- NULL
  tots2 <- t(foreach(subA = iter(Avalue, by = "row", chunksize = chunksize),
                     .combine = "cbind", .packages = "bumphunter") %dorng%
  {
    greaterOrEqual <- function(x,y) {
      precision <- sqrt(.Machine$double.eps)
      (x >= y) | (abs(x-y) <= precision)
    }
    
    
    toHorizontalMatrix <- function(L)
    {
      if (is.vector(L))
        return(matrix(L, nrow=1))
      else 
        return(as.matrix(L))
    }
    
    toHorizontalMatrix(sapply(subA, function(x) {
      return(sapply(seq(along = A), function(i) {
        sum(greaterOrEqual(A[[i]], x[1]))
      }))
    }))
  })
  
  if (is.vector(tots2)) {
    tots2 <- matrix(tots2, nrow = 1)
  }
  
  rate2 <- rowMeans(tots2 > 0)
  pvalues2 <- rowSums(tots2)/sum(sapply(nulltabs, nrow))
  return(list(rate2 = rate2, pvalues2 = pvalues2))
}







computation.tots.jointly <- function(tabs, V, L, A, D, maxGap, chr, pos, mat, beta, workers,
                                     nulltabs, coef,
                                     verbose = F, 
                                     bumpDirections = NULL, SamplesToDetermineDirection = NULL,
                                     SamplesContraintedByDistribution = NULL, distribution = NULL,
                                     ks.threshold = NULL) {
  if (verbose)
    message("[bumphunterEngine] Joining bumps.")
  
  #joined_tabs <- findIntersection(tabs, maxGap = maxGap)
  joined_tabs <- tabs[order(tabs$chr, tabs$start, tabs$end),]
  
  all_cpg = cbind.data.frame(CHR=unlist(lapply(chr, as.character)), MAPINFO=pos, stringsAsFactors=F)
  cpgs = find_cpg_in_table(joined_tabs, all_cpg, region_names=1:nrow(joined_tabs), returnIndices=T)
  
  
  # Samples specified by SamplesContraintedByDistribution should come from a specified distribution
  if (is.null(SamplesContraintedByDistribution) || is.null(distribution))
  {
    controls_come_from_distribution <- matrix(TRUE, nrow=nrow(cpgs))
  } else {
    if (verbose)
      message("[bumphunterEngine] Checking that controls come from the distribution.")
    ptime1 <- proc.time()
    
    chunksize <- ceiling(length(cpgs$indices)/workers) 

    controls_come_from_distribution.distances <- foreach(subMat = iter(mat[cpgs$indices, SamplesContraintedByDistribution], by = "row", chunksize = chunksize),
                   .combine = "c") %dorng% {
                     apply(subMat, 1, function(x) {suppressWarnings(ks.test(x, distribution )$statistic)})
                   }
    attributes(controls_come_from_distribution.distances)[["rng"]] <- NULL
    
    if (is.null(ks.threshold))
      ks.threshold <- quantile(controls_come_from_distribution.distances)[2]
    
    controls_come_from_distribution <- (controls_come_from_distribution.distances < ks.threshold)
    if (verbose)
      print(proc.time() - ptime1)
  }
  
  cpgs <- cpgs[controls_come_from_distribution, ]
  
  ind <- cpgs$indices
  cluster <- cpgs$name
  Indexes <- vector("list", length=2)
  Indexes[[1]] <- list()
  Indexes[[2]] <- aggregate(1:length(cluster), by=list(cluster), FUN=list)[,2]
  names(Indexes[[2]]) <- unique(cluster)
  
  joined_tabs <- data.frame(chr=sapply(Indexes[[2]],function(Index) chr[ind[Index[1]]]),
                            start=sapply(Indexes[[2]],function(Index) min(pos[ind[Index]])),
                            end=sapply(Indexes[[2]], function(Index) max(pos[ind[Index]])))
  
  cpgs = find_cpg_in_table(joined_tabs, all_cpg, region_names=1:nrow(joined_tabs), returnIndices=T)
  
  ind <- cpgs$indices
  cluster <- cpgs$name
  Indexes <- vector("list", length=2)
  Indexes[[1]] <- list()
  Indexes[[2]] <- aggregate(1:length(cluster), by=list(cluster), FUN=list)[,2]
  names(Indexes[[2]]) <- unique(cluster)
  
  if (verbose)
    message("[bumphunterEngine] Recomputing parameters for joined tabs.")
  joined_tabs_with_values <- list()
  LvalueList <- list()
  AvalueList <- list()
  for (j in 1:length(coef))
  {
    joined_tabs_with_values[[j]] <- regionFinder(x = beta[,j], chr = chr, pos = pos, cluster=cluster,
                                                 cutoff = cutoff, ind = ind, verbose = FALSE, 
                                                 addMeans = T, mat=mat, design=design,
                                                 Indexes=Indexes, clusterInSelectedPositions=T)
    joined_tabs_with_values[[j]] <- joined_tabs_with_values[[j]][order(joined_tabs_with_values[[j]]$chr, joined_tabs_with_values[[j]]$start,joined_tabs_with_values[[j]]$end),]
    # The elements of the lists are Lvalues for each covariate of interest
    LvalueList[[j]] <- cbind.data.frame(L=joined_tabs_with_values[[j]]$L, value=abs(joined_tabs_with_values[[j]]$value))
    AvalueList[[j]] <- cbind.data.frame(area=joined_tabs_with_values[[j]]$area)
  }
  
  Diff <- joined_tabs_with_values[[1]][,grep("covariate.diff*", names(joined_tabs_with_values[[1]]))]
  
  # Control that bumps have one of the directions specified by bumpDirections
  if (is.null(SamplesToDetermineDirection) || is.null(bumpDirections))
  {
    correct_direction <- matrix(TRUE, nrow=nrow(joined_tabs_with_values[[1]]))
  } else {
    correct_direction <- matrix(FALSE, nrow=nrow(joined_tabs_with_values[[1]]))
    for (j in 1:length(bumpDirections))
    {
      correct_direction[which(rowSums(t(bumpDirections[[j]] * t(Diff)) < 0) == 0)] <- TRUE
    }
  }
  
  for (j in 1:length(coef))
  {
    joined_tabs_with_values[[j]] <- joined_tabs_with_values[[j]][correct_direction, ]
    LvalueList[[j]] <- LvalueList[[j]][correct_direction, ]
    AvalueList[[j]] <- cbind.data.frame(area=AvalueList[[j]][correct_direction, ])
  }
  Diff <- Diff[correct_direction,]
  
  Ls <- sapply(LvalueList, function(x) {x$L})
  values <- sapply(LvalueList, function(x) {x$value})
  #Lvalue <- cbind(Ls, values)
  LDiff <- cbind(Ls, abs(Diff))
  Avalue <- sapply(AvalueList, function(x) {x$area})
  
  nullmodel_opposite_directions <- vector("list", length=length(D))
  for (i in 1:length(nullmodel_opposite_directions))
  {
    nullmodel_opposite_directions[[i]] <- matrix(FALSE, nrow=nrow(D[[i]]))
    for (j in 1:length(bumpDirections))
    {
      nullmodel_opposite_directions[[i]][which(rowSums(t(bumpDirections[[j]] * t(D[[i]])) < 0) == 0)] <- TRUE
    }
  }

  if (verbose)
    message("[bumphunterEngine] Computing p-values....")
  ptime1 <- proc.time()
  
  chunksize <- ceiling(nrow(LDiff)/workers)
  subL <- NULL
  tots <- foreach(subL = iter(LDiff, by = "row", chunksize = chunksize),
                  .combine = "cbind", .packages = "bumphunter") %dorng%
  {
    
    greaterOrEqual <- function(x,y) {
      precision <- sqrt(.Machine$double.eps)
      (x >= y) | (abs(x-y) <= precision)
    }
    
    toHorizontalMatrix <- function(L)
    {
      if (is.vector(L))
        return(matrix(L, nrow=1))
      else 
        return(as.matrix(L))
    }
    
    toHorizontalMatrix(apply(subL, 1, function(x) {
      FoundBumpsLs <- x[1:(length(x)/2)]
      FoundBumpsDiff  <- x[((length(x)/2)+1):length(x)]
      
      BiggerNullBumpIndicesForPermutation <- function(list_for_comparison, value, pairwise=F)
      {
        tmp <- sapply(1:length(value), function(j) {
          if (pairwise & !is.vector(list_for_comparison))
            to_compare <- list_for_comparison[,j]
          else
            to_compare <- list_for_comparison
          greaterOrEqual(to_compare, value[j])
        } )
        return(LogicalAnd(tmp))
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

      res <- sapply(seq(along = D), function(i) {
        sum(BiggerNullBumpIndicesForPermutation(L[[i]], FoundBumpsLs) & 
              BiggerNullBumpIndicesForPermutation(abs(D[[i]]), FoundBumpsDiff, pairwise=T) &
              nullmodel_opposite_directions[[i]])
      })
      c(mean(res > 0), sum(res))
    }))
  }
  attributes(tots)[["rng"]] <- NULL
  rate1 <- tots[1, ]
  pvalues1 <- tots[2, ]/sum(sapply(nulltabs, nrow))
  if (verbose)
    print(proc.time() - ptime1)
  
  
  if (verbose)
    message("[bumphunterEngine] Computing area p-values...")
  ptime1 <- proc.time()
  
  chunksize <- ceiling(nrow(Avalue)/workers)
  subA <- NULL
  tots2 <- foreach(subA = iter(Avalue, by = "row", chunksize = chunksize),
                   .combine = "rbind", .packages = "bumphunter") %dorng%
  {
    greaterOrEqual <- function(x,y) 
    {
      precision <- sqrt(.Machine$double.eps)
      (x >= y) | (abs(x-y) <= precision)
    }
    
    toHorizontalMatrix <- function(L)
    {
      if (is.vector(L))
        return(matrix(L, nrow=1))
      else 
        return(as.matrix(L))
    }
    
    BiggerNullBumpIndicesForPermutation <- function(list_for_comparison, value)
    {
      tmp <- sapply(1:length(value), function(j) {
        greaterOrEqual(list_for_comparison, value[j])
      } )
      return(LogicalAnd(tmp))
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
    
    t(toHorizontalMatrix(apply(subA, 1, function(x) {
      controls_come_from_distribution <- TRUE
      #control <- melt(mat[,])
      correct_direction <- TRUE
      
      sapply(seq(along = A), function(i) {
        sum(BiggerNullBumpIndicesForPermutation(A[[i]], x))
      })
    })))
  }
  
  if (is.vector(tots2)) {
    tots2 <- matrix(tots2, nrow = 1)
  }
  
  if (verbose)
    print(proc.time() - ptime1)

  rate2 <- rowMeans(tots2 > 0)
  pvalues2 <- rowSums(tots2)/sum(sapply(nulltabs, nrow))
  
  return(list(joined_tabs_with_values = joined_tabs_with_values, 
              rate1 = rate1, pvalues1 = pvalues1,
              rate2 = rate2, pvalues2 = pvalues2))
}
