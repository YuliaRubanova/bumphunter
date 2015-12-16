.getSegments_old <- function(x, factor, cutoff=quantile(abs(x),0.99), verbose=TRUE){
    Indexes <- split(seq(along=x),factor)
    regionID <- vector("numeric", length(x))
    LAST <- 0
  
    segmentation <- vector("numeric", length(x))
    type <- vector("numeric", length(x))

    for (i in seq(along = Indexes)) {
        if (verbose) if (i%%1000 == 0) cat(".")
        Index <- Indexes[[i]]
        y <- x[Index]
        z <- sign(y) * as.numeric(abs(y) > cutoff)
        w <- cumsum(c(1, diff(z) != 0)) + LAST
        segmentation[Index] <- w
        type[Index] <- z
        LAST <- max(w)
    }
    ##add a vector of the pns
    res <- list(upIndex = split(which(type>0), segmentation[type>0]),
                dnIndex = split(which(type<0), segmentation[type<0]),
                zeroIndex = split(which(type==0), segmentation[type==0]))
    names(res[[1]]) <- NULL
    names(res[[2]]) <- NULL
    names(res[[3]]) <- NULL
    return(res)
}

getSegments <- function(x, f = NULL, cutoff=quantile(abs(x), 0.99), assumeSorted = FALSE, verbose=FALSE){
    if(is.null(f))
        f <- rep(1L, length(x))
    stopifnot(getLengthMatrixOrVector(x) == length(f))
    stopifnot(length(cutoff) <= 2)
    if(is.character(f))
        f <- as.factor(f)
    if(is.numeric(f))
        f <- as.integer(f)
    stopifnot(is.factor(f) || is.integer(f))
    if(length(cutoff) == 1)
        cutoff <- c(-cutoff, cutoff)
    cutoff <- sort(cutoff)

    reordered <- FALSE
    if(!assumeSorted && is.unsorted(f)) {
	od <- order(f)
	x <- x[od,]
	f <- f[od]
	reordered <- TRUE
    }
        
    if(verbose) message("getSegments: segmenting")
    Indexes <- split(seq(getLengthMatrixOrVector(x)), f)
    # !!! Check that this works fine
    direction <- as.integer(greaterOrEqual(x, cutoff[2]))
    direction[x <= cutoff[1]] <- -1L

    ## We now need to convert f into cid
    if(verbose) message("getSegments: splitting")
    segments <- cumsum(c(1, diff(direction) != 0)) +
        cumsum(c(1, diff(f) != 0))
    names(segments) <- NULL

    res <- list(upIndex = split(which(direction>0), segments[direction>0]),
                dnIndex = split(which(direction<0), segments[direction<0]),
                zeroIndex = split(which(direction==0), segments[direction==0]))
    names(res[[1]]) <- NULL
    names(res[[2]]) <- NULL
    names(res[[3]]) <- NULL

    if(reordered) {
        res <- lapply(res, function(sp) lapply(sp, function(xx) od[xx]))
    }
    res
}

clusterMaker <- function(chr, pos, assumeSorted = FALSE, maxGap=300){
    
    nonaIndex <- which(!is.na(chr) & !is.na(pos))
    Indexes <- split(nonaIndex, chr[nonaIndex])
    clusterIDs <- rep(NA, length(chr))
    LAST <- 0
    for(i in seq(along = Indexes)){
        Index <- Indexes[[i]]
        x <- pos[Index]
        if(!assumeSorted){
            Index <- Index[order(x)]
            x <- pos[Index]
        }
        y <- as.numeric(diff(x) > maxGap)
        z <- cumsum(c(1, y))
        clusterIDs[Index] <- z + LAST
        LAST <- max(z) + LAST
    }
    clusterIDs
}

##you can pass cutoff through the ...
regionFinder <- function(x, chr, pos, cluster=NULL, y=x, summary=mean,
                         ind=seq(along=getLengthMatrixOrVector(x)),order=TRUE, oneTable=TRUE,
                         maxGap=300, cutoff=quantile(abs(x), 0.99),
                         assumeSorted = FALSE, verbose = TRUE, 
                         addMeans = F, mat=NULL, design=NULL, null_model_coef=NULL){
    x <- as.matrix(x)
    if(any(is.na(x[ind,]))){
        warning("NAs found and removed. ind changed.")
        ind <- intersect(which(!is.na(x)),ind)
    } 
    if(is.null(cluster))
        cluster <- clusterMaker(chr, pos, maxGap=maxGap, assumeSorted = assumeSorted)
    Indexes <- getSegments(x = x[ind,], f = cluster[ind], cutoff = cutoff,
                           assumeSorted = assumeSorted, verbose = verbose)
    clusterN <- table(cluster)[as.character(cluster)]
    
    res <- vector("list",2)
    for(i in 1:2){
      res[[i]]<-
        data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
                   start=sapply(Indexes[[i]],function(Index) min(pos[ind[Index]])),
                   end=sapply(Indexes[[i]], function(Index) max(pos[ind[Index]])),
                   value=sapply(Indexes[[i]],function(Index)summary(y[ind[Index],])),
                   area=sapply(Indexes[[i]],function(Index)abs(sum(y[ind[Index],]))),
                   cluster=sapply(Indexes[[i]],function(Index)cluster[ind[Index]][1]),
                   indexStart=sapply(Indexes[[i]], function(Index) min(ind[Index])),
                   indexEnd = sapply(Indexes[[i]], function(Index) max(ind[Index])),
                   stringsAsFactors=FALSE)
      
      res[[i]]$L <- res[[i]]$indexEnd - res[[i]]$indexStart+1
      res[[i]]$clusterL <- sapply(Indexes[[i]], function(Index) clusterN[ind[Index]][1])
      
      takeRowMeans <- function(x)
      {
        if (is.null(dim(x)))
        {
          mean(x)
        } else {
          rowMeans(x)
        }
      }
      
      if (addMeans & !is.null(mat))
      {
        if (!is.null(nullmodel_coef))
        {
          res_null <- sapply(Indexes[[i]],function(Index) median(matrix(mat[ind[Index],nullmodel_coef], nrow=1)))
          res_null <- matrix(res_null, ncol=1)
          colnames(res_null)[1] <- "null_model.median" 
          res[[i]] <- cbind(res[[i]], res_null)
        }
        
        if (!is.null(design) && length(Indexes[[i]]) != 0)
        {
          res_design <- apply(design, 2, function(col) {
            sapply(Indexes[[i]],function(Index) median(matrix(mat[ind[Index],which(col > 0)], nrow=1)))
            })
          colnames(res_design) <- paste0("covariate.median", 1:ncol(design))
          res[[i]] <- cbind(res[[i]], res_design)
        }
      }
    }
    names(res) <- c("up","dn")
    if(order & !oneTable){
        if(nrow(res$up)>0) res$up <- res$up[order(-res$up$area),]
        if(nrow(res$dn)>0) res$dn <- res$dn[order(-res$dn$area),]
    }
    if(oneTable){
        res <- rbind(res$up,res$dn)
        if(order & nrow(res)>0) res <- res[order(-res$area),]
    }
    return(res)
}

boundedClusterMaker <- function(chr, pos, assumeSorted = FALSE,
                                maxClusterWidth = 1500, maxGap = 500) {
    nonaIndex <- which(!is.na(chr) & !is.na(pos))
    Indexes <- split(nonaIndex, chr[nonaIndex])
    clusterIDs <- rep(NA, length(chr))
    LAST <- 0
    for (i in seq(along = Indexes)) {
        Index <- Indexes[[i]]
        x <- pos[Index]
        if (!assumeSorted) {
            Index <- Index[order(x)]
            x <- pos[Index]
        }
        y <- as.numeric(diff(x) > maxGap)
        z <- cumsum(c(1, y))
        startPosition <- tapply(x, list(z), min)
        startPositions <- startPosition[as.character(z)]

        offset <- x - startPositions
        addToZ <- offset %/% maxClusterWidth
        addToZ <- cumsum(addToZ * c(0, as.numeric(diff(addToZ) > 0)))
        z <- z + addToZ
        clusterIDs[Index] <- z + LAST
        LAST <- max(z) + LAST
    }
    clusterIDs
}

intersectChromosomalRegions(regions1, regions2)
{
  if (unique(regions1$chr) != 1 || unique(regions2$chr) != 1)
  {
    stop("intersectChromosomalRegions works with regions only from one chromosome. Please separate data per-chromosome.")
  }
  stopifnot(unique(regions1$chr) == unique(regions2$chr))
  
  regions1 <- regions1[,c("chr", "start", "end")]
  regions1 <- regions1[order(regions1$start, regions1$end),]
  regions2 <- regions2[,c("chr", "start", "end")]
  regions2 <- regions2[order(regions2$start, regions2$end),]
  
  regions1.vector <- as.vector(t(regions1[,c("start", "end")]))
  regions2.vector <- as.vector(t(regions2[,c("start", "end")]))
  
  clusterStarts <- findInterval(regions1$start, regions2.vector, rightmost.closed = T)
  # If the position is equal to the end of the region, the findInterval function will put it into even interval,
  # but we want this position to be in odd interval, as we want all odd intervals to be closed because those are the regions with which we intersect.
  not_in_interval <- which(clusterStarts %% 2 == 0)
  should_be_in_interval <- not_in_interval[regions1$start[not_in_interval] %in% regions2$end]
  clusterStarts[should_be_in_interval] <- clusterStarts[should_be_in_interval]-1
  
  clusterEnds <- findInterval(regions1$end, regions2.vector, rightmost.closed = T)
  not_in_interval <- which(clusterEnds %% 2 == 0)
  should_be_in_interval <- not_in_interval[regions1$end[not_in_interval] %in% regions2$end]
  clusterEnds[should_be_in_interval] <- clusterEnds[should_be_in_interval]-1
  
  # The only case when region1 interval does not overlap any regions2 interval is when regions1 intervals turns out to be inside the even interval
  overlap <- !((clusterStarts %% 2 == 0) & (clusterEnds %% 2 == 0) & clusterStarts == clusterEnds)
  # which(overlap) -- rows in regions1 which have at least one overlap with regions2
  # intersecting_regions -- List of intervals in regions2 which have overlap and row in regions1 which correspond to this overlap
  intersecting_regions <- do.call(rbind, lapply(which(overlap), function(x) { cbind(clusterStarts[x]:clusterEnds[x], x) }))
  intersecting_regions <- intersecting_regions[intersecting_regions[,1] %% 2 == 1,]
  # first column of intersecting_rows -- rows in regions2 with have at least one overlap with regions1
  intersecting_rows <- cbind(ceiling(intersecting_regions[,1]/2), intersecting_regions[,2])
  intersecting_rows_regions2 <- regions2[intersecting_rows[,1],]
  intersecting_rows_regions1 <- regions1[intersecting_rows[,2],]
  colnames(intersecting_rows) <- c("regions2.row_index", "regions1.row_index")

  # Index of first row in regions2 which corresponds to a an overlap particular interval from regions1
  #first_element_of_intersection <- aggregate(1:nrow(intersecting_rows), by=list(intersecting_rows[,2]), function(x) {x[1]})[,2]
  #last_element_of_intersection <- aggregate(1:nrow(intersecting_rows), by=list(intersecting_rows[,2]), function(x) {x[length(x)]})[,2]
  
  # For each interval in regions1 this is start and end in regions1
  regions1.first_start <- intersecting_rows_regions1$start
  regions1.last_end <- intersecting_rows_regions1$end
  
  # For each interval in regions1 this is start and end of edge intervals in regions2
  regions2.first_start <- intersecting_rows_regions2$start
  regions2.last_end <- intersecting_rows_regions2$end
  
  intersecting_rows_regions2$start <- rowMax(cbind(regions1.first_start, regions2.first_start))
  intersecting_rows_regions2$end <- rowMin(cbind(regions1.last_end, regions2.last_end))
  
  intersection <- intersecting_rows_regions2

  return(intersection)
}

findIntersection <- function(tabs)
{ 
  chrs <- do.call(union, sapply(tabs, function(x) {x$chr}))
  
  tabs.splitted.per_chr <- as.list(vector(length = length(chrs)))
  names(tabs.splitted.per_chr) <- chrs  
  for (i in 1:length(tabs))
  {
    tmp.splitted <- split(tabs[[i]], tabs[[i]]$chr)
    
    for (chr in names(tmp.splitted[[i]]))
    { 
      tabs.splitted.per_chr[[chr]] <- c(tabs.splitted.per_chr[[chr]],  tmp.splitted[[i]][[chr]])
    }
  }
  
  joined_tabs.splitted <- as.list(vector(length = length(chrs)))
  names(joined_tabs.splitted) <- chrs  
  for (chr in names(joined_tabs.splitted))
  {
    intersection <- tabs.splitted.per_chr[[chr]][[1]]
    for (j in 2:length(tabs.splitted.per_chr[[chr]]))
    {
      intersection <- intersectChromosomalRegions(intersection, tabs.splitted.per_chr[[chr]][[j]])
    }
    joined_tabs.splitted[[chr]] <- intersection
  }
  
  joined_tabs <- unsplit(joined_tabs, chrs)

  return(joined_tabs)
}