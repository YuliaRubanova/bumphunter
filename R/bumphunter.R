setMethod("bumphunter", signature(object = "matrix"),
          function(object, design, chr=NULL, pos, cluster=NULL,
                   coef=2, cutoff=NULL, pickCutoff=FALSE, pickCutoffQ=0.99,
                   maxGap=500, smooth=FALSE, smoothFunction=locfitByCluster,
                   useWeights=FALSE,  B=ncol(permutations), permutations=NULL,
                   verbose=TRUE, ...){
            if(missing(design)) stop("design must be specified")
            if(missing(pos)) stop("If object is a matrix, pos must be specified")
            bumphunterEngine(object, design=design, chr=chr, pos, cluster=cluster,
                             coef=coef,
                             cutoff=cutoff, pickCutoff=pickCutoff, pickCutoffQ=pickCutoffQ,
                             maxGap=maxGap,smooth=smooth,
                             smoothFunction=smoothFunction,
                             useWeights=useWeights, B=B,
                             permutations=NULL,
                             verbose=verbose, ...)
          })

bumphunterEngine <- function(mat, design, chr=NULL, pos, cluster=NULL,
                             coef=2,
                             cutoff=NULL, pickCutoff=FALSE, pickCutoffQ=0.99, 
                             maxGap=500, smooth=FALSE,
                             smoothFunction=locfitByCluster,
                             useWeights=FALSE, B=ncol(permutations), permutations=NULL,
                             verbose=TRUE, ...){
    ## cutoff = c(L, U)

    if(is.null(B)) B=0
    if(!is.matrix(permutations) & !is.null(permutations))
        stop("permutations must be NULL or a matrix.")
    if(!is.null(permutations)){
        if(nrow(design)!=nrow(permutations))
            stop("Number of rows of 'design' must match number of rows of 'permutations'.")
        if(B!=ncol(permutations)){
            warning("Ignoring the supplied B. Using 'ncol(permutations)' instead.")
            B=ncol(permutations)
        }
    }
    
    if(ncol(design)>2 & B>0)
        warning("The use of the permutation test (B>0), is not recommended with multiple covariates, (ncol(design)>2). See vignette for more information.")
    
    if(!is.matrix(mat))
        stop("'mat' must be a matrix.")
    if(ncol(mat) != nrow(design))
        stop("Number of columns of 'mat' must  match number of rows of 'design'")
    if(B < 0) stop("'B' has to be an integer greater or equal to zero")

    if(! (is.null(cutoff) || length(cutoff) %in% 1:2))
        stop("'cutoff' has to be either NULL or a vector of length 1 or 2")

    if(length(cutoff) == 2)
        cutoff <- sort(cutoff)

    if(is.null(cutoff) && !pickCutoff)
        stop("Must pick a cutoff, either using 'cutoff' or 'pickCutoff'")
    if(!is.null(cutoff) && pickCutoff) {
        pickCutoff <- FALSE
        warning("'cutoff' was supplied so ignoring 'pickCutoff=TRUE'")
    }
    
    if(pickCutoff && (length(pickCutoffQ) != 1 || pickCutoff < 0 || pickCutoffQ > 1))
        stop("Using `pickCutoff = TRUE' requires that 'pickCutoffQ' is a single number between 0 and 1")
    if(pickCutoff && B<1) stop("Must do at least one permution to pick a cutoff")

    if (!getDoParRegistered())
        registerDoSEQ()
    workers <- getDoParWorkers()
    backend <- getDoParName()
    version <- getDoParVersion()
    subverbose <- max(as.integer(verbose) - 1L, 0)
    if (verbose) {
        if (workers == 1) {
            mes <- "[bumphunterEngine] Using a single core (backend: %s, version: %s)."
            message(sprintf(mes, backend, version))
        } else {
            mes <- "[bumphunterEngine] Parallelizing using %s workers/cores (backend: %s, version: %s)."
            message(sprintf(mes, workers, backend, version))
        }
    }
    
    ##B is the number of random samples we take
    if (is.null(chr))
        chr <- rep("Unspecified", length(pos))
    if(is.factor(chr))
        chr <- as.character(chr)
    if(is.null(cluster))
        cluster <- clusterMaker(chr, pos, maxGap=maxGap)
    
    ## if(B>0 & B<workers) stop("B must be bigger than workers (or 0)")
    
    if(verbose) message("[bumphunterEngine] Computing coefficients.")
    
    if(useWeights & smooth){
        tmp <- .getEstimate(mat=mat, design=design, coef=coef, full=TRUE)
        rawBeta <- tmp$coef
        weights <- tmp$sigma
        rm(tmp)
    } else{
        rawBeta <- .getEstimate(mat=mat, design=design, coef=coef, full=FALSE)
        weights <- NULL
    }
  
    if(smooth){
        if(verbose) message("[bumphunterEngine] Smoothing coefficients.")
        beta <- smoother(y=rawBeta, x=pos, cluster=cluster, weights=weights,
                         smoothFunction=smoothFunction, verbose=subverbose, ...) 
        Index <- which(beta$smoothed)
        beta <- beta$fitted
    } else {
        beta <- rawBeta
        Index <- seq(along=beta)
    }
    
    if(B>0){
        if (verbose) message("[bumphunterEngine] Performing ", B, " permutations.")
        if(useWeights && smooth) {
            tmp <- .getEstimate(mat, design, coef, B, permutations, full=TRUE)
            permRawBeta <- tmp$coef
            weights <- tmp$sigma
            rm(tmp)
        } else {
            permRawBeta <- .getEstimate(mat, design, coef, B, permutations, full=FALSE)
            weights <- NULL
        }
        
        ## Get individual p-values based on permutation of samples
        ## For each permutation we consider whether the absolute value of
        ## the observed beta is 
        if(verbose) message("[bumphunterEngine] Computing marginal permutation p-values.")
        
        sumGreaterOrEqual <- rowSums(greaterOrEqual(abs(permRawBeta), abs(as.vector(rawBeta))))
        pvs <- (sumGreaterOrEqual + 1L) / (B + 1L)
        
        if(smooth){
            if(verbose) message("[bumphunterEngine] Smoothing permutation coefficients.")
            permBeta <- smoother(y=permRawBeta, x=pos, cluster=cluster, weights=weights,
                                 smoothFunction=smoothFunction,
                                 verbose=subverbose,...)$fitted
        } else permBeta <- permRawBeta
        
        if(is.null(cutoff))
            cutoff <- quantile(abs(permBeta), pickCutoffQ, na.rm=TRUE)
        if(verbose) message(sprintf("[bumphunterEngine] cutoff: %s", round(cutoff,3)))
    } ## Done with permutations

    if(verbose) message("[bumphunterEngine] Finding regions.")
    tab <- regionFinder(x=beta, chr=chr, pos=pos, cluster=cluster,
                        cutoff=cutoff, ind=Index, verbose=FALSE)
    if (nrow(tab)==0) {
        if (verbose) message("[bumphunterEngine] No bumps found!")
        return(list(table=NA, coef=rawBeta, fitted=beta, pvaluesMarginal=NA))
    } else {
        if (verbose) message(sprintf("[bumphunterEngine] Found %s bumps.", nrow(tab)))
    }

    
    if (B<1) {
        return(list(table=tab, coef=rawBeta, fitted=beta, pvaluesMarginal=NA))
    }
    
    if (verbose) message("[bumphunterEngine] Computing regions for each permutation.")
    chunksize <- ceiling(B/workers)
    subMat <- NULL
    nulltabs <- foreach(subMat=iter(permBeta, by="col", chunksize=chunksize),
                        .combine="c", .packages = "bumphunter") %dorng% {
        apply(subMat, 2, regionFinder, chr=chr, pos=pos, cluster=cluster,
              cutoff=cutoff, ind=Index, verbose=FALSE)
    }
    attributes(nulltabs)[["rng"]] <- NULL
    
    if (verbose) message("[bumphunterEngine] Estimating p-values and FWER.")

    L <- V <- A <- as.list(rep(0, B))
    for(i in 1:B) {
        nulltab <- nulltabs[[i]]
        if (nrow(nulltab)>0) { 
            L[[i]] <- nulltab$L
            V[[i]] <- nulltab$value
            A[[i]] <- nulltab$area
        }
    }
    ## for observed length and height
    ## compute the total compute total number of times
    ## it is seen in permutations

    computation.tots <- function(tab, V, L) {
        Lvalue <- cbind(tab$L, abs(tab$value))
        chunksize <- ceiling(nrow(Lvalue)/workers)
        subL <- NULL
        tots <- foreach(subL=iter(Lvalue, by="row", chunksize=chunksize),
                        .combine="cbind", .packages = "bumphunter") %dorng% {
                            apply(subL,1,function(x) {
                                res <- sapply(seq(along=V), function(i) {
                                    sum(greaterOrEqual(L[[i]], x[1]) &
                                        greaterOrEqual(abs(V[[i]]), x[2]))
                                })
                                c(mean(res>0),sum(res))
                            })
                        }
        attributes(tots)[["rng"]] <- NULL
        rate1 <- tots[1,]
        pvalues1 <- tots[2,] / sum(sapply(nulltabs,nrow))
        return(list(rate1 = rate1, pvalues1 = pvalues1))
    }

    ptime1 <- proc.time()
    comp <- computation.tots(tab = tab, V=V, L=L)
    rate1 <- comp$rate1
    pvalues1 <- comp$pvalues1
    ptime2 <- proc.time()
    ## cat(sprintf("tots timing: %s\n", (ptime2 - ptime1)[3]))

    ## ## Old tots, keep around since it is easier to read, for vectorizing
    ## tots <- sapply(seq(along=V), function(i) {
    ##     apply(Lvalue, 1, function(x) {
    ##         sum(greaterOrEqual(L[[i]], x[1]) &
    ##             greaterOrEqual(abs(V[[i]]), x[2]))
    ##     })
    ## })
    ## if (is.vector(tots)) {
    ##     tots <- matrix(tots, nrow=1)
    ## }
    ## ## This is like a FWER
    ## rate1 <- rowMeans(tots>0)
    
    ## ## Now compute pvalues by assuming everything is exchangeable
    ## pvalues1 <- rowSums(tots) / sum(sapply(nulltabs,nrow))
    
    computation.tots2 <- function(tab, A) {
        Avalue <- matrix(tab$area, ncol=1)
        chunksize <- ceiling(nrow(Avalue)/workers)
        subA <- NULL
        tots2 <- t(foreach(subA=iter(Avalue,by="row",chunksize=chunksize),
                           .combine="cbind",.packages="bumphunter") %dorng% {
                               sapply(subA, function(x) {
                                   return(sapply(seq(along = A), function(i) {
                                       sum(greaterOrEqual(A[[i]], x[1]))
                                   }))    
                               })
                           })
        if (is.vector(tots2)) {
            tots2 <- matrix(tots2, nrow=1)
        }
        rate2 <- rowMeans(tots2>0)
        pvalues2 <- rowSums(tots2)/sum(sapply(nulltabs,nrow))
        return(list(rate2 = rate2, pvalues2 = pvalues2))
    }

    ptime1 <- proc.time()
    comp <- computation.tots2(tab = tab, A = A)
    rate2 <- comp$rate2
    pvalues2 <- comp$pvalues2
    ptime2 <- proc.time()
    ## cat(sprintf("tots2 timing: %s\n", (ptime2 - ptime1)[3]))

    ##    tots2 <- sapply(seq(along=A), function(i) {
    ##    sapply(tab$area, function(x) {
    ##        sum(greaterOrEqual(A[[i]], x[1]))
    ##    })
    ## })

    tab$p.value <- pvalues1
    tab$fwer <- rate1
    
    tab$p.valueArea <- pvalues2
    tab$fwerArea <- rate2
    
    tab <- tab[order(tab$fwer,-tab$area),]

    algorithm <- list(version = packageDescription("bumphunter")$Version, coef = coef,
                     cutoff = cutoff, pickCutoff = pickCutoff, pickCutoffQ = pickCutoffQ,
                     smooth = smooth, maxGap = maxGap, B = B,
                      permutations=permutations,
                      useWeights = useWeights,
                      smoothFunction = deparse(substitute(smoothFunction)))
    ret <- list(table=tab, coef=rawBeta, fitted=beta,
                pvaluesMarginal=pvs, null=list(value=V,length=L), algorithm = algorithm)
    class(ret) <- "bumps"
    return(ret)
}

print.bumps <- function(x, ...) {
    cat(sprintf("a 'bumps' object with %s bumps\n", nrow(x$table)))
}
    