MultiTargetGetEstimate <- function(mat, design, coef, B=NULL, permutations=NULL,
                         full=FALSE) {
    v <- design[,coef]
    A <- design[,-coef, drop=FALSE]
    qa <- qr(A)
    S <- diag(nrow(A)) - tcrossprod(qr.Q(qa))

    vv <- if (is.null(B)) as.matrix(v) else{
        if(is.null(permutations)){
            if (length(coef) == 1)
            {
              replicate(B, sample(v))
            } else {
              do.call(cbind, replicate(B, v[sample(nrow(v)),], simplify=F))
            }
        } else{
            apply(permutations,2,function(i) v[i])
        }
    }
    
    # Make a matrix with normalizing coefficients. Normalizing coefficients should be coherent with design matrix
    multipliers <- apply(vv, 2, function(col)
                    { 
                      cases <- which(col > 0)
                      controls <- which(col < 0)
                      
                      col_with_multipliers <- col
                      col_with_multipliers[cases] <- 1/length(cases)
                      col_with_multipliers[controls] <- -1/length(controls)
                      col_with_multipliers
                    })
    
    
    #sv <- S %*% vv
    #vsv <- diag(crossprod(vv,sv))
    #b <- (mat %*% crossprod(S, vv)) / vsv
    
    b <- (mat %*% crossprod(S, multipliers))
    b <- abs(b)
    
    # If there are several case types, we search for bumps in them simultaneously. 
    # Meth differences in different case types are summed.
    if (length(coef) > 1)
    {
      new_b <- b[,seq(1,ncol(b),length(coef))]
      
      for (i in 2:length(coef))
      {
        new_b <- new_b + b[,seq(i,ncol(b),length(coef))]
      }
    }
    
    b <- new_b
    
    if(!is.matrix(b))
        b <- matrix(b, ncol = 1)
    if (full) {
        sy <- mat %*% S
        df.residual <- ncol(mat) - qa$rank - 1
        if (is.null(B)) {
            sigma <- matrix(sqrt(rowSums((sy - tcrossprod(b, sv))^2) / df.residual), ncol=1)
        } else {
            sigma <- b
            tmp <- sy
            for (j in 1:B) {
                tmp <- tcrossprod(b[,j], sv[,j])
                sigma[,j] <- rowSums((sy-tmp)^2)
            }
            sigma <- sqrt(sigma/df.residual)
        }
        out <- list(coef=as.matrix(b),
                    sigma=sigma,
                    stdev.unscaled=sqrt(1/vsv),
                    df.residual=df.residual)
        if (is.null(B)) out$stdev <- as.numeric(out$stdev)
    } else {
        out <- as.matrix(b)
    }
  return(out)
}

MultiTargetGetModT <- function(obj) {
    s2 <- apply(obj$sigma^2, 2, limma::squeezeVar, obj$df.residual)
    out <- obj$coef
    for (j in 1:ncol(out)) {
        out[,j] <- out[,j] / obj$stdev.unscaled[j] / sqrt(s2[[j]]$var.post)
    }
    df.total <- obj$df.residual + sapply(s2,"[[","df.prior")
    return(list(t=out, df.total=df.total))
}
