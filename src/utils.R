# common functions used in different analyses
compute_PCA <- function(exp, topN = 50) {
  top <- head(names(sort(apply(exp, 1, sd), decreasing = T)), 5000)
  M <- t(exp[top,])
  pca <- stats::prcomp(M, rank. = topN)
  return(pca)
}


# given a matrix of measurements (rows: patients, columns: different conditions/covariates)
# and a factor vector (representing cluster membership of each patient), 
# find out which covariates show differential values between pairs of factors 
get_differential_factors <- function(M, factors, Nodes = 10) {
  cl <- parallel::makeForkCluster(nnodes = Nodes)
  parallel::clusterExport(cl = cl, varlist = c('M', 'factors'), envir = environment())
  res <- do.call(rbind, pbapply::pblapply(cl = cl, X = colnames(M), FUN = function(x) {
    l <- split(M[,x], factors)
    do.call(rbind, lapply(names(l), function(i) {
      do.call(rbind, lapply(names(l), function(j) {
        g1 <- l[[i]]
        g2 <- l[[j]]
        t <- wilcox.test(g1, g2, alternative = 'greater')
        fc <- log2((mean(g1)+1) / (mean(g2)+1))
        data.frame("variable" = x, "ref_cl" = i, "target_cl" = j, "pval" = t$p.value, "log2fc" = fc)
      }))
    }))
  }))
  parallel::stopCluster(cl)
  res$padj <- p.adjust(res$pval, method = 'BH')
  # find markers that are really specific for each cluster and differential compared to every other cluster
  res <- data.table::as.data.table(res)
  k <- length(unique(factors))
  # find those markers that are differential compared to all target clusters
  res_sig <- res[padj < 0.05, length(unique(target_cl)), by = c('variable', 'ref_cl')][V1 == (k-1)]
  
  res_sig <- lapply(split(res_sig, res_sig$ref_cl), function(m) {
    m <- merge(m[,1:2], res, by = c('variable', 'ref_cl'))[order(padj)]
    m <- m[ref_cl != target_cl, .SD[which.min(padj)],by = variable][order(padj)][,-3]
  })
  
  return(res_sig)
}

# M: matrix, samples on rows, features on columns 
plot_umap <- function(M, factors, returnData = FALSE) {
  require(umap)
  umap.df <- as.data.frame(umap::umap(M)[['layout']])
  umap.df$group <- factors
  colnames(umap.df)[1:2] <- c('UMAP1', 'UMAP2')
  if(returnData == TRUE) {
    return(umap.df)
  }
  ggplot(umap.df, aes(x = UMAP1, y = UMAP2)) + 
    geom_point(aes(color = group))
}

# get the center coordinate of samples in a matrix group by a factor vector
get_centroid <- function(M, factors) {
  centroids <- data.frame(
      pbapply::pbapply(M, 2, function(x) {
        tapply(X = x, factors, FUN = median)
      })
    )
  colnames(centroids) <- c('x', 'y')
  centroids$text.label <- as.factor(rownames(centroids))
  return(centroids)
}
# get a basis matrix that is the mean value per factor of variable
get_basis_matrix <- function(M, factors) {
  B <- pbapply::pbapply(M, 2, function(x) {
    tapply(X = x, factors, FUN = mean)
  })
  return(B)
}

# M: matrix, samples on rows, features on columns 
plot_tsne <- function(M, factors = NULL, show.labels = F,
		      label.size = 4, perplexity = 30, 
                      jitter = FALSE, returnData = FALSE) {
  require(Rtsne)
  if(jitter == TRUE) {
    M <- jitter(M)
  }
  to_keep <- which(apply(M, 1, sd) != 0)
  if(length(to_keep) < nrow(M)) {
    warning("Removing 0 sd samples N=",nrow(M) - length(to_keep))
    M <- M[to_keep, ]
    factors <- factors[to_keep]
  }
  
  tsne <- Rtsne(M, perplexity = perplexity) #as.data.frame(umap::umap(M)[['layout']])
  df <- as.data.frame(tsne$Y)
  rownames(df) <- rownames(M)
  df$group <- factors
  colnames(df)[1:2] <- c('tSNE1', 'tSNE2')
  if(returnData == TRUE) {
    return(df)
  } 
  if(is.null(factors)) {
    return(df)
  }
  p <- ggplot(df, aes(x = tSNE1, y = tSNE2)) + 
    geom_point(aes(color = group))
  if(show.labels == TRUE) {
    centroids <- data.frame(get_basis_matrix(tsne$Y, factors))
    colnames(centroids) <- c('x', 'y')
    centroids$group <- as.factor(rownames(centroids))
    p <- p + geom_text(data = centroids, aes(x = x, y = y, 
                                        label = group), 
                       size = label.size)
  }
  return(p)
}

# M: matrix, samples on columns, features on rows
plot_pca <- function(M, factors) {
  pca <- compute_PCA(M, topN = 2)
  df <- data.frame(pca$x, check.names = FALSE, row.names = rownames(pca$x))
  df$group <- factors

  var_exp <- round(diag(cov(pca$x))/sum(diag(cov(pca$x))) * 100, 1)
  ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = group), size = 5, alpha = 0.5) +
    theme_bw(base_size = 8) +
    labs(x = paste0('PC1 (',var_exp[['PC1']],'%)'),
         y = paste0('PC2 (',var_exp[['PC2']],'%)')) 
}

# get pipeline output for a given prefix
get_data <- function(prefix, task, workdir, data_types = NULL) {
  # get sample embeddings from the best flexynesis model
  Imp <- data.table::fread(file.path(workdir, 'results', paste0(prefix, ".feature_importance.csv")))
  
  # embeddings
  E <- rbind(
    read.csv(file.path(workdir, 'results', paste0(prefix, ".embeddings_train.csv")), row.names = 1),
    read.csv(file.path(workdir, 'results', paste0(prefix, ".embeddings_test.csv")), row.names = 1)
  )
  # metadata
  colData <- rbind(
    read.csv(file.path(workdir, 'data', task, 'train', 'clin.csv')),
    read.csv(file.path(workdir, 'data', task, 'test', 'clin.csv'))
  )[rownames(E),]
  
  # predicted labels 
  pred_labels <- data.table::fread(file.path(workdir, 'results',
                                             paste0(prefix, ".predicted_labels.csv")))
  
  return(list('Imp' = Imp, 'E' = E, 'colData' = colData, 'pred_labels' = pred_labels))
}

get_pca <- function(exp, colData, topN = 50) {
  top <- head(names(sort(apply(exp, 1, sd), decreasing = T)), 5000)
  M <- t(exp[top,])
  fit <- stats::prcomp(M, rank. = topN)
  df <- as.data.frame(fit[['x']])
  df <- merge(df, colData, by = 'row.names')
  var_exp <- round(diag(cov(fit$x))/sum(diag(cov(fit$x))) * 100, 1)
  var_exp <- sapply(names(var_exp), function(x) {
    paste0(x, " (", var_exp[[x]], "%)")
  })
  return(list('dat' = df, 'var_exp' = var_exp))
}


split_dat <- function(dat, samples = NULL, ratio = 0.7) {
  set.seed(42)
  if(is.null(samples)) {
    samples <- rownames(dat$clin)
  }
  cat(length(samples), "\n")
  train <- sample(samples, round(ratio * length(samples)))
  test <- setdiff(samples, train)
  dat.train <- sapply(simplify = F, names(dat), function(x) {
    if(x == 'clin') {
      dat[[x]][train,]
    } else {
      dat[[x]][,train]
    }
  })
  dat.test <- sapply(simplify = F, names(dat), function(x) {
    if(x == 'clin') {
      dat[[x]][test,]
    } else {
      dat[[x]][,test]
    }
  })
  return(list('train' = dat.train, 'test' = dat.test))
} 

# print train/test folders 
print_dataset <- function(dat, outdir) {
  if(!dir.exists(outdir)) {
    dir.create(outdir)
  }
  lapply(names(dat), function(s) {
    p <- file.path(outdir, s)
    if(!dir.exists(p)) {
      dir.create(p)
    }
    lapply(names(dat[[s]]), function(f) {
      write.table(dat[[s]][[f]], 
                  file = file.path(p, paste0(f, ".csv")), sep = ',')
    })
  })
}
