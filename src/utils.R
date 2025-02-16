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
      dat[[x]][train,,drop=F]
    } else {
      dat[[x]][,train,drop=F]
    }
  })
  dat.test <- sapply(simplify = F, names(dat), function(x) {
    if(x == 'clin') {
      dat[[x]][test,,drop=F]
    } else {
      dat[[x]][,test,drop=F]
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


# plot bi-partite visualisation of labels (useful for comparing how different cluster memberships look)
plot_cluster_comparison <- function(df1, df2) {
  df1 <- data.frame(df1, check.names = F)
  df2 <- data.frame(df2, check.names = F)
  if(sum(is.na(df1[,1])) > 0) {
    warning("Converting NA values to 'Undefined'")
    df1[is.na(df1[,1]),1] <- 'Undefined'
  }
  if(sum(is.na(df2[,1])) > 0) {
    warning("Converting NA values to 'Undefined'")
    df2[is.na(df2[,1]),1] <- 'Undefined'
  }
  df1[,1] <- as.factor(df1[,1])
  df2[,1] <- as.factor(df2[,1])

  dt <- data.table(merge(df1[,1,drop=F], df2[,1,drop = F], by = 'row.names'))
  labels <- colnames(dt[,2:3])
  colnames(dt) <- c('rn', 'g1', 'g2')
  dt1 <- dt[order(g1), c('rn', 'g1')]
  dt2 <- dt[order(g2), c('rn', 'g2')]
  dt1$r1 <- 1:nrow(dt1)
  dt2$r2 <- 1:nrow(dt2)
  dt <- merge(dt1, dt2, by = 'rn')[order(rn)]
  ami <- aricode::AMI(dt$g1, dt$g2)
  label_pos_left <- -0.05
  label_pos_right <- 1.05
  # plot segments
  ggplot(dt[order(g1)]) +
    geom_point(aes(x = 0, y = r1, color = g1))  +
    geom_point(aes(x = 1, y = r2, color = g2)) +
    geom_segment(aes(x = 0, xend = 1, y = r1, yend = r2, color = g1), alpha = 0.25) +
    geom_label(data = dt[,median(r1), by = g1], aes(x = label_pos_left, y = V1, label = g1, color = g1), hjust = 1) +
    geom_segment(data = dt[,list('max' = max(r1), 'min' = min(r1), 'median' = median(r1)), by = g1],
                 aes(x = label_pos_left, y = median, xend = 0, yend = max, color = g1)) +
    geom_segment(data = dt[,list('max' = max(r1), 'min' = min(r1), 'median' = median(r1)), by = g1],
                 aes(x = label_pos_left, y = median, xend = 0, yend = min, color = g1)) +
    geom_label(data = dt[,median(r2), by = g2], aes(x = label_pos_right, y = V1, label = g2, color = g2), hjust = 0) +
    geom_segment(data = dt[,list('max' = max(r2), 'min' = min(r2), 'median' = median(r2)), by = g2],
                 aes(x = label_pos_right, y = median, xend = 1, yend = max, color = g2)) +
    geom_segment(data = dt[,list('max' = max(r2), 'min' = min(r2), 'median' = median(r2)), by = g2],
                 aes(x = label_pos_right, y = median, xend = 1, yend = min, color = g2)) +
    annotate('label', x = 0, y = max(dt$r1)+1, label = labels[1], vjust = 0) +
    annotate('label', x = 1, y = max(dt$r2)+1, label = labels[2], vjust = 0) +
    ggtitle(label = paste0(labels, collapse = " <-> "),
            subtitle = paste0("Adjusted Mutual Information: ",round(ami,2))) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
          legend.position = 'none',
          plot.margin = unit(c(0.1, 1.2, 0, 1.2), units = 'in'),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    coord_cartesian(clip = 'off')
}
