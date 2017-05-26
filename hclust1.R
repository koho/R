hclust1 = function(d, method = "complete") {
  method = match.arg(method, c("ward.D", "single", "complete", "average", 
                               "mcquitty", "median", "centroid", "ward.D2"))
  size = attr(d, "Size")
  tind = rep(1:(size - 1), (size - 1):1) + size * unlist(lapply(1:(size - 1), `:`, size - 1))
  ind2sub = function(index, n) { c(tind[index] %% n, tind[index] %/% n + 1) }
  sub2ind = function(sub1, sub2, n = size) {
    sub = apply(rbind(sub1, sub2), 2, sort);
    sapply((sub[2,] - 1) * n + sub[1,], function(index) { which(tind == index) })
  }
  
  cdist = function(a, b, method) {
    distmat = outer(a, b, function(x, y) { d[sub2ind(x, y)] })
    switch(method, 
           "complete" = max(distmat),
           "single" = min(distmat),
           "average" = sum(distmat) / (length(a) * length(b)))
  }
  
  merged = matrix(nrow = size - 1, ncol = 2)
  height = numeric(size - 1)
  clusters = lapply(1:size, eval)
  
  for(k in seq_len(size - 1)) {
    mincluster = c(0, 0, Inf)
    for(i in seq_len(length(clusters) - 1)) {
      for( j in (i + 1):length(clusters)) {
        proximity = cdist(clusters[[i]], clusters[[j]], method)
        if(proximity < mincluster[3]) {
          mincluster = c(i, j, proximity)
        }
      }
    }
    
    merged[k,] = sort(sapply(mincluster[-3], function(ind) { ifelse(length(clusters[[ind]]) > 1, 
                                                  as.integer(names(clusters[ind])), 
                                                        -clusters[[ind]]) }))
    if(all(merged[k,] < 0)) merged[k,] = rev(merged[k,])
    height[k] = mincluster[3]
    clusters[[mincluster[2]]] = c(clusters[[mincluster[1]]], clusters[[mincluster[2]]])
    names(clusters)[mincluster[2]] = k
    clusters[mincluster[1]] = NULL
  }
  structure(list(merge = merged, height = height, labels = attr(d, "Labels"), method = method, 
            call = match.call(), dist.method = attr(d, "method")),
            class = "hclust1")
}

print.hclust1 = function(hc) { getS3method("print", "hclust")(hc) }

cutree1 = function(hc, k) {
  if(k < 1) stop("k must be greater than or equal to 1")
  cl = list()
  tree = hc$merge
  for(i in seq_len(nrow(tree) + 1 - k)) {
    if(prod(tree[i,]) < 0) {
      oldcl = max(tree[i,])
      cl[[as.character(oldcl)]] = c(-min(tree[i,]), cl[[as.character(oldcl)]])
      names(cl)[names(cl) == as.character(oldcl)] = as.character(i)
    } else if(all(tree[i,] > 0)) {
      cl[[as.character(tree[i,2])]] = c(cl[[as.character(tree[i,1])]], cl[[as.character(tree[i,2])]])
      names(cl)[names(cl) == as.character(tree[i,2])] = as.character(i)
      cl[as.character(tree[i,1])] = NULL
    } else {
      cl = c(cl, list(-tree[i,]))
      names(cl)[length(cl)] = as.character(i)
    }
  }
  ans = numeric(nrow(tree) + 1)
  cl = c(cl, lapply(setdiff(seq_len(length(ans)), unlist(cl)), eval))
  clabels = rep(rank(sapply(cl, min)), sapply(cl, length))
  ans[unlist(cl)] = clabels
  names(ans) = hc$labels
  ans
}

plot.hclust1 = function(x, k = 1, hang = 0.1,
                        main = "Cluster Dendrogram", labels = NULL,
                        sub = NULL, xlab = NULL, ylab = "Height", ...) {
  tree = x$merge
  uppernodes = if(k == 1) integer(0) else (nrow(tree) + 2 - k):nrow(tree)
  cols = matrix(nrow = length(uppernodes), ncol = 2)
  curcol = 1
  for(i in seq_along(uppernodes)) {
    cols[i,] = sapply(tree[uppernodes[i],], function(node) {
      if(node < 0) curcol <<- curcol + 1 else if(node %in% uppernodes) 1 else curcol <<- curcol + 1
    })
  }
  if(is.null(labels)) { 
    labels = if(is.null(x$labels)) as.character(1:(nrow(tree) + 1)) else x$labels
  } else {
    if(is.logical(labels) && length(labels) == 1) {
      labels = if(labels) if(is.null(x$labels)) as.character(1:(nrow(tree) + 1)) else x$labels else character(nrow(tree) + 1)
    } else {
      labels = as.character(labels)
    }
  }
  nodeheight = hang * diff(range(x$height))
  txtmar = 0.03 * diff(range(x$height))
  nodes = 0
  plot.new()
  plot.window(xlim = c(0.8, nrow(tree) + 1), ylim = c(ifelse(hang < 0, -max(nchar(labels)) * txtmar, min(x$height) - nodeheight - max(nchar(labels)) * txtmar), max(x$height)))
  axis(side = 2, at = pretty(range(x$height)))
  #par("usr")
  rp = function(l, r, x, height, col) {
    tree = x$merge
    lpos = if(l > 0) rp(tree[l, 1], tree[l, 2], x, x$height[l], 
                        if(l %in% uppernodes) cols[l - uppernodes[1] + 1,] else rep(col[1], 2)) else c(nodes <<- nodes + 1, ifelse(hang < 0, 0, height - nodeheight))
    rpos = if(r > 0) rp(tree[r, 1], tree[r, 2], x, x$height[r],
                        if(r %in% uppernodes) cols[r - uppernodes[1] + 1,] else rep(col[2], 2)) else c(nodes <<- nodes + 1, ifelse(hang < 0, 0, height - nodeheight))
    #lines(rep(c(lpos[1], rpos[1]), c(2, 2)), c(lpos[2], height, height, rpos[2]))
    lines(c(lpos[1], lpos[1], (lpos[1] + rpos[1]) / 2), c(lpos[2], height, height), col = col[1])
    lines(c((lpos[1] + rpos[1]) / 2, rpos[1], rpos[1]), c(height, height, rpos[2]), col = col[2])
    if(l < 0) text(lpos[1], lpos[2], labels[-l], srt = 90, col = col[1], adj = c(1.1, 0.5))
    if(r < 0) text(rpos[1], rpos[2], labels[-r], srt = 90, col = col[2], adj = c(1.1, 0.5))
    c((lpos[1] + rpos[1]) / 2, height)
  }
  rp(tree[nrow(tree), 1], tree[nrow(tree), 2], x, x$height[nrow(tree)], if(k == 1) c(1, 1) else cols[length(uppernodes),])
  if(is.null(xlab)) xlab = deparse(x$call[[2L]])
  if(is.null(sub)) sub = paste0(deparse(x$call[[1L]])," (*, \"", x$method,"\")")
  title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
  invisible()
}
