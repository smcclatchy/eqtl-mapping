### QTL heatmap
qtl_heatmap = function(qtl, map, low.thr = 3, high.thr = 8) {

  # Set values below low.thr = 0 and values above high.thr = high.thr.
  # This prevents the color scale from being dominated by high QTL and
  # minimizes the low end noise in the plot.
  qtl[qtl < low.thr] = 0
  qtl[qtl > high.thr] = high.thr
  
  # Set up breaks and colors.
  breaks = 0:800/100
  colors = colorRampPalette(brewer.pal(9, "YlOrBr"))(length(breaks) - 1)
  
  # Set up genome Mb positions of markers and Chr names.
  gmb = data.frame(chr = rep(names(map), sapply(map, length)), pos = unlist(map))
  chrmax = sapply(map, max)
  chrmid = cumsum(c(0, chrmax[-length(chrmax)])) + chrmax / 2
  for(i in 2:length(chrmax)) {
    wh = which(gmb$chr == names(chrmax)[i])
    gmb$pos[wh] = gmb$pos[wh] + sum(chrmax[1:(i-1)])
  }
  
  # Cluster phenotypes.
  cl = hclust(as.dist(1.0 - cor(qtl)), method = "average")
  qtl = qtl[,cl$order]
  
  # Plot heatmap.
  layout(matrix(1:2, 1, 2), widths = c(0.9, 0.1))
  par(plt = c(0.15, 0.99, 0.1, 0.95))
  image(gmb$pos, 1:ncol(qtl), qtl, col = colors, breaks = breaks, axes = F, ann = F)
  mtext(side = 1, line = 0.5, at = chrmid, text = names(chrmax), cex = 1.5)
  mtext(side = 2, line = 0.5, at = 1:ncol(qtl), text = colnames(qtl), las = 2)
  abline(h = 1:nrow(qtl) + 0.5, col = "grey50")
  abline(v = cumsum(chrmax) + 0.5, col = "grey50")
  box()
  
  # Plot color scale.
  scale.mat = matrix(breaks, nrow = 1)
  par(plt = c(0.4, 0.8, 0.1, 0.95))
  image(1:nrow(scale.mat), 1:ncol(scale.mat), scale.mat, breaks = breaks, col = colors,
        ann= F, axes = F)
  axis(side = 2, at = 0:8 * 100, labels = c(0:7, "> 8"), las = 2)
  
} # qtl_heatmap()