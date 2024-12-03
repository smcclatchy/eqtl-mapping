################################################################################
# Transcriptome map using ggplot framework.
# Given a set of eQTL locations and (optionally) gene locations, plot the QTL
# location versus the gene location. Genes must be identified by Ensembl ID.
# The function can look up gene locations by Ensembl ID or the user may provide
# them.
# Points may be optionally colored by LOD score.
#
# Daniel Gatti
# Yuka Takemon
# dan.gatti@jax.org
# Nov. 14, 2017
################################################################################
library(AnnotationHub)
library(rtracklayer)
library(tidyverse)

get_ensembl_genes = function() {

    hub = AnnotationHub()
    hub = query(hub, c("ensembl", "gtf", "mus musculus", "GRCm39"))
    hub = hub[grep("^Mus_musculus\\.GRCm39\\.[0-9]+\\.gtf$", hub$title)]
    latest.gtf = sort(hub$title)[length(hub)]
    ensembl = hub[[names(hub)[hub$title == latest.gtf]]]

    return(ensembl[ensembl$type == "gene"])

} # get_ensembl_genes()

# I had hoped to get the chr lengths from seqlengths(ensembl), but they
# don't always fill it in. I'm getting it from the maximum gene location
# on each chromosome.
get_chr_length = function(ensembl) {

  tmp = data.frame(chr = seqnames(ensembl), end = end(ensembl))
  tmp = tmp %>% filter(!(substring(chr, 1, 1) %in% c("G", "J"))) %>%
          group_by(chr) %>%
          summarize(len = max(end))
  chrlen = tmp$len * 1e-6
  names(chrlen) = tmp$chr

  return(chrlen)

} # get_chr_length()

# Get


# Arguments:
# data: data.frame (or tibble) with the following columns:
#       gene_id: (required) character string containing the Ensembl gene ID.
#       qtl_chr: (required) character string containing QTL chromsome.
#       qtl_pos: (required) floating point number containing the QTL position
#                in Mb.
#       qtl_lod: (optional) floating point number containing the LOD score.
#       gene_chr:  (optional) character string containing transcript chromosome.
#       gene_start: (optional) character string containing transcript start
#                 postion in Mb.
#       gene_end:  (optional) character string containing transcript end
#                position in Mb.
# color.points: logical that is TRUE if the points should be colored by LOD.
# cis.points: logical that is TRUE if the points should be colored if they
#             are with in cis.
# cis.radius: numeric value containing the radius in Mb between a gene and a cis-eQTL.
#             Optional.
# cis.color: color for cis QTL. Optional.
# Returns:
# a plot of the QTL and gene location for each gene.
ggtmap = function(data, color.points = FALSE, cis.points = FALSE, cis.radius = 2,
         cis.color = "#4286f4") {

  # Check for required column names.
  required.colnames = c("gene_id", "qtl_chr", "qtl_pos")

  if(all(!required.colnames %in% colnames(data))) {
    stop(paste("colnames must contain the following columns:",
         paste(required.colnames, collapse = ",")))
  } # if(!required.colnames %in% colnames(data))

  # Make sure that columns are not factors.
  data$gene_id = as.character(data$gene_id)
  data$qtl_chr = as.character(data$qtl_chr)

  gene.position.colnames = c("gene_chr", "gene_start", "gene_end")
  if(!all(gene.position.colnames %in% colnames(data))) {

    message(paste("Using Ensembl gene locations because optional gene",
            "position columns (", paste0(gene.position.colnames, collapse = ","),
            ") not found."))

	# Get the latest Ensembl GTF.
    ensembl = get_ensembl_genes()

    id    = ensembl$gene_id
    chr   = seqnames(ensembl)
    start = start(ensembl) * 1e-6
    end   = end(ensembl)   * 1e-6

    df = data.frame(gene_id = id, gene_chr = chr, gene_start = start,
                    gene_end = end, stringsAsFactors = F)
    data = left_join(data, df, by = "gene_id")

  } # if(gene.position.colnames %in% colnames(data))

  # Make sure that columns are not factors.
  data$gene_chr = as.character(data$gene_chr)

  # Get the gene mid-point.
  data = data %>% mutate(gene_pos = (gene_end + gene_start) * 0.5)

  # Fix the factor levels for the chr.
  all.chr = data %>% select(qtl_chr, gene_chr) %>%
            gather(k, v) %>%
            select(v) %>%
            distinct() %>%
            arrange(v)
  all.chr = all.chr$v[!is.na(all.chr$v)]

  if(length(grep("M", all.chr)) > 0) {
    wh = grep("M", all.chr)
    all.chr = all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
  }

  # Remove any NAs.
  data = data %>% na.omit

  data$qtl_chr  = factor(data$qtl_chr,  levels = all.chr[order(as.numeric(all.chr))])
  data$gene_chr = factor(data$gene_chr, levels = rev(all.chr[order(as.numeric(all.chr))]))

  # If we're plotting cis points, then add a cis-QTL column.
  if(cis.points) {

    data = data %>% mutate(cis = (gene_chr == qtl_chr) & (abs(gene_start - qtl_pos) <= cis.radius))
	cis.colors = c("black", cis.color)
	names(cis.colors) = c("FALSE", "TRUE")
    print(ggplot(data, aes(x = qtl_pos, y = gene_pos), alpha = 0.5) +
      geom_point(aes(color = cis), alpha = 0.5) +
      scale_color_manual(values = cis.colors) +
      facet_grid(gene_chr ~ qtl_chr, scales = "free", shrink = TRUE) +
      labs(x = "QTL Position", y = "Gene Position") +
      theme(panel.background = element_blank(),
	        panel.border = element_rect(fill = 0, color = "grey70"),
	        panel.grid.minor = element_blank(),
            panel.spacing = unit(0.05, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 1)))

  } else {

    print(ggplot(data, aes(x = qtl_pos, y = gene_pos)) +
      geom_point(aes(color = qtl_lod, alpha = 0.5)) + {
        if(color.points) scale_color_continuous(low = "grey50", high = "red")
      } +
      facet_grid(gene_chr ~ qtl_chr, scales = "free", shrink = TRUE) +
      labs(x = "QTL Position", y = "Gene Position") +
      theme(panel.background = element_blank(),
	  	    panel.border = element_rect(fill = 0, color = "grey70"),
	  	    panel.grid.minor = element_blank(),
            panel.spacing = unit(0.05, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 1)))

   } # else

} # ggtmap()


# Plot the density of eQTL along the chromosomes.
# Arguments:
# data: data.frame (or tibble) with the following columns:
#       ensembl: (required) character string containing the Ensembl gene ID.
#       qtl_chr: (required) character string containing QTL chromsome.
#       qtl_pos: (required) floating point number containing the QTL position
#                in Mb.
#       qtl_lod: (optional) floating point number containing the LOD score.
#       gene_chr:  (optional) character string containing transcript chromosome.
#       gene_start: (optional) character string containing transcript start
#                 postion in Mb.
#       gene_end:  (optional) character string containing transcript end
#                position in Mb.
# lod_thr: numeric value that is the LOD above which QTL will be retained.
#          Default = 7.
eqtl_density_plot = function(data, lod_thr = 7) {

  # Create a set of rolling breakpoints, 4 Mb apart.
  breaks = matrix(c(seq(0, 200, 4), seq(1, 201, 4), seq(2, 202, 4), seq(3, 203, 4)), ncol = 4)

  tmp = as.list(1:ncol(breaks))

  for(i in 1:ncol(breaks)) {
    tmp[[i]] = data %>%
                 filter(qtl_lod >= lod_thr) %>%
                 arrange(qtl_chr, qtl_pos) %>%
                 group_by(qtl_chr) %>%
                 mutate(win = cut(qtl_pos, breaks = breaks[,i])) %>%
                 group_by(qtl_chr, win) %>%
                 summarize(cnt = n()) %>%
                 separate(win, into = c("other", "prox", "dist")) %>%
                 mutate(prox = as.numeric(prox),
                        dist = as.numeric(dist),
                        mid  = 0.5 * (prox + dist)) %>%
                 dplyr::select(qtl_chr, mid, cnt)
  } # for(i)

  trans = bind_rows(tmp[[1]], tmp[[1]], tmp[[3]], tmp[[4]])
  rm(tmp)

  ggplot(trans, aes(mid, cnt)) +
    geom_line() +
    geom_hline(aes(yintercept = 100), linetype = 2, color = "grey50") +
    facet_grid(.~qtl_chr, scales = "free") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "grey70"),
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(angle = 90)) +
    labs(title = "eQTL Density Plot", x = "Mb", y = "Number of Transcripts")

} # eqtl_density_plot()
