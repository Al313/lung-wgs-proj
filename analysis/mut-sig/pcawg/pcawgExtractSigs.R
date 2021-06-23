
# Loading libraries
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(ggdendro)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(MutationalPatterns)

# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   
#   sbs_sig_profile_lung <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/mut-sig-ref/sbs/sbs-sig-profiles.rds")
#   dbs_sig_profile_lung <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/HMF_PCAWG/sigs_denovo/02_mergedCancerTypes/nmf_out/Lung.dbs/matrices/3/sig_profiles.txt.gz",
#                                    header = T, sep = "\t", stringsAsFactors = F)
#   id_sig_profile_lung <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/HMF_PCAWG/sigs_denovo/02_mergedCancerTypes/nmf_out/Lung.indel/matrices/06/sig_profiles.txt.gz",
#                                      header = T, sep = "\t", stringsAsFactors = F)
#   
#   } else {
#   
#   sbs_sig_profile_lung <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/mut-sig-ref/sbs/sbs-sig-profiles.rds")
#   dbs_sig_profile_lung <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/HMF_PCAWG/sigs_denovo/02_mergedCancerTypes/nmf_out/Lung.dbs/matrices/3/sig_profiles.txt.gz",
#                                    header = T, sep = "\t", stringsAsFactors = F)
#   id_sig_profile_lung <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/HMF_PCAWG/sigs_denovo/02_mergedCancerTypes/nmf_out/Lung.indel/matrices/06/sig_profiles.txt.gz",
#                                      header = T, sep = "\t", stringsAsFactors = F)
# 
# }
# 
# order_colnames <- order(as.numeric(sapply(strsplit(colnames(sbs_sig_profile_lung), split = "_", 2), tail, 1)))
# sbs_sig_profile_lung <- sbs_sig_profile_lung[,order_colnames]
# 
# order_colnames <- order(as.numeric(sapply(strsplit(colnames(dbs_sig_profile_lung), split = "_", 2), tail, 1)))
# dbs_sig_profile_lung <- dbs_sig_profile_lung[,order_colnames]
# 
# order_colnames <- order(as.numeric(sapply(strsplit(colnames(id_sig_profile_lung), split = "_", 2), tail, 1)))
# id_sig_profile_lung <- id_sig_profile_lung[,order_colnames]



if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
} else {
  
  library(mutSigExtractor)
}


if (dir.exists("/hpc/cuppen/")){
  pcawg_meta <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
} else {
  pcawg_meta <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
}


pcawg_lung_meta <- pcawg_meta[pcawg_meta$organ_system == "LUNG & BRONCHUS",]
pcawg_lung_meta <- pcawg_lung_meta[pcawg_lung_meta$icgc_donor_id != "DO23717",]



# 
# sig_cont_pcawg <- list()
# sbs_mut_context_mat_pcawg <- matrix(nrow = nrow(pcawg_lung_meta), ncol = nrow(sbs_sig_profile_lung))
# rownames(sbs_mut_context_mat_pcawg) <- pcawg_lung_meta$icgc_donor_id
# colnames(sbs_mut_context_mat_pcawg) <- rownames(sbs_sig_profile_lung)
# 
# dbs_mut_context_mat_pcawg <- matrix(nrow = nrow(pcawg_lung_meta), ncol = nrow(dbs_sig_profile_lung))
# rownames(dbs_mut_context_mat_pcawg) <- pcawg_lung_meta$icgc_donor_id
# colnames(dbs_mut_context_mat_pcawg) <- rownames(dbs_sig_profile_lung)
# 
# id_mut_context_mat_pcawg <- matrix(nrow = nrow(pcawg_lung_meta), ncol = nrow(id_sig_profile_lung))
# rownames(id_mut_context_mat_pcawg) <- pcawg_lung_meta$icgc_donor_id
# colnames(id_mut_context_mat_pcawg) <- rownames(id_sig_profile_lung)
# 
# 
# 
# for (i in 1:nrow(pcawg_lung_meta)){
#   print(i)
#   print(pcawg_lung_meta$icgc_donor_id[i])
#   sample_id <- pcawg_lung_meta$icgc_donor_id[i]
#   
#   
#   if (dir.exists("/hpc/cuppen/")){
#     path <- paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                         sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz")
#   } else {
#   
#   path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                  sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz")
#   }
#   
#     snv_mut_matrix <- extractSigsSnv(path, output='contexts', vcf.filter='PASS', sample.name=sample_id)
#     sbs_mut_context_mat_pcawg[sample_id,] <- as.numeric(snv_mut_matrix)
#     
#     snv_mut_vector <- as.vector(snv_mut_matrix)
#     names(snv_mut_vector) <- rownames(snv_mut_matrix)
#     
#     snv_sig_cont <- fitToSignatures(snv_mut_vector, sbs_sig_profile_lung)
#     
#     
#     
#     dbs_mut_matrix <- extractSigsDbs(path, output='contexts', vcf.filter='PASS', sample.name=sample_id)
#     dbs_mut_context_mat_pcawg[sample_id,] <- as.numeric(dbs_mut_matrix)
#     
#     dbs_mut_vector <- as.vector(dbs_mut_matrix)
#     names(dbs_mut_vector) <- rownames(dbs_mut_matrix)
#     dbs_sig_cont <- fitToSignatures(dbs_mut_vector, dbs_sig_profile_lung)
#     
#    
#      # We will use the PCAWG set of ID mutation types which has 83 elements
#     indel_mut_matrix <- extractSigsIndel(path, output='contexts', vcf.filter='PASS', method = "PCAWG", sample.name=sample_id)
#     id_mut_context_mat_pcawg[sample_id,] <- as.numeric(indel_mut_matrix)
#     
#     indel_mut_vector <- as.vector(indel_mut_matrix)
#     names(indel_mut_vector) <- rownames(indel_mut_matrix)
#     indel_sig_cont <- fitToSignatures(indel_mut_vector, id_sig_profile_lung)
# 
#     
#     sig_cont_pcawg[[sample_id]] <- list(SBS = snv_sig_cont, DBS = dbs_sig_cont, ID = indel_sig_cont)
# }
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(sig_cont_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/list-sig-cont-lung-pcawg.rds")
#   saveRDS(sbs_mut_context_mat_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/sbs_matrix-sig-context-pcawg.rds")
#   saveRDS(dbs_mut_context_mat_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/dbs_matrix-sig-context-pcawg.rds")
#   saveRDS(id_mut_context_mat_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/id_matrix-sig-context-pcawg.rds")
# } else {
#   saveRDS(sig_cont_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/list-sig-cont-lung-pcawg.rds")
#   saveRDS(sbs_mut_context_mat_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/sbs_matrix-sig-context-pcawg.rds")
#   saveRDS(dbs_mut_context_mat_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/dbs_matrix-sig-context-pcawg.rds")
#   saveRDS(id_mut_context_mat_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/id_matrix-sig-context-pcawg.rds")
# } 
# 




if (dir.exists("/hpc/cuppen/")){
  sig_cont_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/list-sig-cont-lung-pcawg.rds")
  sbs_mut_context_mat_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/sbs_matrix-sig-context-pcawg.rds")
  dbs_mut_context_mat_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/dbs_matrix-sig-context-pcawg.rds")
  id_mut_context_mat_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/id_matrix-sig-context-pcawg.rds")
  
  } else {
  sig_cont_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/list-sig-cont-lung-pcawg.rds")
  sbs_mut_context_mat_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/sbs_matrix-sig-context-pcawg.rds")
  dbs_mut_context_mat_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/dbs_matrix-sig-context-pcawg.rds")
  id_mut_context_mat_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/id_matrix-sig-context-pcawg.rds")
  
}


sbs_sig_cont_pcawg_df <- as.data.frame(matrix(unlist(lapply(sig_cont_pcawg, function(x) x[["SBS"]])), ncol = length(lapply(sig_cont_pcawg, function(x) x[["SBS"]])[[1]]), byrow = T))
colnames(sbs_sig_cont_pcawg_df) <- names(lapply(sig_cont_pcawg, function(x) x[["SBS"]])[[1]])
rownames(sbs_sig_cont_pcawg_df) <- names(sig_cont_pcawg)

dbs_sig_cont_pcawg_df <- as.data.frame(matrix(unlist(lapply(sig_cont_pcawg, function(x) x[["DBS"]])), ncol = length(lapply(sig_cont_pcawg, function(x) x[["DBS"]])[[1]]), byrow = T))
colnames(dbs_sig_cont_pcawg_df) <- names(lapply(sig_cont_pcawg, function(x) x[["DBS"]])[[1]])
rownames(dbs_sig_cont_pcawg_df) <- names(sig_cont_pcawg)

id_sig_cont_pcawg_df <- as.data.frame(matrix(unlist(lapply(sig_cont_pcawg, function(x) x[["ID"]])), ncol = length(lapply(sig_cont_pcawg, function(x) x[["ID"]])[[1]]), byrow = T))
colnames(id_sig_cont_pcawg_df) <- names(lapply(sig_cont_pcawg, function(x) x[["ID"]])[[1]])
rownames(id_sig_cont_pcawg_df) <- names(sig_cont_pcawg)




# modified functions for plotting from MutPat package


my_plot_contribution <- function(contribution,
                              signatures = NA,
                              index = NA,
                              coord_flip = FALSE,
                              mode = c("relative", "absolute"),
                              palette = NA) {
  # check mode parameter
  mode <- match.arg(mode)
  
  # optional subsetting if index parameter is provided
  if (!is.na(index)) {
    contribution <- contribution[, index, drop = FALSE]
  }
  
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  Sample <- Contribution <- Signature <- NULL
  
  # When working on NMF results, the contribution needs to be multiplied by the signature colSums.
  if (mode == "absolute" & !is.na(signatures)) {
    # calculate signature contribution in absolute number of signatures
    total_signatures <- colSums(signatures)
    abs_contribution <- contribution * total_signatures
  }
  
  # Make data long. Also create factors for ordering.
  tb <- contribution %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Signature") %>%
    tidyr::pivot_longer(-Signature, names_to = "Sample", values_to = "Contribution") %>%
    dplyr::mutate(
      Sample = factor(Sample, levels = unique(Sample)),
      Signature = factor(Signature, levels = unique(Signature))
    )
  
  # Different plotting between absolute and relative
  if (mode == "absolute") {
    bar_geom <- geom_bar(stat = "identity", colour = "black")
    y_lab <- "Absolute contribution \n (no. mutations)"
  } else if (mode == "relative") {
    bar_geom <- geom_bar(position = "fill", stat = "identity", colour = "black")
    y_lab <- "Relative contribution"
  }
  
  # Determine what signatures are present for the legend.
  present_sigs <- tb %>% 
    dplyr::filter(Contribution != 0) %>% 
    dplyr::pull(Signature) %>% 
    unique()
  
  #Create plot
  plot <- ggplot(tb, aes(x = Sample, y = Contribution, fill = Signature)) +
    bar_geom +
    labs(x = "", y = y_lab) +
    scale_fill_discrete(breaks = present_sigs) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    )
    
  # Allow custom color palettes.
  if (!is.na(palette)) {
    plot <- plot + scale_fill_manual(name = "Signature", values = palette)
  }
  
  # Handle coord_flip.
  if (coord_flip) {
    plot <- plot + coord_flip() + xlim(rev(levels(factor(tb$Sample)))) + theme(axis.ticks.y = element_blank()) + theme(axis.text.y = element_blank())
    
  } else {
    plot <- plot + xlim(levels(factor(tb$Sample))) + theme(axis.ticks.x = element_blank()) + theme(axis.text.x = element_blank())
    
  }
  
  return(plot)
}



my_plot_contribution_heatmap <- function(contribution, sig_order = NA, sample_order = NA, cluster_samples = TRUE,
                                      cluster_sigs = FALSE, method = "complete", plot_values = FALSE) {
  
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  Signature <- Sample <- Contribution <- x <- y <- xend <- yend <- NULL
  
  # check contribution argument
  if (!inherits(contribution, "matrix")) {
    stop("contribution must be a matrix")
  }
  # check if there are signatures names in the contribution matrix
  if (is.null(row.names(contribution))) {
    stop("contribution must have row.names (signature names)")
  }
  
  # transpose
  contribution <- t(contribution)
  # relative contribution
  contribution_norm <- contribution / rowSums(contribution)
  
  # If cluster_samples is TRUE perform clustering. Else use supplied sample_order or
  # the current column order.
  if (!is.na(sample_order) & cluster_samples == TRUE) {
    stop("sample_order can only be provided when cluster_samples is FALSE", call. = FALSE)
  } else if (!is.na(sample_order)) {
    # check sample_order argument
    if (!inherits(sample_order, "character")) {
      stop("sample_order must be a character vector", call. = FALSE)
    }
    if (length(sample_order) != nrow(contribution_norm)) {
      stop("sample_order must have the same length as the number
          of samples in the explained matrix", call. = FALSE)
    }
  } else if (cluster_samples == TRUE) {
    # cluster samples based on eucledian distance between relative contribution_norm
    hc.sample <- hclust(dist(contribution_norm), method = method)
    # order samples according to clustering
    sample_order <- rownames(contribution_norm)[hc.sample$order]
    
    dhc <- as.dendrogram(hc.sample)
    # rectangular lines
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_y_reverse(expand = c(0.2, 0)) +
      ggdendro::theme_dendro()
  }
  else {
    sample_order <- rownames(contribution_norm)
  }
  
  
  # If cluster_sigs is TRUE perform clustering. Else use supplied sig_order or
  # the current column order.
  if (!is.na(sig_order) & cluster_sigs == TRUE) {
    stop("sig_order can only be provided when cluster_sigs is FALSE", call. = FALSE)
  } else if (!is.na(sig_order)) {
    # check sig_order argument
    if (!inherits(sig_order, "character")) {
      stop("sig_order must be a character vector", call. = FALSE)
    }
    if (length(sig_order) != ncol(contribution_norm)) {
      stop("sig_order must have the same length as the number
           of signatures in the explained matrix", call. = FALSE)
    }
  } else if (cluster_sigs == TRUE) {
    # Cluster cols
    hc.sample2 <- contribution_norm %>%
      t() %>%
      dist() %>%
      hclust(method = method)
    sig_order <- colnames(contribution_norm)[hc.sample2$order]
    
    dhc <- as.dendrogram(hc.sample2)
    # rectangular lines
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_cols <- ggplot(ggdendro::segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      ggdendro::theme_dendro() +
      scale_y_continuous(expand = c(0.2, 0))
  } else {
    sig_order <- colnames(contribution_norm)
  }
  
  # Make matrix long and set factor levels, to get the correct order for plotting.
  contribution_norm.m <- contribution_norm %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "Signature", values_to = "Contribution") %>%
    dplyr::mutate(
      Signature = factor(Signature, levels = sig_order),
      Sample = factor(Sample, levels = sample_order)
    )
  
  # plot heatmap
  heatmap <- ggplot(contribution_norm.m, aes(x = Signature, y = Sample, fill = Contribution, order = Sample)) +
    geom_raster() +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Relative \ncontribution", limits = c(0, 1)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Signatures", y = "PCAWG Lung Samples")
  # if plot_values is TRUE, add values to heatmap
  if (plot_values) {
    heatmap <- heatmap + geom_text(aes(label = round(Contribution, 2)), size = 3)
  }
  
  return(heatmap)
}

# plot_96_profile

# example
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
                               package = "MutationalPatterns"
))


sbs_mut_context_mat_pcawg <- t(sbs_mut_context_mat_pcawg)
arr <- plot_96_profile(sbs_mut_context_mat_pcawg[,1:5], colors = NA, ymax = 0.2, condensed = FALSE)

# compare two pofiles! Could be handy to compare the average of primary and metastatic profiles

plot_compare_profiles(profile1 = sbs_mut_context_mat_pcawg[,1], profile2 = sbs_mut_context_mat_pcawg[,2], profile_names = c("a", "b"))




# plot contributions

nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
                               package = "MutationalPatterns"
))


# I need a matrix where rows are signature names and columns are sample names.


sbs_sig_cont_pcawg_mat <- as.matrix(sbs_sig_cont_pcawg_df)
sbs_sig_cont_pcawg_mat <- t(sbs_sig_cont_pcawg_mat)

dbs_sig_cont_pcawg_mat <- as.matrix(dbs_sig_cont_pcawg_df)
dbs_sig_cont_pcawg_mat <- t(dbs_sig_cont_pcawg_mat)

id_sig_cont_pcawg_mat <- as.matrix(id_sig_cont_pcawg_df)
id_sig_cont_pcawg_mat <- t(id_sig_cont_pcawg_mat)






# Unsupervised clustering and plotting the dendogram


#### SBS


sbs_sig_cont_pcawg_mat_norm <- sbs_sig_cont_pcawg_mat / rowSums(sbs_sig_cont_pcawg_mat)
hc.sample <- hclust(dist(t(sbs_sig_cont_pcawg_mat_norm)), method = "complete")
sample_order <- colnames(sbs_sig_cont_pcawg_mat_norm)[hc.sample$order]
sbs_sig_cont_pcawg_mat <- sbs_sig_cont_pcawg_mat[,sample_order]

cols <- append(brewer.pal(n = 12, name = "Paired"), brewer.pal(n = 3, name = "Dark2"))
sbs_stacked_bar_cont <- my_plot_contribution(sbs_sig_cont_pcawg_mat, signatures = NA,
                  index = NA,
                  coord_flip = T,
                  mode = "relative",
                  palette = cols)

sbs_stacked_bar_cont



dhc <- as.dendrogram(hc.sample)

ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
# plot dendrogram of hierachical clustering
dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdendro::theme_dendro() +
  coord_flip() +
  scale_y_reverse(expand = c(0.1, 0)) +
  scale_x_reverse(expand = c(0,0.8)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


empty_fig <- ggplot() +
  theme_void()
sbs_stacked_bar_cont_clustering <- cowplot::plot_grid(empty_fig, dendrogram_rows, sbs_stacked_bar_cont, align = "hv", rel_widths = c(1, 0.5, 5), nrow = 1)



# plot contributions in a heatmap format


sbs_heatmap_count <- my_plot_contribution_heatmap(sbs_sig_cont_pcawg_mat, sig_order = NA, sample_order = NA, cluster_samples = F,
                          cluster_sigs = FALSE, method = "complete", plot_values = FALSE)


## Saving




for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/png/sbs-rel-cont-stacked-barplot.png", width = 960, height = 960)
    print(sbs_stacked_bar_cont_clustering)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/png/sbs-rel-cont-stacked-heatmap.png", width = 660, height = 660)
    print(sbs_heatmap_count +
            ggtitle("Single-base Substitution Contribution Heatmap \n \n") +
            theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) 
    )
    dev.off()
    
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/pdf/sbs-rel-cont-stacked-barplot.pdf", width = 14, height = 14)
    print(sbs_stacked_bar_cont_clustering)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/pdf/sbs-rel-cont-stacked-heatmap.pdf", width = 10, height = 10)
    print(sbs_heatmap_count +
            ggtitle("Single-base Substitution Contribution Heatmap \n \n") +
            theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) 
    )
    dev.off()
  }
}



#### DBS




dbs_sig_cont_pcawg_mat_norm <- dbs_sig_cont_pcawg_mat / rowSums(dbs_sig_cont_pcawg_mat)
hc.sample <- hclust(dist(t(dbs_sig_cont_pcawg_mat_norm)), method = "complete")
sample_order <- colnames(dbs_sig_cont_pcawg_mat_norm)[hc.sample$order]
dbs_sig_cont_pcawg_mat <- dbs_sig_cont_pcawg_mat[,sample_order]

cols <- brewer.pal(n = 3, name = "Paired")
dbs_stacked_bar_cont <- my_plot_contribution(dbs_sig_cont_pcawg_mat[,colSums(dbs_sig_cont_pcawg_mat) != 0], signatures = NA,
                  index = NA,
                  coord_flip = T,
                  mode = "relative",
                  palette = cols)




dhc <- as.dendrogram(hc.sample)

ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
# plot dendrogram of hierachical clustering
dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdendro::theme_dendro() +
  coord_flip() +
  scale_y_reverse(expand = c(0.1, 0)) +
  scale_x_reverse(expand = c(0,0.8)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


empty_fig <- ggplot() +
  theme_void()
dbs_stacked_bar_cont_clustering <- cowplot::plot_grid(empty_fig, dendrogram_rows, dbs_stacked_bar_cont, align = "hv", rel_widths = c(1, 0.5, 5), nrow = 1)



# plot contributions in a heatmap format
dbs_sig_cont_pcawg_mat <- dbs_sig_cont_pcawg_mat[,colSums(dbs_sig_cont_pcawg_mat) != 0]

dbs_heatmap_count <- my_plot_contribution_heatmap(dbs_sig_cont_pcawg_mat, sig_order = NA, sample_order = NA, cluster_samples = F,
                          cluster_sigs = FALSE, method = "complete", plot_values = FALSE)





## Saving



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/png/dbs-rel-cont-stacked-barplot.png", width = 960, height = 960)
    print(dbs_stacked_bar_cont_clustering)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/png/dbs-rel-cont-stacked-heatmap.png", width = 660, height = 660)
    print(dbs_heatmap_count +
            ggtitle("Double-base Substitution Contribution Heatmap \n \n") +
            theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) 
    )
    dev.off()
    
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/pdf/dbs-rel-cont-stacked-barplot.pdf", width = 14, height = 14)
    print(dbs_stacked_bar_cont_clustering)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/pdf/dbs-rel-cont-stacked-heatmap.pdf", width = 10, height = 10)
    print(dbs_heatmap_count +
            ggtitle("Double-base Substitution Contribution Heatmap \n \n") +
            theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) 
    )
    dev.off()
  }
}




#### ID



id_sig_cont_pcawg_mat_norm <- id_sig_cont_pcawg_mat / rowSums(id_sig_cont_pcawg_mat)
hc.sample <- hclust(dist(t(id_sig_cont_pcawg_mat_norm)), method = "complete")
sample_order <- colnames(id_sig_cont_pcawg_mat_norm)[hc.sample$order]
id_sig_cont_pcawg_mat <- id_sig_cont_pcawg_mat[,sample_order]

cols <- brewer.pal(n = 6, name = "Dark2")
id_stacked_bar_cont <- my_plot_contribution(id_sig_cont_pcawg_mat[,colSums(id_sig_cont_pcawg_mat) != 0], signatures = NA,
                  index = NA,
                  coord_flip = T,
                  mode = "relative",
                  palette = cols)

id_sig_cont_pcawg_mat[,"DO23648"]


dhc <- as.dendrogram(hc.sample)

ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
# plot dendrogram of hierachical clustering
dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdendro::theme_dendro() +
  coord_flip() +
  scale_y_reverse(expand = c(0.1, 0)) +
  scale_x_reverse(expand = c(0,0.8)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


empty_fig <- ggplot() +
  theme_void()
id_stacked_bar_cont_clustering <- cowplot::plot_grid(empty_fig, dendrogram_rows, id_stacked_bar_cont, align = "hv", rel_widths = c(1, 0.5, 5), nrow = 1)



# plot contributions in a heatmap format
id_sig_cont_pcawg_mat <- id_sig_cont_pcawg_mat[,colSums(id_sig_cont_pcawg_mat) != 0]

id_heatmap_count <- my_plot_contribution_heatmap(id_sig_cont_pcawg_mat, sig_order = NA, sample_order = NA, cluster_samples = TRUE,
                          cluster_sigs = FALSE, method = "ward.D2", plot_values = FALSE)





## Saving




for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/png/id-rel-cont-stacked-barplot.png", width = 960, height = 960)
    print(id_stacked_bar_cont_clustering)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/png/id-rel-cont-stacked-heatmap.png", width = 660, height = 660)
    print(id_heatmap_count +
            ggtitle("Small INDEL Contribution Heatmap \n \n") +
            theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) 
    )
    dev.off()
    
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/pdf/id-rel-cont-stacked-barplot.pdf", width = 14, height = 14)
    print(id_stacked_bar_cont_clustering)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/pdf/id-rel-cont-stacked-heatmap.pdf", width = 10, height = 10)
    print(id_heatmap_count +
            ggtitle("Small INDEL Contribution Heatmap \n \n") +
            theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) 
    )
    dev.off()
  }
}




# Combine


combined_sig_cont_norm <- rbind(sbs_sig_cont_pcawg_mat_norm, dbs_sig_cont_pcawg_mat_norm, id_sig_cont_pcawg_mat_norm)
hc.sample.combined <- hclust(dist(t(combined_sig_cont_norm)), method = "complete")
sample_order_comb_sbs <- colnames(sbs_sig_cont_pcawg_mat_norm)[hc.sample.combined$order]

sbs_sig_cont_pcawg_mat <- sbs_sig_cont_pcawg_mat[,sample_order_comb_sbs]
dbs_sig_cont_pcawg_mat <- dbs_sig_cont_pcawg_mat[,sample_order_comb_sbs]
id_sig_cont_pcawg_mat <- id_sig_cont_pcawg_mat[,sample_order_comb_sbs]




cols <- append(brewer.pal(n = 12, name = "Paired"), brewer.pal(n = 3, name = "Dark2"))
sbs_stacked_bar_cont <- my_plot_contribution(sbs_sig_cont_pcawg_mat, signatures = NA,
                                             index = NA,
                                             coord_flip = F,
                                             mode = "relative",
                                             palette = cols)





cols <- brewer.pal(n = 3, name = "Paired")
dbs_stacked_bar_cont <- my_plot_contribution(dbs_sig_cont_pcawg_mat[,colSums(dbs_sig_cont_pcawg_mat) != 0], signatures = NA,
                                             index = NA,
                                             coord_flip = F,
                                             mode = "relative",
                                             palette = cols)





cols <- brewer.pal(n = 6, name = "Dark2")
id_stacked_bar_cont <- my_plot_contribution(id_sig_cont_pcawg_mat[,colSums(id_sig_cont_pcawg_mat) != 0], signatures = NA,
                                            index = NA,
                                            coord_flip = F,
                                            mode = "relative",
                                            palette = cols)




# Adding dendogram     # Failed to align them perfectly! 

dhc <- as.dendrogram(hc.sample.combined)

ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
# plot dendrogram of hierachical clustering

ggdendrogram(segment(ddata))

dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdendro::theme_dendro()


dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdendro::theme_dendro() +
  coord_flip() +
  scale_y_reverse(expand = c(0.1, 0)) +
  scale_x_reverse(expand = c(0,0.8)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

empty_fig <- ggplot() +
  theme_void()
sbs_stacked_bar_cont_clustering <- cowplot::plot_grid(dendrogram_rows, sbs_stacked_bar_cont, align = "hv", rel_heights = c(0.5,2), nrow = 2, ncol = 1, scale=c(0.675,1))




pcawg_mut_sig_spectrum <- cowplot::plot_grid(sbs_stacked_bar_cont, dbs_stacked_bar_cont, id_stacked_bar_cont, align = "v", nrow = 3, ncol = 1)




for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/png/combined-rel-cont-stacked-barplots.png", width = 960, height = 960)
    print(pcawg_mut_sig_spectrum)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/figs/pdf/combined-rel-cont-stacked-barplots.pdf", width = 14, height = 14)
    print(pcawg_mut_sig_spectrum)
    dev.off()
  }
}
