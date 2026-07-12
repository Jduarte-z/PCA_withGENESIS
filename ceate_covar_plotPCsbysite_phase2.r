pcair_r1 <- readRDS("mypcair_1round_results_rds_1stRound")
pcs <- as.data.frame(pcair_r1$vectors)
str(pcair_r1)

colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
pcs <- cbind(sample.id = rownames(pcs), pcs)
write.table(pcs, "pcair_r1_covariates.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

pcs <- read.table("pcair_r1_covariates.tsv", sep = "\t", header = TRUE)
head(pcs)

psam <- read.table("/home/duartej3/beegfs/JF/LPD_DATA/phase2_newQC/qc_new_samplesUpdated/update_psam/phase2.psam",
                   sep = "\t", header = TRUE, comment.char = "", check.names = FALSE)
colnames(psam)[1] <- "IID"

head(psam)

# check IDs match between the two files
# how many psam samples have PCs
sum(psam$IID %in% pcs$sample.id)   
nrow(psam)  

merged <- merge(psam, pcs, by.x = "IID", by.y = "sample.id", all.x = TRUE)
head(merged)
colSums(is.na(merged))
colMeans(is.na(merged))

library(ggplot2)
library(cowplot)

pairs <- list(c("PC1","PC2"), c("PC3","PC4"), c("PC5","PC6"), c("PC7","PC8"), c("PC9","PC10"))

#by site
plots <- lapply(pairs, function(p) {
  ggplot(merged, aes_string(x = p[1], y = p[2], color = "SITE")) +
    geom_point(alpha = 0.5, size = 0.8) +
    theme_bw() +
    theme(legend.position = "none")
})

legend <- get_legend(
  ggplot(merged, aes(x = PC1, y = PC2, color = SITE)) +
    geom_point() + theme_bw()
)

final_plot <- plot_grid(plot_grid(plotlist = plots, ncol = 3), legend, rel_widths = c(1, 0.15))

ggsave("pcair_r1_bySite.png", final_plot, width = 15, height = 8, dpi = 300)

#by country 

plots <- lapply(pairs, function(p) {
  ggplot(merged, aes_string(x = p[1], y = p[2], color = "COUNTRY")) +
    geom_point(alpha = 0.5, size = 0.8) +
    theme_bw() +
    theme(legend.position = "none")
})

legend <- get_legend(
  ggplot(merged, aes(x = PC1, y = PC2, color = COUNTRY)) +
    geom_point() + theme_bw()
)

final_plot <- plot_grid(plot_grid(plotlist = plots, ncol = 3), legend, rel_widths = c(1, 0.15))

ggsave("pcair_r1_byCountry.png", final_plot, width = 15, height = 8, dpi = 300)



