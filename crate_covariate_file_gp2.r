pcair_r1 <- readRDS("mypcair_r2_results_rds")
pcs <- as.data.frame(pcair_r1$vectors)
str(pcair_r1)

colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
pcs <- cbind(sample.id = rownames(pcs), pcs)
head(pcs)

names(pcs)[names(pcs) == 'sample.id'] <- 'IID'
write.table(pcs, "pcair_r2_covariates.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

pcs <- read.table("pcair_r2_covariates.tsv", sep = "\t", header = TRUE)
head(pcs)

#to plot pcs by site or by specific covariates stored in the psam file, not applicable here


psam <- read.table("../../EURdownSampled_noMatchAge.covar.tsv",
                   sep = "\t", header = TRUE, check.names = FALSE)
colnames(psam)[1] <- "IID"

head(psam)

# check IDs match between the two files
sum(psam$IID %in% pcs$IID)   # how many psam samples have PCs
nrow(psam)  

merged <- merge(psam, pcs, by.x = "IID", by.y = "IID", all.x = TRUE)
head(merged)
colSums(is.na(merged))
colMeans(is.na(merged))

write.table(merged[, !names(merged) %in% c("AAO", "AAD")],
            "pcair_r2_covariates_merged.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


#the following not needed for now:

# library(ggplot2)
# library(cowplot)

# pairs <- list(c("PC1","PC2"), c("PC3","PC4"), c("PC5","PC6"), c("PC7","PC8"), c("PC9","PC10"))

# #by site
# plots <- lapply(pairs, function(p) {
#   ggplot(merged, aes_string(x = p[1], y = p[2], color = "SITE")) +
#     geom_point(alpha = 0.5, size = 0.8) +
#     theme_bw() +
#     theme(legend.position = "none")
# })

# legend <- get_legend(
#   ggplot(merged, aes(x = PC1, y = PC2, color = SITE)) +
#     geom_point() + theme_bw()
# )

# final_plot <- plot_grid(plot_grid(plotlist = plots, ncol = 3), legend, rel_widths = c(1, 0.15))

# ggsave("pcair_r2_bySite.png", final_plot, width = 15, height = 8, dpi = 300)

# #by country 

# plots <- lapply(pairs, function(p) {
#   ggplot(merged, aes_string(x = p[1], y = p[2], color = "COUNTRY")) +
#     geom_point(alpha = 0.5, size = 0.8) +
#     theme_bw() +
#     theme(legend.position = "none")
# })

# legend <- get_legend(
#   ggplot(merged, aes(x = PC1, y = PC2, color = COUNTRY)) +
#     geom_point() + theme_bw()
# )

# final_plot <- plot_grid(plot_grid(plotlist = plots, ncol = 3), legend, rel_widths = c(1, 0.15))

# ggsave("pcair_r2_byCountry.png", final_plot, width = 15, height = 8, dpi = 300)



