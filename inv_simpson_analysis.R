metadata <- data.frame(
  Sample = c(
    "SRR3721799_clones_TRD",
    "SRR3721799_clones_TRA",
    "SRR3721797_clones_TRB",
    "SRR3721797_clones_TRA",
    "SRR3721796_clones_TRD",
    "SRR3721796_clones_TRA"
  )
)

# Save the metadata.csv file to the 'clones' subfolder relative to your project directory
write.csv(metadata, "clones/metadata.csv", row.names = FALSE)

# Load the immune data from the 'clones' folder
library(immunarch)

immdata <- repLoad("clones/")


head(immdata$meta)

# View the data for a specific sample
head(immdata$data$SRR3721799_clones_TRA)

# Clonotype count comparison
repExplore(immdata$data, .method = "volume")


# Overlap between samples
repOverlap(immdata$data, .method = "public")



# Clone assembly heatmap across samples
overlap_heatmap<- repOverlap(immdata$data, .method = "public", .verbose = TRUE) %>%
  vis()

#save heatmap
ggsave('overlap_heatmap.pdf',overlap_heatmap,width = 10, height = 8, dpi = 300)

# Calculate Inverse Simpson diversity
inv_simp <- repDiversity(immdata$data,.method="inv.simp")

write.table(inv_simp,'InverseSimpson.tsv')