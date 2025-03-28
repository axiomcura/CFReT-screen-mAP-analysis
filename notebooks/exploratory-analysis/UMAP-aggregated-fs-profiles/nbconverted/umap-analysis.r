suppressMessages(library(umap))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggthemes))

# set file path to the concatenated aggregated feature selected profiles
path.data = file.path("./results/concat_data/batch_1_concat_agg_fs.csv")
if (!file.exists(path.data)){
    stop("Concatenated aggregated feature selected profiles not found. Please run the 'preprocessing-data.ipynb' script first.")
}

# setting output path
path.output = file.path("./results/umap")
if (!dir.exists(path.output)){
    dir.create(path.output, recursive = TRUE)
}


# loading the data
concat.agg_fs_df = read.csv(path.data)
head(concat.agg_fs_df)

# setting metadata to keep
sel_metadata <- c("Metadata_plate_barcode", "Metadata_control_type", "Metadata_Pathway", "Metadata_treatment")

# separating metadata and morphology data
metadata_df <- concat.agg_fs_df[,sel_metadata]
morphology_df <- concat.agg_fs_df[,!grepl("^Metadata_", colnames(concat.agg_fs_df))]

# if Treatment is DMSO and Metadata_Pathway is empty, set Metadata_Pathway to DMSO
metadata_df$Metadata_Pathway[metadata_df$Metadata_treatment == "DMSO-positive" & metadata_df$Metadata_Pathway == ""] <- "DMSO"
metadata_df$Metadata_Pathway[metadata_df$Metadata_treatment == "DMSO-negative" & metadata_df$Metadata_Pathway == ""] <- "DMSO"

# if Treatment is not DMSO and Metadata_Pathway is empty, set Metadata_Pathway to unknown
metadata_df$Metadata_Pathway[metadata_df$Metadata_treatment != "DMSO" & metadata_df$Metadata_Pathway == ""] <- "Unknown"

# setting seed
set.seed(0)

# generated a control_df where the Metadata_control_type is positive or negative.
# only select rows where the Metadata_control_type is positive or negative
control_df <- concat.agg_fs_df[concat.agg_fs_df$Metadata_control_type %in%
							   c("positive", "negative"), ]
control_morphology_df <- control_df[,!grepl("^Metadata_", colnames(control_df))]
cntrls_umap_result <- umap(control_morphology_df, n_components = 2)
umap_control_df <- data.frame(
  Plate = control_df$Metadata_plate_barcode,
  Pathway = control_df$Metadata_Pathway,
  ControlType = control_df$Metadata_control_type,
  UMAP1 = cntrls_umap_result$layout[,1],
  UMAP2 = cntrls_umap_result$layout[,2]
)

# setting seed
set.seed(0)

# Extract only morphology data (exclude metadata)
morphology_features <- morphology_df

# Ensure row names in metadata match morphology_features
rownames(metadata_df) <- rownames(morphology_features)

# Perform UMAP
umap_result <- umap(morphology_features, n_components = 2)

# Create a data frame with UMAP results and metadata
umap_df <- data.frame(
  Plate = metadata_df$Metadata_plate_barcode,
  Pathway = metadata_df$Metadata_Pathway,
  Treatment = metadata_df$Metadata_treatment,
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2]
)

# make figure larger
options(repr.plot.width=15, repr.plot.height=11)

# Generate a color palette dynamically for ControlType (colorblind-friendly)
control_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(umap_control_df$ControlType))),
  unique(umap_control_df$ControlType)
)

# Define distinct filled shapes for Plate (21–25 allow fill)
plate_shapes <- c(21, 22, 23, 24, 25)  # These allow filling
names(plate_shapes) <- unique(umap_control_df$Plate)

# plotting the UMAP
ggplot(umap_control_df, aes(x = UMAP1, y = UMAP2, shape = Plate, fill = ControlType)) +
  geom_point(size = 4, stroke = 1, color = "black", alpha = 0.9) +  # Larger points with black border
  scale_shape_manual(values = plate_shapes) +  # Assign shapes to Plate
  scale_fill_manual(values = control_colors) +  # Fill shapes by ControlType
  theme_bw(base_size = 18) +  # Increase base font size
  labs(
    title = "UMAP of control profiles",
    x = "UMAP 1",
    y = "UMAP 2",
    fill = "Control type",
    shape = "Plate"
  ) +
  guides(fill = guide_legend(override.aes = list(color = control_colors))) +  # Add color to legend
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16), # Place legend on the right
    legend.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )


# save figure into the output folder
path.figure <- file.path(path.output, "umap_control_profiles.png")
ggsave(path.figure, width = 15, height = 11, units = "in", dpi = 300)
print(paste("UMAP of control profiles saved at", path.figure))

# Generate a color palette dynamically for unique plates
plate_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(umap_df$Plate))),
  unique(umap_df$Plate)
)

# Generate a color palette dynamical for unique treatments (total of 51)
treatment_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(umap_df$Treatment))),
  unique(umap_df$Treatment)
)

# set seed
set.seed(0)

# Plot UMAP with enhanced aesthetics for publication
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Plate)) +
  geom_point(alpha = 0.7, size = 5, shape = 16) +
  scale_color_manual(values = plate_colors) +
  theme_bw(base_size = 18) +
  labs(
    title = "UMAP of profiles across all plates",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Plate"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

  # save figure into the output folder
path.figure <- file.path(path.output, "umap_all_plate_profiles.png")
ggsave(path.figure, width = 15, height = 11, units = "in", dpi = 300)
print(paste("UMAP of all plate profiles saved at", path.figure))


# Generate a color palette dynamically for unique pathways (colorblind-friendly)
pathway_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(umap_df$Pathway))),
  unique(umap_df$Pathway)
)

# Plotting UMAP coloring based on the Pathway
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Pathway)) +
  geom_point(alpha = 0.7, size = 5, shape = 16) +
  scale_color_manual(values = pathway_colors) +
  theme_bw( base_size = 18) +
  labs(
    title = "UMAP of treated profiles and targeted associated pathways",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pathway"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

# save figure into the output folder
path.figure <- file.path(path.output, "umap_pathway_profiles.png")
ggsave(path.figure, width = 15, height = 11, units = "in", dpi = 300)
print(paste("UMAP of pathway profiles saved at", path.figure))


# update figure size
options(repr.plot.width = 15, repr.plot.height = 30)

# using the same UMAP coordinates, create a new dataframe as "background_df" in order to plot the background
background_df <- data.frame(
  Plate = umap_df$Plate,
  UMAP1 = umap_df$UMAP1,
  UMAP2 = umap_df$UMAP2
)

# create a facet grid plot where the rows represent Plate and the columns represent Pathway
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Pathway)) +

  # generate a background plot in light gray
  geom_point(data = background_df, aes(x = UMAP1, y = UMAP2), color = "lightgray", alpha = 0.7, size = 5) +

  # generate forground plot with color based on Pathway
  geom_point(alpha = 0.7, size = 5, shape = 16) +
  scale_color_manual(values = pathway_colors) +
  theme_bw(base_size = 18) +
  labs(
    title = "UMAP of profiles across all plates and pathways",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pathway"
  ) +
  facet_grid(Pathway ~ Plate, scales = "free") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

# save plot
path.figure <- file.path(path.output, "facetplot_umap_all_plate_pathway_profiles.png")
ggsave(path.figure, width = 15, height = 30, units = "in", dpi = 300)
print(paste("UMAP of all plate and pathway profiles saved at", path.figure))


# filter this to only where Pathwat is DMSO
umap_df_dmsos <- umap_df[umap_df$Pathway == "DMSO", ]

umap_df_dmsos
