#!/usr/bin/env Rscript
# This is an R script for visualizing the relationships between FAPROTAX2 and KEGG functional traits detected
# in microbiome samples.
#
#
#
# Rob Porch
# Updated 2026 Mar 25
#

##############################
# Load necessary packages
library(tidyverse)
library(RColorBrewer)
library(pals)
library(data.table)
library(readxl)
library(optparse)
library(showtext)
library(extrafont)
library(data.table)
###############################

# Parse command-line arguments
option_list <- list(
  make_option(c("-a", "--in_asv_table"), type="character", help="Path to input ASV table."),
  make_option(c("-t", "--in_taxa_table"), type="character", help="Path to input taxa table."),
  make_option(c("-f", "--in_function_table"), type="character", help="Path to FAPROTAX2-mapper output function table."),
  make_option(c("-s", "--in_sample_metadata"), type="character", help="Path to microbiome sample metadata."),
  make_option(c("-o", "--out_figure"), type="character", help="Path to store output stacked barplot showing the functional distribution for each microbiome sample.")
)

opt <- parse_args(OptionParser(option_list=option_list))


#######################
# MAIN SCRIPT BODY
#######################

asv_table <- data.table::fread(file=opt$in_asv_table, header=TRUE)
taxa_table <- data.table::fread(file=opt$in_taxa_table, header=TRUE)
names(taxa_table) <- c("ASV", "taxonomy")

function_table <- data.table::fread(file=opt$in_function_table, header=TRUE)

sample_metadata <- readxl::read_xlsx(opt$in_sample_metadata)
names(sample_metadata) <- c("SampleID", "sample_type")
sample_metadata$sample_type <- gsub("_", " ", sample_metadata$sample_type)

sample_metadata$sample_type <- gsub("Freshwater and sediment", "Freshwater\nand sediment",
									sample_metadata$sample_type)
									
sample_metadata$sample_type <- gsub("Hypersaline lagoon sediment", "Hypersaline\nlagoon\nsediment",
									sample_metadata$sample_type)
									
sample_metadata$sample_type <- gsub("Marine water and sediment", "Marine water\nand sediment",
									sample_metadata$sample_type)

func_ids <- function_table[[1]]
func_data <- function_table[, -1]

func_df <- cbind(`function` = func_ids, func_data)
superfunctions <- c("chemoheterotrophy", "phototrophy", "photoautotrophy")
func_df <- func_df[!`function` %in% superfunctions]

# Calculate relative abundances of functions after filtering
func_rel <- func_df[, lapply(.SD, function(x) x / sum(x, na.rm=T)), .SDcols = -1]
func_rel_with_ids <- cbind(`function` = func_df$`function`, func_rel)


### Convert dataframes to long format
func_long <- func_rel_with_ids %>%
  pivot_longer(cols = -`function`, names_to="SampleID", values_to="relative_abundance")

func_long$`function` <- gsub("_", " ", func_long$`function`)

### Create large color palette for functions
# base_function_palette <- c(
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# )
base_function_palette <- as.vector(polychrome(31))
function_palette <- rep_len(base_function_palette, 31)

func_long <- func_long %>%
  left_join(sample_metadata, by="SampleID")

top_30_traits <- func_long %>%
  group_by(`function`) %>%
  summarize(mean_abund = mean(relative_abundance)) %>%
  slice_max(mean_abund, n = 30) %>%
  pull(`function`)

func_order <- func_long %>%
  group_by(`function`) %>%
  summarize(mean_abund = mean(relative_abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(`function`)

func_plot_data <- func_long %>%
  mutate(func_group = if_else(`function` %in% top_30_traits, `function`, "other"))

func_order_grouped <- func_plot_data %>%
	group_by(func_group) %>%
	summarize(mean_abund = mean(relative_abundance)) %>%
	arrange(mean_abund) %>%
	pull(func_group)

func_plot_data$func_group <- factor(func_plot_data$func_group, levels = func_order_grouped)

# generate colors for non-"other" functions
n_funcs <- length(func_order_grouped) - 1  # exclude "other"
set.seed(145)
base_colors <- sample(rep_len(base_function_palette, n_funcs))

# name them according to function order (excluding "other")
color_mapping <- setNames(base_colors, func_order_grouped[func_order_grouped != "other"])

# add "other" manually
color_mapping["other"] <- "grey70"

function_stacked_barplot <- ggplot(func_plot_data, 
									aes(x = SampleID, y = relative_abundance, 
										fill = func_group)) +
  geom_bar(stat="identity", width = 1, position=position_stack()) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values=color_mapping) + facet_grid(~sample_type, scales="free_x", space="free_x", switch="x",
                                                     labeller=(sample_type = label_wrap_gen(width=20))) +
  scale_x_discrete(expand = expansion(add = 0.1)) + coord_cartesian(clip="off") +
  theme_classic() +
  labs(
    x = "",
    y = "Relative abundance",
    fill = "function"
  ) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.line.y = element_blank(),
    legend.text = element_text(size=6),
    legend.title = element_text(size=8),
    legend.key.size = unit(0.35, "cm"),
    strip.clip="off",
    
    # Style facet grids
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(angle=45, hjust=1, vjust=1, size=6, margin = margin(b = 2)),
    strip.switch.pad.grid = unit(0.05, "cm"),
    panel.spacing.x = unit(0.25, "cm")
  )

ggsave(opt$out_figure, width=8, height=3)
