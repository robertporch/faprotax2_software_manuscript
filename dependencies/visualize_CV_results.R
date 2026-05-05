#### /usr/bin/Rscript
#
#
#
# This is an R script for visualizing the cross-validation results
# between FAPROTSAX v1 and v2.
#
#
# Rob Porch
# last updated 2026 March 11

###################
# Load packages
library(tidyverse)
library(data.table)
library(showtext)
library(optparse)
library(patchwork)
library(extrafont)
####################


####################
# Set up command-line arguments

option_list <- list(
	make_option(c("-s", "--in_sequence_CV_results"), type="character", default="", help="Path to CV results table, using sequence-based placement in FAPROTAX2-mapper."),
	make_option(c("-t", "--in_taxonomy_CV_results"), type="character", default="", help="Path to CV results table, using taxonomy-based placement in FAPROTAX2-mapper"),
	make_option(c("-o", "--out_figure"), type="character", default="", help="Path to store output figure showing CV results."),
	make_option("--out_summary_table", type="character", default="", help="Path to store output summary table, listing average scores separately for each version."),
	make_option("--out_full_table", type="character", default="", help="Path to store output full summary table, listing average scores, separately for each version, and each function.")
)

# Parse command-line arguments
opt <- parse_args(optparse::OptionParser(option_list = option_list))



###################
# Main script body
###################

# Load in CV results tables
sequence_CV_results <- data.table::fread(file=opt$in_sequence_CV_results, header=TRUE)
taxonomy_CV_results <- data.table::fread(file=opt$in_taxonomy_CV_results, header=TRUE)

# Create df for boxplots
sequence_CV_long <- sequence_CV_results %>%
# pivot the mean columns into long format
	pivot_longer(
		cols=starts_with("mean_"),
		names_to=c("Version", "Metric"),
		names_pattern="mean_(v1|v2)_(BA|TPR|TNR)",
		values_to="Score"
	) %>%
	mutate(
		Version=ifelse(Version=="v1", "FAPROTAX v1", "FAPROTAX2-mapper\n(sequence-based)"),
		Metric_Full = case_when(
			Metric == "BA" ~ "Balanced accuracy",
			Metric == "TPR" ~ "True positive rate",
			Metric == "TNR" ~ "True negative rate"
		)
	)
	
taxonomy_CV_long <- taxonomy_CV_results %>%
# pivot the mean columns into long format
	pivot_longer(
		cols=starts_with("mean_"),
		names_to=c("Version", "Metric"),
		names_pattern="mean_(v1|v2)_(BA|TPR|TNR)",
		values_to="Score"
	) %>%
	mutate(
		Version=ifelse(Version=="v1", "FAPROTAX v1", "FAPROTAX2-mapper\n(taxonomy-based)"),
		Metric_Full = case_when(
			Metric == "BA" ~ "Balanced accuracy",
			Metric == "TPR" ~ "True positive rate",
			Metric == "TNR" ~ "True negative rate"
		)
	)

taxonomy_only_CV_long <- taxonomy_CV_long %>%
	filter(Version=="FAPROTAX2-mapper\n(taxonomy-based)")
	

full_boxplot_CV <- rbind(sequence_CV_long, taxonomy_only_CV_long)

# Create df for scatterplots
sequence_scatterplot_df <- sequence_CV_long %>%
	pivot_wider(names_from="Version",
				values_from="Score")
				
taxonomy_scatterplot_df <- taxonomy_CV_long %>%
	pivot_wider(names_from="Version",
				values_from="Score")
				
summary_accuracy_table <- full_boxplot_CV %>%
	group_by(Version, Metric) %>%
	summarise(mean_score = mean(Score, na.rm=T))
	
full_accuracy_table <- full_boxplot_CV %>%
	group_by(Version, Metric, trait) %>%
	summarise(mean_score = mean(Score, na.rm=T))

summary_accuracy_table$Version <- gsub("\n", " ", summary_accuracy_table$Version)
full_accuracy_table$Version <- gsub("\n", " ", full_accuracy_table$Version)


full_boxplot_CV$Version <- factor(full_boxplot_CV$Version, levels=c("FAPROTAX v1", "FAPROTAX2-mapper\n(taxonomy-based)", "FAPROTAX2-mapper\n(sequence-based)"))

## Create boxplots
BA_boxplot <- ggplot(full_boxplot_CV %>% filter(Metric == "BA"), aes(x = Version, y = Score)) +
geom_boxplot(aes(fill=Version), outlier.shape = NA, width=0.4) + labs(x = "", y = "Balanced accuracy") + 
  geom_jitter(alpha = 0.5, width=0.03, size=1.2) + ylim(0,1) +
  theme_classic() + scale_fill_manual(values=c("FAPROTAX v1" = "grey50", "FAPROTAX2-mapper\n(taxonomy-based)" = "#ed8824", 
  												"FAPROTAX2-mapper\n(sequence-based)" = "#2880a8")) +
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill=NA, color="black", linewidth=1),
        legend.position = "none",
        axis.text.x=element_text(angle=45, hjust=1, size=8, color="black"))

TPR_boxplot <- ggplot(full_boxplot_CV %>% filter(Metric == "TPR"), aes(x = Version, y = Score)) +
  geom_boxplot(aes(fill=Version), outlier.shape = NA, width=0.4) + labs(x = "", y = "True positive rate") + 
  geom_jitter(alpha = 0.5, width=0.03, size=1.2) + ylim(0,1) +
  theme_classic() + scale_fill_manual(values=c("FAPROTAX v1" = "grey50", "FAPROTAX2-mapper\n(taxonomy-based)" = "#ed8824", 
  												"FAPROTAX2-mapper\n(sequence-based)" = "#2880a8")) +
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill=NA, color="black", linewidth=1),
        legend.position = "none",
        axis.text.x=element_text(angle=45, hjust=1, size=8, color="black"))

TNR_boxplot <- ggplot(full_boxplot_CV %>% filter(Metric == "TNR"), aes(x = Version, y = Score)) +
  geom_boxplot(aes(fill=Version), outlier.shape = NA, width=0.4) + labs(x = "", y = "True negative rate", fill="version") + 
  geom_jitter(alpha = 0.5, width=0.03, size=1.2) + ylim(0,1) +
  theme_classic() + scale_fill_manual(values=c("FAPROTAX v1" = "grey50", "FAPROTAX2-mapper\n(taxonomy-based)" = "#ed8824", 
  												"FAPROTAX2-mapper\n(sequence-based)" = "#2880a8")) +
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill=NA, color="black", linewidth=1),
        axis.text.x=element_text(angle=45, hjust=1, size=8, color="black"),
        legend.box=element_blank(),
        legend.position="none")

full_accuracy_boxplots <- BA_boxplot + TPR_boxplot + TNR_boxplot



## Create taxonomy-based scatterplots
taxonomy_BA_scatterplot <- ggplot(taxonomy_scatterplot_df %>% filter(Metric=="BA"), aes(x=`FAPROTAX v1`, y=`FAPROTAX2-mapper\n(taxonomy-based)`)) +
	geom_point(size=2, color="#ed8824", alpha=0.8) + geom_abline(color="grey70") + xlim(0,1) + ylim(0,1) +
	theme_classic() + labs(title="Balanced accuracy") + 
	theme(
		axis.line=element_blank(),
		panel.border=element_rect(fill=NA, color="black", linewidth=1),
		axis.title.y=element_text(size=9),
		title=element_text(size=8)
	)

taxonomy_TPR_scatterplot <- ggplot(taxonomy_scatterplot_df %>% filter(Metric=="TPR"), aes(x=`FAPROTAX v1`, y=`FAPROTAX2-mapper\n(taxonomy-based)`)) +
	geom_point(size=2, color="#ed8824", alpha=0.8) + geom_abline(color="grey70") + xlim(0,1) + ylim(0,1) +
	theme_classic() + labs(title="True positive rate") + 
	theme(
		axis.line=element_blank(),
		panel.border=element_rect(fill=NA, color="black", linewidth=1),
		axis.title.y=element_blank(),
		title=element_text(size=8)
	)
	
taxonomy_TNR_scatterplot <- ggplot(taxonomy_scatterplot_df %>% filter(Metric=="TNR"), aes(x=`FAPROTAX v1`, y=`FAPROTAX2-mapper\n(taxonomy-based)`)) +
	geom_point(size=2, color="#ed8824", alpha=0.8) + geom_abline(color="grey70") + xlim(0,1) + ylim(0,1) +
	theme_classic() + labs(title="True negative rate") + 
	theme(
		axis.line=element_blank(),
		panel.border=element_rect(fill=NA, color="black", linewidth=1),
		axis.title.y=element_blank(),
		title=element_text(size=8)
	)

full_accuracy_scatterplots_taxonomy <- taxonomy_BA_scatterplot + taxonomy_TPR_scatterplot + taxonomy_TNR_scatterplot

## Create sequence-based scatterplots
sequence_BA_scatterplot <- ggplot(sequence_scatterplot_df %>% filter(Metric=="BA"), aes(x=`FAPROTAX v1`, y=`FAPROTAX2-mapper\n(sequence-based)`)) +
	geom_point(size=2, color="#2880a8", alpha=0.8) + geom_abline(color="grey70") + xlim(0,1) + ylim(0,1) +
	theme_classic() + labs(title="Balanced accuracy") + 
	theme(
		axis.line=element_blank(),
		panel.border=element_rect(fill=NA, color="black", linewidth=1),
		axis.title.y=element_text(size=9),
		title=element_text(size=8)
	)

sequence_TPR_scatterplot <- ggplot(sequence_scatterplot_df %>% filter(Metric=="TPR"), aes(x=`FAPROTAX v1`, y=`FAPROTAX2-mapper\n(sequence-based)`)) +
	geom_point(size=2, color="#2880a8", alpha=0.8) + geom_abline(color="grey70") + xlim(0,1) + ylim(0,1) +
	theme_classic() + labs(title="True positive rate") + 
	theme(
		axis.line=element_blank(),
		panel.border=element_rect(fill=NA, color="black", linewidth=1),
		axis.title.y=element_blank(),
		title=element_text(size=8)
	)
	
sequence_TNR_scatterplot <- ggplot(sequence_scatterplot_df %>% filter(Metric=="TNR"), aes(x=`FAPROTAX v1`, y=`FAPROTAX2-mapper\n(sequence-based)`)) +
	geom_point(size=2, color="#2880a8", alpha=0.8) + geom_abline(color="grey70") + xlim(0,1) + ylim(0,1) +
	theme_classic() + labs(title="True negative rate") + 
	theme(
		axis.line=element_blank(),
		panel.border=element_rect(fill=NA, color="black", linewidth=1),
		axis.title.y=element_blank(),
		title=element_text(size=8)
	)

full_accuracy_scatterplots_sequence <- sequence_BA_scatterplot + sequence_TPR_scatterplot + sequence_TNR_scatterplot



full_figure <- full_accuracy_boxplots / full_accuracy_scatterplots_taxonomy / full_accuracy_scatterplots_sequence + plot_layout(heights = c(1.4, 1, 1))
ggsave(opt$out_figure, full_figure, width=8, height=9)

# Save output tables
data.table::fwrite(full_accuracy_table, opt$out_full_table, sep="\t")
data.table::fwrite(summary_accuracy_table, opt$out_summary_table, sep="\t")