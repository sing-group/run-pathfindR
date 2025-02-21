## Script input parameters:
##  1.- input file: path to the input file with the DEA results. It must be a CSV file with 
##		at least these four columns: gtf_ensembl_id, gene_name, logFC, and padj
##  2.- counts file: path to the CSV file with the counts for the samples (rows are genes and
##		columns are samples; genes are specified as Ensembl IDs).
##  3.- metadata file: path to the TSV file with the samples metadata.
##  4.- reference file: path to the file with the reference condition (only one line with the condition name).
##  5.- output directory: path to the directory where results should be stored.
##  6.- gene sets: the gene sets to be used for enrichment analysis, one of:
##      KEGG, Reactome, BioCarta, GO-All, GO-BP, GO-CC and GO-MF (all for Homo sapiens)
##  7.- pin: the protein interaction network to be used for enrichment analysis, one of:
##		Biogrid, STRING, GeneMania, IntAct, KEGG, and mmu_STRING

suppressPackageStartupMessages(library(pathfindR))

args <- commandArgs(TRUE)

## Parse input parameters
input_file <- args[1];
if (!file(input_file)) {
	stop(paste("Input file not found. Provided path: ", input_file, sep=""))
}

input_counts_file <- args[2];
if (!file(input_counts_file)) {
	stop(paste("Input counts file not found. Provided path: ", input_counts_file, sep=""))
}

input_metadata_file <- args[3];
if (!file(input_metadata_file)) {
	stop(paste("Input metadata file not found. Provided path: ", input_metadata_file, sep=""))
}

input_reference_file <- args[4];
if (!file(input_reference_file)) {
	stop(paste("Input reference file not found. Provided path: ", input_reference_file, sep=""))
}

outputDirectory <- args[5];
if(substring(outputDirectory, nchar(outputDirectory)) != "/") {
	outputDirectory <- paste(outputDirectory, "/", sep="");
}
if (!dir.exists(outputDirectory)) {
	dir.create(outputDirectory, recursive = TRUE)
}

gene_sets = args[6]
VALID_GENE_SET_NAMES = c("KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC", "GO-MF")
if (!(gene_sets %in% VALID_GENE_SET_NAMES)) {
	stop(paste("Invalid gene set name. Provided gene set: ", gene_sets, ". It must be one of: ", paste(VALID_GENE_SET_NAMES, collapse=", "), sep=""))
}

pin = args[7]

VALID_PIN_NAMES = c("Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING")
if (!(pin %in% VALID_PIN_NAMES)) {
	stop(paste("Invalid PIN name. Provided PIN: ", pin, ". It must be one of: ", paste(VALID_PIN_NAMES, collapse=", "), sep=""))
}

tsv <- read.table(file = input_file, sep = ",", header = TRUE)
pathfindR_input <- tsv[, c('gene_name', 'logFC', 'padj')]

#
# 1. Use the run_pathfindR() wrapper function for the active-subnetwork-oriented enrichment analysis
#
RA_output <- run_pathfindR(
	pathfindR_input,
	output = paste(outputDirectory,"results", sep=""),
	gene_sets = gene_sets,
	plot_enrichment_chart = FALSE,
	pin_name_path = pin
)
write.table(RA_output, row.names = FALSE, paste(outputDirectory, "enriched_terms.tsv", sep=""), sep="\t")

if(nrow(RA_output) == 0){
	print("No enriched terms found, finishing analysis here.")
	q("no")
}

# 1.1 Plot the enrichment summary
pdf(paste(outputDirectory, 'enriched_terms_chart.pdf', sep=""))
enrichment_chart(RA_output, plot_by_cluster = FALSE, top_terms = 20)
dev.off()

#
# 2. Use the wrapper function cluster_enriched_terms() to perform clustering of enriched terms and 
# partitioning the terms into biologically-relevant groups. 
#
pdf(paste(outputDirectory, 'clustering.pdf', sep=""))
example_pathfindR_output_clustered <- cluster_enriched_terms(RA_output, plot_dend = TRUE, plot_clusters_graph = TRUE, plot_hmap = TRUE)
dev.off()

write.table(example_pathfindR_output_clustered, row.names = FALSE, paste(outputDirectory, "clustered_terms.tsv", sep=""), sep="\t")

representative <- example_pathfindR_output_clustered[example_pathfindR_output_clustered$Status == "Representative", ]
write.table(representative, row.names = FALSE, paste(outputDirectory, "clustered_terms_only_representative.tsv", sep=""), sep="\t")

# 2.1 Plot the enrichment summary
pdf(paste(outputDirectory, 'clustered_terms_chart.pdf', sep=""))
enrichment_chart(example_pathfindR_output_clustered, plot_by_cluster = TRUE, top_terms = 20)
dev.off()

pdf(paste(outputDirectory, 'enrichment_summary_clustered_representative.pdf', sep=""))
enrichment_chart(representative, plot_by_cluster = TRUE, top_terms = 20)
dev.off()

#
# 3. Aggregated Term Scores per Sample
#
counts <- read.table(input_counts_file, head= TRUE, sep = "\t")

# 
# 3.1 Prepare the input data
# The score_matrix function requires that exp_mat contains the experiment (e.g., gene expression) matrix.
# In this matrix, columns are samples and rows are genes. 
# Also, column names must contain sample names and row names must contain the gene symbols.
# Therefore, first we must assign gene symbols to the rows of the experiment matrix which are currently Ensembl IDs.
# 

mapping <- tsv[, c('gtf_ensembl_id', 'gene_name')]

# From counts we only keep those rows that have a gene_name in mapping
merged <- merge(counts, mapping, by.x = "gene", by.y = "gtf_ensembl_id", all.x = FALSE, all.y = TRUE)

# Remove all rows with gene = NA
merged <- merged[!is.na(merged$gene_name),]
rownames(merged) <- merged$gene_name

# And now remove gene and gene_name columns
merged <- merged[, !colnames(merged) %in% c("gene", "gene_name")]
merged <- as.matrix(merged)

# Load the metadata file to extract classes
metadata <- read.table(input_metadata_file, head= TRUE, sep = "\t")
reference_condition <- as.character(read.table(input_reference_file, head= FALSE, sep = "\t")[1,1])
case_condition <- as.character(setdiff(metadata$class, c(reference_condition))[1])
cases <- metadata[metadata$class == case_condition, 'sample']

# 3.2 Calculate scores for all terms and plot heat map using term descriptions
pdf(paste(outputDirectory, 'heatmap_aggregated_term_scores_per_sample.pdf', sep=""))
score_matrix <- score_terms(
  enrichment_table = RA_output,
  exp_mat = merged,
  cases = cases,
  use_description = TRUE,
  label_samples = FALSE,
  case_title = case_condition,
  control_title = reference_condition,
  low = "#f7797d",
  mid = "#fffde4",
  high = "#1f4037",
  plot_hmap = TRUE
)
dev.off()

# 3.2 Calculate scores for representative terms and plot heat map using term descriptions
representative_df <- example_pathfindR_output_clustered[example_pathfindR_output_clustered$Status == "Representative", ]
pdf(paste(outputDirectory, 'heatmap_aggregated_representative_term_scores_per_sample.pdf', sep=""))
score_matrix <- score_terms(
  enrichment_table = representative_df,
  exp_mat = merged,
  cases = cases,
  use_description = TRUE,
  label_samples = FALSE,
  case_title = case_condition,
  control_title = reference_condition,
  low = "#f7797d",
  mid = "#fffde4",
  high = "#1f4037",
  plot_hmap = TRUE
)
dev.off()

#
# 4. Enriched Term Diagrams
#

input_processed <- input_processing(pathfindR_input, pin_name_path = pin)

# 4.1 Generate interaction diagrams first for non-KEGG enrichment analyses
gg_list <- visualize_terms(
  result_df = RA_output,
  input_processed = input_processed,
  is_KEGG_result = FALSE,
  pin_name_path = pin
)

terms_directory = paste0(outputDirectory, "terms/")
if (!dir.exists(terms_directory)) {
	dir.create(terms_directory, recursive = TRUE)
}

for (term_id in names(gg_list)) {
	ggplot2::ggsave(
		filename = paste0(terms_directory, term_id, "_diagram.pdf"),
		plot = gg_list[[term_id]],
		width = 12,
		height = 12
	)
}

# For KEGG enrichment analyses, visualize_terms() can be used to generate KEGG pathway diagrams that are returned as a list of ggraph objects (using ggkegg):

if (gene_sets == "KEGG") {
	gg_list_kegg <- visualize_terms(
  		result_df = RA_output,
  		input_processed = input_processed,
  		is_KEGG_result = TRUE,
		pin_name_path = pin
	)
	kegg_terms_directory = paste0(outputDirectory, "terms_kegg/")
	if (!dir.exists(kegg_terms_directory)) {
		dir.create(kegg_terms_directory, recursive = TRUE)
	}

	for (term_id in names(gg_list_kegg)) {
		ggplot2::ggsave(
			filename = paste0(kegg_terms_directory, term_id, "_diagram.pdf"),
			plot = gg_list_kegg[[term_id]],
			width = 12,
			height = 12
		)
	}
}
