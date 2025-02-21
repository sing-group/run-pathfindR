# Run pathfindR

The `run-pathfindR.R` script works in [pathfindR](https://github.com/egeulgen/pathfindR) v2.4.2 and implements the most common pathfindR workflow using these four input files:
1. DEA results: path to the input file with the DEA results. It must be a CSV file with at least these four columns: gtf_ensembl_id, gene_name, logFC, and padj.
2. Counts file: path to the CSV file with the counts for the samples (rows are genes and columns are samples; genes are specified as Ensembl IDs).
3. Metadata file: path to the TSV file with the samples metadata.
4. Reference file: path to the file with the reference condition (only one line with the condition name).

The script also takes three additional parameters:

5. output directory: path to the directory where results should be stored.
6. gene sets: the gene sets to be used for enrichment analysis, one of: KEGG, Reactome, BioCarta, GO-All, GO-BP, GO-CC and GO-MF (all for Homo sapiens)
7. pin: the protein interaction network to be used for enrichment analysis, one of: Biogrid, STRING, GeneMania, IntAct, KEGG, and mmu_STRING

## Motivation

This script was created in the context of a [Compi pipeline for RNA-Seq data analysis](https://github.com/sing-group/compi-rnaseq-pipeline). In this pipeline, DElite was used for performing DEA and so this is the input CSV file to the script. As the original GTF annotation uses EnsemblIDs, this file was postprocessed to add gene names. These two fields are used to convert the counts file gene EnsemblIDs into gene names before using it in pathfindR. Hopefully this script can be reused for similar purposes!

## Test

The example files are provided in the `test` folder. They can be used to run the script using the `pegi3s/r_pathfindr:2.4.2` Docker image from the [Bioinformatics Docker Images Project](http://bdip.i3s.up.pt/container/r_pathfindr):
```shell
docker run --rm -v $(pwd):$(pwd) -w $(pwd) \
    pegi3s/r_pathfindr:2.4.2 \
        Rscript run-pathfindR.R \
            test/test.csv \
            test/counts.tsv \
            test/metadata.tsv \
            test/reference.txt \
            test/results \
            KEGG \
            STRING
```
