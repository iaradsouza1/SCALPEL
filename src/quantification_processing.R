#!/usr/bin/env Rscript

# =============================================================================
# Script: quantification_processing.R
# Description: Merges Salmon quantifications and GTF for transcript abundance
# Inputs:
#   1. Genome GTF file path
#   2. Output file path for merged TPM table
# Outputs:
#   - Merged quantification table with TPM percentages (to args$output_path)
#   - Split GTFs per chromosome in the current directory
# =============================================================================


#libs
suppressWarnings({
  library(argparse, quietly=T, verbose=F)
  library(dplyr, quietly=T, verbose=F)
})

#- Argument parsing
parser <- ArgumentParser(description = 'Process Salmon quantifications and GTF')
parser$add_argument('genome', type = "character", help = 'GTF file path')
parser$add_argument('output_path', type = "character", help = 'Output path for merged TPM summary')
parser$add_argument('pseudocount', type = "double", nargs = '?', default = 1,
                    help = 'added to each TPM before averaging; guards the within-gene division and regularises lowly-expressed genes')
args <- parser$parse_args()

#- Opening
#gtf.file...
gtf = rtracklayer::import(args$genome)
gtf = data.frame(gtf[gtf$type=="exon",])
gtf.attributes = c("seqnames", "start", "end", "width", "strand", "gene_id", "gene_name", "transcript_id", "transcript_name")
gtf = gtf[, gtf.attributes]

#bulk quantification...
bulk.files = list.files(path = ".", pattern = "*.sf", full.names = T)

#- Processing bulk quantification sample files... and write merge_quant file
RECEPT = lapply(bulk.files, function(x){
    #reading bulk counts
    curr = data.table::fread(x, nThread=1)
    curr$sample = basename(x)
    return(curr)
}) %>% data.table::rbindlist() %>%
    # Pseudocount: guards the within-gene division below, and pulls lowly-expressed
    # genes' priors toward uniform. At the original value of 1 that pull was strong --
    # on the B32 pilot, 28.9% of multi-isoform genes ended up within 15% of uniform,
    # i.e. the prior told the EM nothing. At 0.1 that falls to 18.6% while only 395
    # genes (1.8%) gain an isoform whose prior drops below 1% of uniform, which would
    # effectively exclude it from the EM. Going to 0.01 pushes that to 18.6% of genes,
    # which asks too much of salmon's zeros: ~20% of reads map here, 3'-tagged against
    # full-length transcript models, so a zero TPM is weak evidence of absence.
    dplyr::mutate(TPM = TPM + args$pseudocount) %>%
    #get MeanTPM & discard null transcrips in bulk
    dplyr::group_by(Name) %>%
    dplyr::reframe(meanTPM = mean(TPM)) %>%
    # dplyr::filter(meanTPM!=0) %>%
    dplyr::rename(transcript_id="Name") %>%
    #include gene metadata
    dplyr::left_join(distinct(gtf, gene_name, transcript_name, transcript_id), by="transcript_id") %>%
    #only considers all the transcript in the annotation GTF...
    na.omit() %>%
    #calculate Transcripts bulk % in gene
    dplyr::group_by(gene_name) %>%
    dplyr::reframe(transcript_name, bulk_TPMperc = meanTPM / sum(meanTPM)) %>%
    dplyr::arrange(gene_name,transcript_name) %>%
    data.table::data.table() %>%
    #write bulk quantification file
    data.table::fwrite(file=args$output_path, nThread=1, row.names=F, sep="\t")

#- Write GTF by chromosomes...
data.table::fwrite(gtf, file = "reduced_annotation.tsv", sep="\t", nThread=1, row.names=F, col.names=F)
