library(dplyr)
library(data.table)
library(argparse)
library(stringr)
library(tidyr)


#Argument Parser
#===============
parser = ArgumentParser(description='Merging DGE & isoforms Prediction')
parser$add_argument('DGE', type="character", help='path of DGE file')
parser$add_argument('PRED', type="character", help='path of predicted isoforms')
parser$add_argument('OUTP', type="character", help="path of output file")
args = parser$parse_args()

#Files opening
preds = fread(args$PRED, col.names=c("bc","gene_name","transcript_name","tr_prob"))
dge = fread(args$DGE)

#renaming
# read_10X.R now writes the DGE long (GENE / bc / counts, nonzero only), so the melt()
# that used to reshape a dense genes x cells table is no longer needed. The left_join
# below is unaffected: absent (bc, gene) pairs become NA and are coerced to 0, exactly
# as a stored zero would have been.
colnames(dge) = c("gene_name","bc","counts")

#left_joining by barcodes & gene_name
dge_preds = left_join(preds, dge, by=c("bc", "gene_name")) %>% 
    mutate(counts=ifelse(is.na(counts), 0, counts)) %>%
    mutate(tr_counts=tr_prob * counts, gene_transcript = paste0(gene_name, "***", transcript_name)) %>% 
    distinct(gene_transcript, bc, tr_counts) %>% 
    pivot_wider(names_from=bc, values_from=tr_counts, values_fill = 0)

print(head(dge_preds))

#writing
fwrite(dge_preds, file=args$OUTP, sep="\t")
