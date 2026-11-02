

library(argparse)
library(data.table)
library(dplyr)
library(Seurat)


#ArgumentParser
#==============
parser = ArgumentParser(description='Process 10X repo file')
parser$add_argument('repo_path', type="character", help='path of 10x repo')
args = parser$parse_args()


#define Sample_name
sample_name = basename(args$repo_path)

#Reading
print("reading file...")
s.mat = Read10X(paste0(args$repo_path,"/outs/filtered_feature_bc_matrix"))


#Writing
# Emit the DGE in LONG form, nonzero entries only, instead of densifying.
#
# data.frame(Read10X(...)) materialised a genes x cells dense table -- 738,608,180
# values for this cohort's largest library (36,601 x 20,180) against 67,544,153
# actually nonzero, 11x. merge_DGE_PRED.R then melt()ed that back to long, so the
# dense table only ever existed to be reshaped away.
#
# Output is unchanged: merge_DGE_PRED.R coerces an unmatched join to zero
# (counts = ifelse(is.na(counts), 0, counts)), so an absent row and a stored zero
# give the same result.
#
# Working from the sparse matrix also drops the str_replace() repair below, which
# only existed to undo the "-" -> "." mangling that data.frame() inflicted on the
# barcode column names.
nz = Matrix::summary(s.mat)
data.table(
  GENE   = rownames(s.mat)[nz$i],
  bc     = colnames(s.mat)[nz$j],
  counts = nz$x
) %>% fwrite(file = paste0(sample_name, ".counts.txt"), sep="\t", row.names=F)
