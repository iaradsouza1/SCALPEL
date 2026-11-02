#!/usr/bin/env python3
"""
Build the APA count matrix for one library, sparse, replacing SCALPEL's R path.

Why this exists: src/merge_DGE_PRED.R pivot_wider()s a dense isoform x cell table, and
tidyr computes `nrow * ncol` as a 32-bit int. Above 2,147,483,647 it overflows to NA and
the task dies:

    Error in `vec_rep_each()`: `times` can't be missing.
      vctrs::vec_rep_each(fill, nrow * ncol)
    Warning: In nrow * ncol : NAs produced by integer overflow

42 of 119 ROSMAP libraries (35%) are above that ceiling -- any library past ~14,600
cells at ~147k isoforms. It is a hard limit, so no amount of memory helps, and the
largest library that did pass cleared it by 0.07%. The matrices are also only ~7.3%
non-zero, so dense was never the right representation: ~31 GB of signal was being
written as ~430 GB of mostly zeros.

Arithmetic is identical to merge_DGE_PRED.R:
    tr_counts = tr_prob * counts        joined on (bc, gene_name)
    counts missing -> 0                 (so an absent join contributes nothing)
    row key      = "<gene_name>***<transcript_name>"

Output is an .h5ad with X as a CSR matrix of cells x isoforms (AnnData's orientation,
transposed relative to the old APADGE which was isoforms x cells).

float64 throughout, deliberately: tr_prob carries 3 decimals and counts are integers, so
the product is exact in double but not in float32 -- at float32 a value of 176.346 came
out 1.5e-5 off R's, the single discrepancy in an otherwise exact 126k-row comparison.
The extra bytes are irrelevant next to the 2.91 GB dense file this replaces.

NaN in tr_prob is still propagated rather than zeroed, matching R, but as of the v48 run
it should no longer occur. Those NaN were an em_algorithm.R bug, not a real quantity: a
gene whose fragments all carry probability 0 is discarded there, which returns an empty
table, and the left_join then left rel_abund NA for that gene's transcripts rather than
dropping them. fwrite wrote the NA as an empty string, fread read it back as NA, and
tr_prob * counts made it NaN -- 46 values in the library this script was validated against
and 35 across the B32 pilot. em_algorithm.R now filters on is.finite(rel_abund), so the
row is absent from the preds file entirely and contributes no COO entry, i.e. an implicit
zero, which is what "discarded" means. The guard here is kept as a tripwire: NaN reaching
this point again means something new, and `nz = tr_counts != 0` below keeps it (NaN != 0),
so it stays visible rather than silently becoming a zero.

Usage: dge_generation.py <isoforms_quantified.txt> <counts.txt> <sample_id> <out.h5ad>
"""

import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_preds(path):
    """EM output: bc, gene_name, transcript_name, rel_abund -- no header."""
    return pd.read_csv(
        path, sep="\t", header=None,
        names=["bc", "gene_name", "transcript_name", "tr_prob"],
        dtype={"bc": "category", "gene_name": "category",
               "transcript_name": "category", "tr_prob": "float64"},
    )


def read_counts(path):
    """read_10X.R output, long: GENE, bc, counts -- with header."""
    df = pd.read_csv(
        path, sep="\t",
        dtype={"GENE": "category", "bc": "category", "counts": "float64"},
    )
    df.columns = ["gene_name", "bc", "counts"]
    return df


def main():
    pred_path, dge_path, sample_id, out_path = sys.argv[1:5]

    preds = read_preds(pred_path)
    print(f"preds  : {len(preds):,} rows", flush=True)
    dge = read_counts(dge_path)
    print(f"dge    : {len(dge):,} rows", flush=True)

    # One shared code space for barcodes and for genes, so the two tables can be
    # joined on a single int64 key instead of a pandas merge over 36M rows.
    bc_cats = pd.Index(sorted(set(preds["bc"].cat.categories) | set(dge["bc"].cat.categories)))
    gn_cats = pd.Index(sorted(set(preds["gene_name"].cat.categories) | set(dge["gene_name"].cat.categories)))
    n_bc = len(bc_cats)

    p_bc = bc_cats.get_indexer(preds["bc"].astype(object)).astype(np.int64)
    p_gn = gn_cats.get_indexer(preds["gene_name"].astype(object)).astype(np.int64)
    d_bc = bc_cats.get_indexer(dge["bc"].astype(object)).astype(np.int64)
    d_gn = gn_cats.get_indexer(dge["gene_name"].astype(object)).astype(np.int64)

    # counts lookup by (gene, bc) via a sorted key array + searchsorted
    d_key = d_gn * n_bc + d_bc
    order = np.argsort(d_key, kind="stable")
    d_key = d_key[order]
    d_val = dge["counts"].to_numpy()[order]
    del dge, d_bc, d_gn, order

    p_key = p_gn * n_bc + p_bc
    pos = np.searchsorted(d_key, p_key)
    pos_clipped = np.clip(pos, 0, len(d_key) - 1)
    hit = d_key[pos_clipped] == p_key
    counts = np.where(hit, d_val[pos_clipped], 0.0).astype(np.float64)
    print(f"join   : {int(hit.sum()):,} / {len(p_key):,} pred rows matched a gene count", flush=True)
    del d_key, d_val, p_key, pos, pos_clipped, hit, p_gn

    tr_counts = preds["tr_prob"].to_numpy() * counts
    del counts

    # Row key "<gene>***<transcript>", matching the old APADGE's gene_transcript
    iso = (preds["gene_name"].astype(str) + "***" + preds["transcript_name"].astype(str))
    iso_codes, iso_cats = pd.factorize(iso, sort=True)
    del iso, preds

    nz = tr_counts != 0
    print(f"nonzero: {int(nz.sum()):,} of {len(tr_counts):,} "
          f"({100*nz.sum()/len(tr_counts):.2f}%)", flush=True)

    # cells x isoforms, so the dense product that broke R is never materialised
    X = sp.coo_matrix(
        (tr_counts[nz], (p_bc[nz], iso_codes[nz])),
        shape=(n_bc, len(iso_cats)), dtype=np.float64,
    ).tocsr()
    X.sum_duplicates()
    print(f"matrix : {X.shape[0]:,} cells x {X.shape[1]:,} isoforms, "
          f"{X.nnz:,} stored ({100*X.nnz/(X.shape[0]*X.shape[1]):.2f}% dense)", flush=True)

    import anndata as ad
    gene = pd.Index(iso_cats).str.split(r"\*\*\*", n=1).str[0]
    tx = pd.Index(iso_cats).str.split(r"\*\*\*", n=1).str[1]
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(bc_cats, name="cell_barcode")),
        var=pd.DataFrame({"gene_name": gene, "transcript_name": tx},
                         index=pd.Index(iso_cats, name="gene_transcript")),
    )
    adata.uns["sample_id"] = sample_id
    adata.uns["source"] = "SCALPEL dge_generation (python, sparse)"
    adata.write_h5ad(out_path, compression="gzip")
    print(f"wrote {out_path} ({os.path.getsize(out_path)/1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
