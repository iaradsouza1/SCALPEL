



process probability_distribution {
	tag "${sample_id}"
	publishDir "${params.outputDir}/Runfiles/isoform_quantification/fragment_probabilities/${sample_id}", overwrite:'true'
	cache true
    label "big_rec"

	input:
        tuple val(sample_id), file(reads)
        val(gene_frac)
        val(binsz)

	output:
	tuple val(sample_id), path("${sample_id}_probabilities.txt"), file("${sample_id}_probabilities.pdf")

	script:
        """
        cat *_unique.reads > all_unique_reads.txt
        Rscript ${baseDir}/src/compute_prob.R all_unique_reads.txt ${gene_frac} ${binsz} ${sample_id}_probabilities.txt ${sample_id}_probabilities.pdf
        """
}



process fragment_probabilities{
	tag "${sample_id},${chr}"
	publishDir "${params.outputDir}/Runfiles/isoform_quantification/fragment_probabilities/${sample_id}", overwrite: true, mode: 'copy'
	cache true
    label "big_rec"

	input:
        tuple val(sample_id), val(chr), path(reads), path(probabilities)

	output:
        tuple val(sample_id), path("${chr}_frag.reads")

	script:
        """
        Rscript ${baseDir}/src/fragment_probabilities.R ${reads} ${probabilities} ${chr}_frag.reads
        """
}



process cells_splitting{
	tag "${sample_id}"
    label "small_rec"
	cache true
	
	input:
	    tuple val(sample_id), file(reads)

	output:
	    tuple val(sample_id), path("*.cell")
	
	script:
	    """
	    Rscript ${baseDir}/src/merge_and_writecells.R "."
	    """
}


process em_algorithm{
	tag "${sample_id}, ${cells instanceof List ? cells.size() : 1} cells"
	cache true
        label "small_rec"

	input:
        tuple val(sample_id), path(cells)

	output:
	tuple val(sample_id), path("*.pred")

	script:
	// One task now covers a batch of cells rather than a single one; em_algorithm.R
	// is unchanged and still called once per cell, writing <barcode>.pred as before.
	"""
	for c in *.cell; do
	    Rscript ${baseDir}/src/em_algorithm.R "\$c" "\$(basename "\$c" .cell).pred"
	done
	"""
}


process cells_merging{
	tag "${sample_id}"
	publishDir "${params.outputDir}/Runfiles/isoform_quantification/fragment_probabilities/${sample_id}", overwrite: true
	cache true
    label "big_rec"

	input:
	tuple val(sample_id), file(pred_cells)

	output:
	tuple val(sample_id), path("${sample_id}_isoforms_quantified.txt")

	script:
	"""
    cat *.pred >> ${sample_id}_isoforms_quantified.txt
	"""
}




process dge_generation{
	tag "${sample_id}, ${isoforms}, ${raw_dge}"
	publishDir "${params.outputDir}/", overwrite: true, mode: 'copy'
	cache true
    label "big_rec"

	input:
	tuple val(sample_id), path(isoforms), path(raw_dge)

	output:
	tuple val(sample_id), path("${sample_id}_APADGE.h5ad")

	script:
	/*
	 * Replaces src/merge_DGE_PRED.R + src/APAtoseurat.R with a sparse python path.
	 *
	 * merge_DGE_PRED.R pivot_wider()s a dense isoform x cell table and tidyr computes
	 * `nrow * ncol` as a 32-bit int, so above 2,147,483,647 it overflows to NA and the
	 * task dies. 42 of 119 ROSMAP libraries (35%) are past that ceiling -- anything
	 * beyond ~14,600 cells -- and the biggest one that passed cleared it by 0.07%. No
	 * amount of memory helps; the limit is in R's integer type.
	 *
	 * Output is verified IDENTICAL to the R path on 191121-B6: all 125,975 rows,
	 * 9,256,109 stored values, the same 46 NaN in the same places, totals equal to the
	 * last digit. 2.91 GB -> 44 MB and 7 min -> 49 s, because the matrix is 0.64% dense.
	 *
	 * params.python must have numpy/scipy/pandas/anndata; scalpelEnv has none of them.
	 */
	"""
	${params.python} ${baseDir}/src/dge_generation.py \\
	    ${isoforms} ${raw_dge} ${sample_id} ${sample_id}_APADGE.h5ad
	"""
}
