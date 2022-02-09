#!/usr/bin/env nextflow
/*================================================================

								---- MODULE PIPELINE ---------

/*================================================================
The Aguilar Lab presents...

- A pipeline to classify SNPs in microRNA regions and provide an overview of
diseases associated with microRNAs that present SNPs

==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)



- Bioinformatics Development
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)


- Nextflow Port
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)

///////////////////////////////////////////////////////////////

  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/

pipeline_name = "nf-miRNA-SNPs-classify.nf"

/*This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*================================================================/*

/* MODULE START */

/* PRE1_CONVERT_GFF_TO_BED */

process CONVERT_GFF_TO_BED {
	tag "$GFF"

	publishDir "${results_dir}/PRE1-CONVERT-GFF-TO-BED /",mode:"copy"

	input:
	file GFF

	output:
	file "*.bed"

	shell:
	"""
	gff2bed < !{GFF} > !{GFF.baseName}.tmp
	less -S !{GFF.baseName}.tmp \
	| tr ";" "\t" \
	| sed -r 's/Name=//g' \
	| awk '{print \$1"\t"\$2"\t"\$3"\t"\$12"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' > !{GFF.baseName}.bed
	rm *.tmp
	"""
	stub:
	"""
	      touch ${GFF.baseName}.bed
	"""
}


process EXTRACT_MIRNA_SEED {
	tag "$BED"

	publishDir "${results_dir}/PRE2-EXTRACT-MIRNA-SEED/",mode:"copy"

	input:
	file BED
	file R_script_1

	output:
	file "*.seed.bed"

	shell:
	"""
	Rscript --vanilla ${R_script_1} ${BED} ${BED.baseName}.seed.bed
	"""
	stub:
	"""
	      touch ${BED.baseName}.seed.bed
	"""
}


process INTERSEPT_BED_VCF {
	tag "$BED, $VCF"

	publishDir "${results_dir}/PRE3_INTERSEPT_BED_VCF/",mode:"copy"

	input:
	file BED
	file VCF

	output:
	file "*.tsv"

	shell:
	"""
	bedtools intersect -a ${VCF} -b ${BED} -wb > ${BED.baseName}.tsv
	"""
	stub:
	"""
	      touch ${BED.baseName}.tsv
	"""
}


process CLASSIFY_SNPs_BY_REGION {
	tag "$MIRNAS, $MIRNAS_SEED"

	publishDir "${results_dir}/CLASSIFY_SNPs_BY_REGION/",mode:"copy"

	input:
	file MIRNAS
	file MIRNAS_SEED
	file R_script_2

	output:
	path "*.tsv", emit: mirna_snps
	file "*.png"

	shell:
	"""
	Rscript --vanilla ${R_script_2} ${MIRNAS} ${MIRNAS_SEED} ${MIRNAS.baseName}_out
	"""
	stub:
	"""
	      touch ${MIRNAS.baseName}_out.tsv
				touch ${MIRNAS.baseName}.png
	"""
}


process PLOT_ASSOCIATED_DISEASES {
	tag "$mirna_snps"

	publishDir "${results_dir}/POS_1_PLOT_ASSOCIATED_DISEASES/",mode:"copy"

	input:
	file mirna_snps
	file R_script_3

	output:
	file "*.tsv"
	file "*.png"

	shell:
	"""
	echo "Dowload database from HMDD https://www.cuilab.cn/hmdd"
	curl -o HMDD_alldata.txt https://www.cuilab.cn/static/hmdd3/data/alldata.txt
	Rscript --vanilla ${R_script_3} ${mirna_snps} ${mirna_snps.baseName}
	"""
	stub:
	"""
	      touch ${mirna_snps.baseName}_out_diseases.tsv
				touch ${mirna_snps.baseName}_primary.png
				touch ${mirna_snps.baseName}_mature.png
				touch ${mirna_snps.baseName}_seed.png
	"""
}
