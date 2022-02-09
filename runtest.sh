mirna_gff="test/data/sample.gff3"
vcf="test/data/sample.vcf.gz"
output_directory="$(dirname $mirna_gff)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run nf-miRNA-SNPs-classify.nf \
	--mirna_gff $mirna_gff \
	--vcf $vcf \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
