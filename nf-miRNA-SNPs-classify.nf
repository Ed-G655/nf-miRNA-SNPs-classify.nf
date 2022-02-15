#!/usr/bin/env nextflow

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
 Israel Aguilar-Ordonez (iaguilaror@gmail)


=============================
Pipeline Processes In Brief:

Pre-processing:
PRE1-CONVERT-GFF-TO-BED
PRE2_EXTRACT_MIRNA_SEED
PRE3_INTERSEPT_BED_VCF

Core-processing
CLASSIFY_SNPs_BY_REGION

Pos-processing
POS_1_PLOT_ASSOCIATED_DISEASES

Analysis:

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
	The miRNA SNPs classify
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name} --mirna_gff <path to input 1> --vcf <path to input 2>

	  --mirna_gff	<- miRNA gff file;

    --vcf <- VCF file;

	  --output_dir     <- directory where results, intermediate and log files will be stored;
	      default: same dir where --query_fasta resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-miRNA-SNPs-classify.nf"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.mirna_gff = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.vcf = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.mirna_gff | !params.vcf ) {
  log.error " Please provide the --mirna_gff AND --vcf \n\n" +
  " For more information, execute: nextflow run  ${pipeline_name} --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.mirna_gff).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The ${pipeline_name} pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['Input miRNA GFF']= params.mirna_gff
pipelinesummary['Input VCF']			= params.vcf
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/* enable DSL2*/
nextflow.enable.dsl=2

/*
	READ GENERAL INPUTS
*/

/* Load VCF file into channel */
Channel
	.fromPath( "${params.vcf}" )
	.set{ vcf_input }

	/* Load GFF file into channel */
Channel
	.fromPath( "${params.mirna_gff}" )
	.set{ gff_input }


/*
 Load R fileS
*/

/* R_script_1 */
Channel
	.fromPath( "./nf-modules/Rscripts/extract_mirna_regions.R" )
	.set{ R_script_1 }

/* R_script_1 */
Channel
	.fromPath( "./nf-modules/Rscripts/miRNA_SNPs.R" )
	.set{ R_script_2 }

/* R_script_1 */
Channel
	.fromPath( "./nf-modules/Rscripts/miRNA_diseases.R" )
	.set{ R_script_3 }

/*	 * Import modules */
include { CONVERT_GFF_TO_BED;
					EXTRACT_MIRNA_SEED;
					INTERSEPT_BED_VCF;
					INTERSEPT_BED_VCF as INTERSEPT_SEED_BED_VCF;
					CLASSIFY_SNPs_BY_REGION;
					PLOT_ASSOCIATED_DISEASES
					} from './nf-modules/modules.nf'

/*  main pipeline logic */
workflow  {
						// PRE 1: CONVERT_GFF_TO_BED
						CONVERT_GFF_TO_BED(gff_input)
						EXTRACT_MIRNA_SEED(CONVERT_GFF_TO_BED.out, R_script_1);
						MATURE_INTERSEPT = INTERSEPT_BED_VCF(CONVERT_GFF_TO_BED.out, vcf_input)
						SEED_INTERSEPT = INTERSEPT_SEED_BED_VCF(EXTRACT_MIRNA_SEED.out, vcf_input)
						CLASSIFY_SNPs_BY_REGION(MATURE_INTERSEPT, SEED_INTERSEPT, R_script_2)
						PLOT_ASSOCIATED_DISEASES(CLASSIFY_SNPs_BY_REGION.out.mirna_snps, R_script_3)

}
