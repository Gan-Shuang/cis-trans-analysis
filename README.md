# cis-trans-analysis
## Usage
Analyse EGFR cis-trans type from mutation reads
## Method
 Calculate persentage of each EGFR mutation in total deepth

 python run_cis-trans.py -h
 optional arguments:
   -h, --help            show this help message and exit
   -v VCF_FILE, --vcf VCF_FILE
                         input vcf
   -b BAM_FILE, --bam BAM_FILE
                         input bam
   -o RESULT_FILE, --output RESULT_FILE
                         cis-trans analysis result
                        
## Test
sh ./test.sh
-sampleid	geneA	phgvsA	geneB	phgvsB	cis_trans	chgvsA	chgvsB	ratio
-test_cis	EGFR	p.T790M	EGFR	p.C797G	顺式	c.2369C>T	c.2389T>G	0.9552238805970149
