# Human Pseudoautosomal Boundary Analysis
**Daniel Cotter**

Arizona State University

Wilson Sayres Lab

2018

## Analysis Steps

### Step 1: Install/Set-up

Clone this repository and make sure [conda](https://conda.io/docs/user-guide/install/index.html) is installed. Then use the provided *PAB_variation.yml* environment file to set up a new conda environment.

```shell
conda env create -f PAB_variation.yml
```

### Step 2: Get data

Data should be downloaded from the provided links and copied into the **data/** directory of the project folder.

#### Variant Data
Variant data is from The 1000 Genomes Project phase3 VCF files for chrX, chrY, & chr8

- **chrX:** [ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz "chrX")

- **chrY:** [ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz "chrY")

- **chr8:** [ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz "chr8")

#### Genome masks
The strict mask as provided by The 1000 genomes Project is used for identifying monomorphic sites

- **Strict Mask chrX:** [20141020.chrX.strict_mask.fasta.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/20141020.chrX.strict_mask.fasta.gz "chrX")

- **Strict Mask chrY:** [20141020.chrY.strict_mask.fasta.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/20141020.chrY.strict_mask.fasta.gz "chrY")

- **Strict Mask chr8:** [20140520.chr8.strict_mask.fasta.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/20140520.chr8.strict_mask.fasta.gz "chr8")

#### Population Lists
Population and subpopulation lists are calculated using [integrated_call_samples_v3.20130502.ALL.panel](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel "Population panel")


We analyze all individuals across the X chromosome and chromosome 8 and we analyze all males across the Y chromosome. We have included the breakdown of the number of males and females in each population.

##### Super Populations

POP | Females | Males
---|:---:|:---:
AFR|342|319
AMR|177|170
EAS|260|244
EUR|263|240
SAS|229|260

##### Populations

POP | Females | Males
---|:---:|:---:
ACB|49|47
ASW|35|26
BEB|44|42
CDX|49|44
CEU|50|49
CHB|57|46
CHS|53|52
CLM|51|43
ESN|46|53
FIN|61|38
GBR|45|46
GIH|47|56
GWD|58|55
IBS|53|54
ITU|43|59
JPT|48|56
KHV|53|46
LWK|55|44
MSL|43|42
MXL|32|32
PEL|44|41
PJL|48|48
PUR|50|54
STU|47|55
TSI|54|53
YRI|56|52

Population codes can be found [here](http://www.internationalgenome.org/faq/which-populations-are-part-your-study/).

### Step 3: Run the analysis

Run the analysis using `snakemake` once all of the raw data files are downloaded. Navigate to the top of the project directory and type the following commands:

```shell
source activate PAB_variation
snakemake
```

## What do the rules in the Snakefile do?
### Rule 1: Parse populations

We use a simple python script to parse the panel file into several lists with the IDs of each individual broken up by **male/female/individual** and by **super-population/population**. The resulting files are simple text-files stored in the results subdirectory of the `01_populations/` folder

### Rule 2: Calculate diversity

Diversity is calculated as **π**, or the average number of pairwise nucleotide differences per site.  is calculated for our purposes as ***π = ( ∑ ( n<sub>i</sub> choose 2 ) ) / ( n choose 2 )***. Diversity is calculated at each site for each population that is supplied in a separate population file. We calculate diversity for each chromosome (chrX, chrY, chr8) for all populations at once.

### Rule 3: Create filter

The filter is created by merging a collection of bed files (for each chromosome) that are taken from the UCSC genome browser. Some filters and how we obtain them are as follows:
- **Coding sequence:** We filter for coding sequence by merging known genes from UCSC, RefSeq, and GENCODE. We use the following commands to accomplish this:
 ```shell
 cat file1 file2 file3 | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' | bedtools merge -i stdin > merged_genes_filter.bed
 ```
- **Repetitive elements:** We get repetitive elements from the simple-repeats track.
- **Segmental duplications:** We get segmental duplications from the Segmental Dups track.
- **Telomeres & Centromere:** We use the Gap table and define the gap types that we want in the filter (*"centromere telomere"*).
