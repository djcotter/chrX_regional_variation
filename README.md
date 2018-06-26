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
The strict mask as provided by The 1000 genomes Project is used for identifying monomorphic sites. It is provided as a whole genome bed file:

- **Strict Mask:** [20141020.strict_mask.whole_genome.bed](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed)

#### Population Lists
Population and subpopulation lists are calculated using:
 - **Population Panel:**  [integrated_call_samples_v3.20130502.ALL.panel](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel "Population panel")


We analyze all individuals across the X chromosome and chromosome 8 and we analyze all males across the Y chromosome. The option to analyze males or females for chrX and chr8 can be changed at the top by altering the global variable, **SEX**. Below, we have included the breakdown of the number of males and females in each population:

<table>
<tr><th> AFR Populations </th><th> AMR Populations </th><th> EAS Populations </th><th> EUR Populations </th><th> SAS Populations </th></tr>
<tr><td>

<table>
<tr><th>POP</th><th>Females</th><th>Males</th></tr>
<tr><td>ACB</td><td>49</td><td>47</td></tr>
<tr><td>YRI</td><td>56</td><td>52</td></tr>
<tr><td>ASW</td><td>35</td><td>26</td></tr>
<tr><td>ESN</td><td>46</td><td>53</td></tr>
<tr><td>MSL</td><td>43</td><td>42</td></tr>
<tr><td>GWD</td><td>58</td><td>55</td></tr>
<tr><td>LWK</td><td>55</td><td>44</td></tr>
<tr><th>TOTAL</th><th>342</th><th>319</th></tr>
</table>
</td><td>
<table>
<tr><th>POP</th><th>Females</th><th>Males</th></tr>
<tr><td>MXL</td><td>32</td><td>32</td></tr>
<tr><td>PUR</td><td>50</td><td>54</td></tr>
<tr><td>CLM</td><td>51</td><td>43</td></tr>
<tr><td>PEL</td><td>44</td><td>41</td></tr>
<tr><th>TOTAL</th><th>177</th><th>170</th></tr>
</table>
</td><td>
<table>
<tr><th>POP</th><th>Females</th><th>Males</th></tr>
<tr><td>CHB</td><td>57</td><td>46</td></tr>
<tr><td>JPT</td><td>48</td><td>56</td></tr>
<tr><td>CHS</td><td>53</td><td>52</td></tr>
<tr><td>CDX</td><td>49</td><td>44</td></tr>
<tr><td>KHV</td><td>53</td><td>46</td></tr>
<tr><th>TOTAL</th><th>260</th><th>244</th></tr>
</table>
</td><td>
<table>
<tr><th>POP</th><th>Females</th><th>Males</th></tr>
<tr><td>CEU</td><td>50</td><td>49</td></tr>
<tr><td>TSI</td><td>54</td><td>53</td></tr>
<tr><td>FIN</td><td>61</td><td>38</td></tr>
<tr><td>GBR</td><td>45</td><td>46</td></tr>
<tr><td>IBS</td><td>53</td><td>54</td></tr>
<tr><th>TOTAL</th><th>263</th><th>240</th></tr>
</table>
</td><td>
<table>
<tr><th>POP</th><th>Females</th><th>Males</th></tr>
<tr><td>GIH</td><td>47</td><td>56</td></tr>
<tr><td>PJL</td><td>48</td><td>48</td></tr>
<tr><td>BEB</td><td>44</td><td>42</td></tr>
<tr><td>STU</td><td>47</td><td>55</td></tr>
<tr><td>ITU</td><td>43</td><td>59</td></tr>
<tr><th>TOTAL</th><th>229</th><th>260</th></tr>
</table>
</td></tr></table>

Population codes can be found [here](http://www.internationalgenome.org/faq/which-populations-are-part-your-study/).

### Step 3: Run the analysis

Run the analysis using `snakemake` once all of the raw data files are downloaded. Navigate to the top of the project directory and type the following commands:

```shell
source activate PAB_variation
snakemake
```

## What do the rules in the Snakefile do?
### Rule: Parse populations

We use a simple python script to parse the panel file into several lists with the IDs of each individual broken up by **male/female/individual** and by **super-population/population**. The resulting files are simple text-files stored in the results subdirectory of the `01_populations/` folder

### Rule: Calculate diversity

Diversity is calculated as **π**, or the average number of pairwise nucleotide differences per site.  is calculated for our purposes as **π = ( ∑ ( n<sub>i</sub> choose 2 ) ) / ( n choose 2 )**, where **n<sub>i</sub>** is the count of allele 1 in the sample and **n** is **∑ n<sub>i</sub>** (in this case, this is n<sub>1</sub> + n<sub>2</sub> as we assume all biallelic sites). Diversity is calculated at each site for each population that is supplied in a separate population file. We calculate diversity for each chromosome (chrX, chrY, chr8) for all populations at once.

### Rule: Create filter

The filter is created by merging a collection of bed files (for each chromosome) that are taken from the UCSC genome browser. Some filters and how we obtain them are as follows:
- **Coding sequence:** We filter for coding sequence by merging known genes from UCSC, RefSeq, and GENCODE. We use the following commands to accomplish this:
 ```shell
 cat file1 file2 file3 | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' | bedtools merge -i stdin > merged_genes_filter.bed
 ```
- **Repetitive elements:** We get repetitive elements from the simple-repeats track.
- **Segmental duplications:** We get segmental duplications from the Segmental Dups track.
- **Telomeres & Centromere:** We use the Gap table and define the gap types that we want in the filter (*"centromere telomere"*).
- **etc...**

We can run the *Snakefile* for any combination of these filter files by defining the combination of filters in config.yml as a list and then providing a list of the filter titles using the global variable **FILTER** at the top of the file.

### Rule: Create windows

We create windows using ```bedtools makewindows``` to create windows for the windowed analysis of diversity across the sex chromosomes. We provide information about the windows in *config.yml* in the following format:
- Window name
  - Window Size: *int*
  - Sliding: *boolean*
  - Overlap size: *int*

We can then select the window(s) that we want to use for the analysis at the top of the *Snakefile*. They can be listed using the global variable **WINDOW**.

### Rule: Split callable sites

The callable sites for the analysis can be generated using the strict whole genome mask from The 1000 Genomes Project:

- **Strict Mask:** [20141020.strict_mask.whole_genome.bed](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed)

This rule uses the following awk command to split the BED file into smaller temp files defining the callable sites for a specific chromosome (where **CHR** = ['chr8', 'chrX', 'chrY']):

```shell
cat INPUT_FILE | awk '$1 == CHR { print }' > OUTPUT_FILE
```

### Rule: Convert diversity to BED

This rule uses a simple python script to convert the output of **Rule: Calculate Diversity** to BED format. The script takes as input (w/ chrY for example) a file that is tab-delimited with the following columns:

<table>
<tr><th>Y</th><th>POS</th><th><i>pi</i></th></tr>
</table>

The output of the script is a new file in BED format with the following tab-delimited columns:

<table>
<tr><th>chrY</th><th> POS - 1 </th><th>POS</th><th><i>pi</i></th></tr>
</table>

We convert the coordinates and alter the chromosome ID's in the same script to make the data compatible with the rest of the analysis.

### Rule: Filter callable sites

We filter the callable sites files that are generated for each chromosome using `bedtools subtract`. We use the following command in order to remove the regions in our filter directly from the callable sites file and keep intervals that are callable and not in the filter. We use the following command:

```shell
bedtools subtract -a CALLABLE_SITES -b FILTER > FILTERED_CALLABLE_SITES
```

### Rule: Filter diversity by site

We filter the diversity files using `bedtools intersect`. Since we have already filtered the callable sites file, we can intersect FILTERED_CALLABLE_SITES with the diversity BED files in order to only get sites that are both callable and that pass the filter. We use the following command:

```shell
bedtools intersect -a DIVERSITY_FILE -b FILTERED_CALLABLE_SITES > FILTERED_DIVERSITY_BY_SITE
```

### Rule: Window analysis

The window analysis uses a simple python script that takes as input three files:
- General interval file (see [above](#rule-create-windows))
- Filtered callable site file (see [above](#rule-filter-callable-sites))
- Filtered diversity by site file (see [above](#rule-calculate-diversity))

These three files will be used to get a calculation of diversity in each window (including monomorphic sites). We sum up diversity in each window and called sites in each window and then simply divide our value of diversity by the number of called sites to get a value for diversity in every window across the genome.

**NOTE: The current window_analysis.py script is more efficient as long as the provided window file is not sliding windows.**

### Rule: Filter Windows by Callable Sites

The output of the window analysis script is a file with the following columns:
<table>
<tr><th>Chr</th><th>Start Position</th><th>End Position</th><th><i>Pi</i></th><th>Sites called in window</th><th>Number of variants used to calculate <i>pi</i> in window</th></tr>
</table>

This rule filters out all windows where the number of sites called in the window is less than 10% of the total sites in that window (e.g. when the window is 100kb, any window with less than 10,000 called sites will be masked out). This rule checks the condition and changes the last three columns to NA if the filter is not passed. 

### Rule: Plot Windowed Diversity

This rule integrates an R script that will take as input the windowed diversity and output a graph of diversity across the chromosome. It takes the following required options:
- Input (-i, --input)
- Output (-o, --output)
- Chromosome (-c, --chrom)

It also takes the following, optional formatting options:
- Units (--units; default is inches)
- Height (--height; default is 7.0)
- Width (--width; default is 12.0)



## TO DO
- [x] permutation tests for chrX regions
- [ ] correct for substitution rate for windowed/byRegion files
- [x] filter for windows that have < 10% called
- [ ] automate the LD analysis
- [ ] ----
