# Analyses from X Chromosome Variation Paper
**Daniel Cotter[1], Timothy Webster[2,3], Melissa Wilson[3,4]**

**2021**

1. Department of Genetics, Stanford University
2. Department of Anthropology, University of Utah
3. School of Life Sciences, Arizona State University
4. Center for Evolution and Medicine, Biodesign Institute, Arizona State University

## Analysis Steps

### Step 1: Install/Set-up

Clone this repository and make sure [conda](https://conda.io/docs/user-guide/install/index.html) is installed. Then use the provided *PAB_variation.yml* environment file to set up a new conda environment.

```shell
conda env create -f environment.yml
```

### Step 2: Get data

Data should be downloaded from the provided links and copied into the **data/** directory of the project folder.

#### Variant Data
Variant data is from The 1000 Genomes Project phase3 VCF files for chrX, chrY, & chr8:

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
conda activate chrX_variation
snakemake
```
