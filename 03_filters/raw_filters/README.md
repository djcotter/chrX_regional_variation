### Rule: Create filter

Here, we have provided a variety of filters used for this analysis. The filter is created by merging a collection of bed files that are taken from the UCSC genome browser.

Some filters and how we obtain them are as follows:
- **Coding sequence:** We filter for coding sequence by merging known genes from UCSC, RefSeq, and GENCODE (v28lift37). We use the following commands to accomplish this:
 ```shell
 cat file1 file2 file3 | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' | bedtools merge -i stdin > merged_genes_filter.bed
 ```
  - To buffer +/- **x** basepairs on each gene we use `bedtools slop` on the result of the above command.
- **Repetitive elements:** We get repetitive elements from the simple-repeats track.
- **Segmental duplications:** We get segmental duplications from the Segmental Dups track.
- **Telomeres & Centromere:** We use the Gap table and define the gap types that we want in the filter (*"centromere telomere"*).

We can run the *Snakefile* for any combination of these filter files by defining the combination of filters in config.yml as a list and then providing a list of the filter titles using the global variable **FILTER** at the top of the file.

--------------

We create the filter for bins that only contain regions a certain distance from genes using the following commands:
```shell
bedtools slop -i MERGED_GENES_FILE -g GENOME_FILE -b BIN_START | bedtools flank -i stdin -g GENOME_FILE -b BIN_LENGTH > DISTANCE_FROM_GENES_ONLY.bed
cat DISTANCE_FROM_GENES_ONLY.bed | sort -k1,1 -k2,2n | bedtools complement -i stdin -g GENOME_FILE > DISTANCE_FROM_GENES_INVERT.bed
```

`GENOME_FILE` points to a file with chromosome coordinate lengths that can be downloaded [here](https://github.com/arq5x/bedtools/blob/master/genomes/human.hg19.genome).

`BIN_START` refers to the start position of the distance we are filtering from genes. For example, in a 1-5kb bin, this is 1000. `BIN_LENGTH` is the length of the bin. For example, for 1-5kb the bin length is 4000.
