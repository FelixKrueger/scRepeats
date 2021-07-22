# scRepeats
Custom script(s) for the processing of files for single-cell repeat analysis


## Repeat alignments of single cells - individual FastQ files

...were carried out using a multi-step approach using on lists of interesting cells based on a clustering analysis in Seurat. Due to limitations of how many filehandles a Linux machine can typically read from or write to concomitantly, this approach works for cell numbers up to several hundreds (up to 500?). For cell numbers higher than that, I might add an additional section further below (which works with one big matrix instead of many individual FastQ files).

#### Extracting FastQ files from possorted_bam.bam files (produced by CellRanger)

To get to FastQ files for individual single cells belonging to various clusters, we first used the cell barcode files from Seurat as input file to [subset10Xbam](https://github.com/s-andrews/subset10xbam) with the option `--add_barcode` as well as the `possorted_bam.bam` file from CellRanger to get subsetted BAM files for sets of interesting cell barcodes.

#### Example cell barcode file
```
AAAGTCCAGCACGTCC-1	4
AACAGGGAGTGCCGAA-1	4
AACCACATCCGCTGTT-1	4
AACGAAAAGATTGTGA-1	4
AACGTCATCTCTGCTG-1	4
...
```
This produces a single BAM file containing all alignments from interesting cells in the cell barcode file (with the CB: cell barcode in the readID).

 
#### BAM to FastQ conversion
... was accomplished using:

```
samtools fastq cells_of_interest.bam > cells_of_interest.fastq
```

This FastQ file was then split into several FastQ files for individual cells, based on the CB cell barcode in the readID (using the attached script `split_FastQ_file_by_CellRanger_barcode.py`, please note the limitation to a few hundred cells mentioned above).

These individual FastQ files were then used as input files for the various repeat alignment scripts. These scripts are specific to the structure of the cluster at the Babraham Institute, but the location of the N-padded repeat genomes can be adapted to work elsewhere.

A typcial usage could be:

```
for i in *fastq.gz; do human_LTR_repeat_family_analysis_elements_of_interest_only.pl --bowtie2 $i; done
```

Once all repeat alignments have finished, the script `make_repeat_family_summary.pl` produces a single summary matrix from all `repeat_family_report.txt` files.
