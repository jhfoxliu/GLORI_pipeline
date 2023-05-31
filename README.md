# GLORI_pipeline

GLORI data processing pipeline, modified from the m5C pipeline (https://github.com/SYSU-zhanglab/RNA-m5C). 

This pipeline is also suitable for eTAM. However, for eTAM, it supposes to use more memory. This is because I load every read with non-converted As into memory, so it is possible to load the whole huge fastq if the conversion rate is low. Be careful!

Some scripts are identical as the original pipeline. here I change some names for better understanding what they are for.

## Usage

### Preparation

1. Generate the metadata files based on that in https://github.com/SYSU-zhanglab/RNA-m5C
2. Generate the hisat2 indexes with the GLORI-specific script: `A2G_hisat2_index.py`
3. Perform analysis step-by-step or using the SJM file generator or using the Jupyter notebook.

### For step-by-step analysis

**Example cmd is presented in the README.md in the `notebook` folder**

1. If appliable, run `UMI-tools` or other software or scripts to extract Unique Molecular Index (UMI) before QC.
2. Run `Cutadapt` to clean the reads, if necessary run `Trimmomatic` for better clean up.
3. Make sure that you have reads which match the real strandness of RNA: for GLORI, A's are converted into G's. Hence, "forward reads" here is defined as "A-less" reads; and the "reverse reads" is defined as the "T-less" reads. If you have only *1* end of the reads, and they are "T-less" reads, please make the bases reverse-complement and quality values reverse first.
4. Run the `A2G_hisat2.py` script to perform hisat2 alignment. Use `-F` with forward reads, and `-R` with reverse reads (optional). 
5. Bowtie2 alignment is omitted. If you need it, contact Jianheng Liu (Fox).
6. Run the `pileup_genome.py` script to perform pileup.
7. Run the `format_pileups.py` script to format the pileups.
8. Repeat 1-8 for multiple samples.
9. When all pileups formations are done, copy or soft-link the formated pileup files to a new folder.
10. Prepare a SampleSheet for the samples. Format for each row: Sample[\t]Formated Pileup[\t]gene[\t]3. No empty rows should occur. the `3` here is the `A-cutoff` used, which means that 
10. Run the `m6A_caller.py` script to extract the possible m6Am sites into one CSV file.
11. Run the `evaluate_sites.py` script to add the TRUE/FALSE label for the reliability of the sites.
12. [optional] If you have replicates to merge, prepare another SampleSheet of it. Format is  Merged_name[\t]Sample[\t]Formated Pileup[\t]gene[\t]3. Only the first two columns are useful.
13. [optional] Run the `merge_replicates.py` script the merge all replicates.

### For SJM file generator (pending)

1. Make sure that you have `SGE`, `SJM` (https://github.com/StanfordBioinformatics/SJM), and `sjm-tools` (https://pypi.org/project/sjm-tools/1.0/) installed. Create a `sjm+` shortcut for `sjm` in `/usr/bin/` for compatibility.
2. Configure the environment variables in the python script.
3. Prepare SampleSheet, you may need to modify the python script to specify the content.
4. Run the python script with your SampleSheet.
5. Wait for the alignment and pileup results.
6. Analyze the pileups according to the step-by-step manual.
