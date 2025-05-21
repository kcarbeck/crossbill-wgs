# crossbill-wgs

Samples for this project were obtained from the https://github.com/tparchman/loxia_wgrs github repo

## Bioinformatics Pipeline

This project includes a bioinformatics pipeline for processing WGS data from 216 crossbill samples, including most representatives of Loxia curvirostra that have been recognized as subspecies, ecotypes, or other populations of interest. 

The pipeline is organized into several steps, with individual scripts for each stage located in the `bioinformatics_pipeline/` directory:

1.  **Quality Control** (`01.0.fastqc.sh`, `01.1.fastqc_missing_files.sh`, `01.2.checkFiles.sh`):
    *   Assess raw read quality using FastQC/MultiQC. MultiQC output = `multiqc_general_stats.txt`
    *   Check integrity of fastq.gz files using `gunzip -t`.
2.  **Adapter and Quality Trimming** (`02.0.adapter_removal_helper_script.py`, `02.0.concatenate_fastq.sh`, `02.1.adapter_removal.sh`, `02.2.check_AR_output.py`):
    *   Remove adapter sequences and low-quality bases using AdapterRemoval (v2.1.1)
    *   Helper script (`02.0.adapter_removal_helper_script.py`) generates commands based on sample metadata. Note: some samples with multiple fastq files were manually concatenated (`02.0.concatenate_fastq.sh`)
    *   Output `.settings` files are parsed by a python script (`02.2.check_AR_output.py`) to generate `adapter_removal_summary.csv`
3.  **Mapping** (`03.0.mapping.sh`, `03.1.sambamba_sort.sh`):
    *   Index the reference genome (`red_crossbill_24Feb2018_V5eIH.fasta`) using `bowtie2-build`.
    *   Align processed reads to the indexed reference genome using Bowtie2
    *   Alignment output is piped directly to `samtools view -bS` to generate BAM files
    *   BAM files are sorted and indexed using Sambamba (v0.7.1) using coordinate sorting
4.  **Mapping QC** (`03.2.qualimap.sh`):
    *   Evaluate the quality of the read alignments using Qualimap (v2.2.1), with the `bamqc` command
5.  **BAM File Preparation** (`04.0.mark_duplicates.sh`, `04.1.validate_bam.sh`, `04.2.index_bam.sh`):
    *   Mark duplicate reads using Picard (v2.8.2) `MarkDuplicates`  
    *   Validate BAM file integrity using picard `ValidateSamFile` 
    *   Index the processed BAM files using Picard `BuildBamIndex`.
