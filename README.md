# Intro to rnaseqbac

This pipeline integrates tools used in high-impact papers to analyze RNA-seq data in bacteria. Briefly, the workflow includes:

:small_blue_diamond: **Genome files download and indexing**: gets FASTA and GTF files from the reference *Escherichia coli* K12 MG1655, [genome assembly ASM584v2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/), and builds aligner indexes 

:small_blue_diamond: **File merging:** concatenates R1 and R2 files from different lanes per sample

:small_blue_diamond: **Quality control**: uses FastQC and MultiQC to assess sequence quality

:small_blue_diamond: **Trimming**: uses fastp to filter low-quality sequences below Q20 and sequencing adapters

:small_blue_diamond: **rRNA filtering**: employs SortMeRNA or Ribodetector to filter ribosomal RNA sequences

:small_blue_diamond: **Aligning**: maps reads against reference genome using HISAT2 or bowtie2

:small_blue_diamond: **Counting transcripts**: uses featureCounts to compute reads assigned to gene features

:small_blue_diamond: **DESeq2**: performs differential expression analysis

## How to get it

You can download the repository as follows:
```
git clone https://github.com/rayomics/rnaseqbac
```

## Exploring the repository

Once downloaded, rnaseqbac comes with the following folders/files:

- **rRNA_databases:** contains SILVA and RFAM rRNA databases to be used in SortMeRNA

- **raw_data:** includes a test dataset encompassing a subset of expression data from [BioProject PRJNA194149](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA194149/) *(McClure R et al., "Computational analysis of bacterial RNA-Seq data.", **Nucleic Acids Res**, 2013 Aug;0041(14):e140)*.

- **scripts:** contains the files:

  - *run_deseq2.R* to run DESeq2 in R

  - *pipeline_functions.sh* with all the pipeline steps to be executed in the main script

  - *slurm_wrapper.sh* allows each pipeline step to be wrapped in a separate script to be submitted to SLURM queue

- *config.yaml*: allows you to customize script configurations such as resource parameters, the genome-associated URLs, the aligner, and the rRNA filtering tool, for example

  :heavy_check_mark: If you are running the pipeline using resources from a personal computer, you can keep `resource_profile` parameter default and adjust `threads` as needed.

  :heavy_check_mark: If you are using a SLURM cluster, you must set `resource_profile = high` and adjust `threads` accordingly

- *environment.yaml*: creates the conda environment to run the analysis

- *run_pipeline.sh*: main pipeline script 

- *sample_sheet.tsv*: has all the sample file and treatment information 

## Dependencies

To run rnaseqbac, you need to have [conda installed](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). To get the minimal conda installer on Linux:

:heavy_check_mark: Download Miniconda3 installer:
  
  ```
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  ```
:heavy_check_mark: Install Miniconda3:
  
  ```
  bash Miniconda3-latest-Linux-x86_64.sh
  ```

:heavy_check_mark: Install mamba inside conda:

  ```
  conda install -n base -c conda-forge mamba conda-libmamba-solver
  ```

## Running the pipeline

### Tutorial using *E. coli* data

After adjusting *config.yaml*, you can execute the pipeline simply by running the following:

```
bash run_pipeline.sh
```

You can find the outputs of all steps in `results` folder. 

Differential expression results are in `results/deseq2`, including Volcano and PCA plots.

### Analyzing your own data

To run rnaseqbac with your own dataset, do the following:

:heavy_check_mark: Set up *config.yaml* according to your computer resources and software preferences

:heavy_check_mark: **Ensure that all path configs are pointing to your data**

- To address that, you can either edit *config.yaml* path variables (e.g. `genome_dir`, `raw_dir`) or replace *E. coli* data with yours (recommended). 

- Make sure that *raw_dir* refers to your RNA-seq raw fastq.gz data **necessarily named in CASAVA format** like so: 
  - `samplename_SX_R[1-2]_001.fastq.gz` as in `EcoliaMG1_S1_R1_001.fastq.gz`, or
  - `samplename_SX_L00[1-4]_R[1-2]_001.fastq.gz` as in `BsubtilisUV_S1_L001_R1_001.fastq.gz` when there are files from different lanes
  - **There must not have any extra _ or spaces in the filenames!** 

- If your genome-related data is deposited on NCBI, replace `fasta_url` and `gtf_url` accordingly. In case you already have the files, comment out (`#`) these lines.

- Guarantee that *genome* folder contains your genome files in FASTA and GTF format. Their filenames **must** be, respectively, *genome.fa* and *annotation.gtf*.

    **Tip:** in case your annotation generated a GFF file, you can use `gffread` from [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/) to convert it to GTF

:heavy_check_mark: Edit *sample_sheet.tsv* to contain R1 and R2 filenames as in `raw_dir` and the corresponding condition **separated by singular tabs**, as in:
    
    file_R1	file_R2	condition
    EcoliaMG1_S1_R1_001.fastq.gz  EcoliaMG1_S1_R2_001.fastq.gz  control
    EcoliaMG2_S1_R1_001.fastq.gz  EcoliaMG2_S1_R2_001.fastq.gz  control	
    EcoliaMG3_S1_R1_001.fastq.gz  EcoliaMG3_S1_R2_001.fastq.gz  control
    EcoliLb1_S1_R1_001.fastq.gz  EcoliLb1_S1_R2_001.fastq.gz  treated
    EcoliLb2_S1_R1_001.fastq.gz  EcoliLb2_S1_R2_001.fastq.gz  treated
    EcoliLb3_S1_R1_001.fastq.gz  EcoliLb3_S1_R2_001.fastq.gz  treated

:heavy_check_mark: Make sure your sample sheet matches `condition_labels` and `group_size` in *config.yaml* (`sample` rule), otherwise adjust it accordingly

:heavy_check_mark: **Run the pipeline as in the tutorial above**

## References

**FastQC:** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

**MultiQC:** Ewels P. et al. *"MultiQC: summarize analysis results for multiple tools and samples in a single report".* **Bioinformatics**, 2016 Oct 1;32(19):3047-8.

**fastp:** Chen S. et al. *"fastp: an ultra-fast all-in-one FASTQ preprocessor".* **Bioinformatics**, 2018 Sep 1;34(17):i884-i890.

**SortMeRNA:** Kopylova E. et al. *"SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data."* **Bioinformatics**, 2012 Dec 15;28(24):3211-7. 

**Ribodetector:** Deng Z. et al. *"Rapid and accurate identification of ribosomal RNA sequences via deep learning."* **Nucleic Acids Res**, 2022 Jun 10;50(10):e60.

**Bowtie 2**: Langmead B. et al. *"Fast gapped-read alignment with Bowtie 2."* **Nat Methods**, 2012 Mar 4;9(4):357-9.

**HISAT2**: Kim D. et al. *"Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype"*. **Nat Biotechnol**, 2019 Aug;37(8):907-915. 

**featureCounts**: Liao Y. et al. *"featureCounts: an efficient general purpose program for assigning sequence reads to genomic features".* **Bioinformatics**, 2014 Apr 1;30(7):923-30. 

**DESeq2**: Love M. et al. *"Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2"*. **Genome Biol**, 2014;15(12):550.



