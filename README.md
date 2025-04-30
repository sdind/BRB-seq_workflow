# BRB‑seq Workflow

## 1. Demultiplexing and Species-Level Reorganization of Multi-Species Libraries

### 1.1 Demultiplex BRBseq library coming from mutiple species

As multiple species were multiplexed into a single library, we first need to **demultiplex the library per sample** and then **re-concatenate the samples coming from the same species**.

Each sample is identified by a **14 bp inline barcode** at the **start of Read 1 (R1)**, followed by a **16bp Unique Molecular Identifier (UMI)**.  **Read 2 (R2)** contains the actual RNA sequence.

To demutiplex the library, I am using cutadapt.
Cutadapt requires a FASTA file to read barcode sequences, e.g: 
```
>Pluc_33_new
CCTTGCCAGAGTGG
>Pluc_34_new
CATTCCTCTAGCCG
>Pluc_35_new
ACTTGAGCAACCGG
```




You can use this command to help convert the CSV barcode file into a FASTA file
```bash
awk -F',' '{print ">"$1"\n"$2}' barcodes.csv > barcodes.fasta
```

Then use the cutadapt command with the follwing paramters:

```
cutadapt -j 4 --no-indels -e 0.1 --action=none \
  -g ^file:barcodes.fasta \
  -o demux/{name}_R1.fastq.gz -p demux/{name}_R2.fastq.gz \
  --untrimmed-output demux/unknown_R1.fastq.gz \
  --untrimmed-paired-output demux/unknown_R2.fastq.gz \
  input_R1.fastq.gz input_R2.fastq.gz
```

**Parameters**:
  - `^`: anchors the barcode to the start of R1.
  - `--action=none`: Do not trim barcodes; just demultiplex
  - `--no-indels`: for stricter matching (no gaps) 
  - `e 0.1` allows up to 1 mismatch for 12 bp barcodes
  - Output goes to demux/{sample}_R1/R2.fastq.gz
  - Unassigned reads go to unknown_R1/R2.fastq.gz
  - `j 2`  Use 2 threads for faster processing


### 1.2 Concatenate Samples by Species

After demultiplexing, individual samples can be grouped by species using a CSV file which contain the species name in the first column and the sample name in the second column. 
E.g:
```csv
Pluchei,Pluc_33_new
Pluchei,Pluc_34_new
Pluchei,Pluc_35_new
Hogna,Hogna 54
Hogna,Hogna 55
Hogna,Hogna 58
...
```

You can use the provided script `concat_by_species.sh` to concatenate the demultiplexed FASTQ files by species as follow:

```./concat_by_species.sh <CSV_FILE> <DEMUX_DIR> <OUTPUT_DIR```

- `CSV_FILE`: CSV file mapping species (1st column) to sample names (2nd column)
- `DEMUX_DIR`: Path to directory containig the demultiplexed FASTQ files
- `OUTPUT_DIR`: Output directory for concatenated species files

---

## 2. BRB-seq Snakemake Workflow 


This Snakemake pipeline processes BRB-seq data through several stages, including quality control with FastQC, species-specific genome indexing, read alignment using STARsolo with integrated barcode error correction, BAM file concatenation and sorting, gene annotation extension with GeneExt, re-alignment against the extended annotation, and final count matrix generation.

### 2.1 Prepare Inputs

1. **Input tables**:
Before running the workflow, you must prepare two files: the BRBseq Library Table and the Species Table.
    - 1. **BRBseq Library Table** (`config["library_table"]`). 
    A CSV table containing information for all your samples.
    Each row correspond to one sample, and the following fields are required:
          - `barcode`: Barcode of each sample.
          -  `library`: Library name.
            > [!WARNING]
            > Library name (must match the raw read filename prefix)
          - `sample`: Sample name *(optional for the workflow, but essential for analysis)*.
          - `anatId`: Anatomical entity (UBERON ID) *(optional)*.
          - `stageId`: Developmental stage (UBERON ID) *(optional)*.
          - `speciesId`: NCBI Taxonomy ID of the species.

    - 2. **Species Table** (`config["species_table"]`)
    A CSV table containing information for each species. Each row should correspond to one species and must include the following fields:
          - `speciesId`: Must match the `speciesId` in the library table.
          - `genome_id`: Genome filename (located in `config["genome_dir"]`).
          - `annotation_id`: Annotation filename (located in `config["annotation_dir"]`).

<br/>

2. **Genome Annotation for GeneExt**:
  The GeneExt tool requires a specific GTF file format.
  BRAKER2 outputs are not directly compatible with GeneExt.
  A provided Python script (`scripts/reformat_gtf_for_GeneExt.py`) can be used to reformat BRAKER2 GTF files to meet GeneExt’s requirements.
  
   **Script usage**:
   ```bash
   ./reformat_gtf_for_GeneExt.py <input.gtf> <output.gtf> <removed_genes.txt>
   ```

    - The script will:
      - Ensures that gene features include a gene_id attribute.
      - Ensures that transcript features include both gene_id and transcript_id attributes.
      - Reorders attributes in exon, CDS, intron (and similar) lines so that gene_id appears before transcript_id.
      - It also removes genes without any exon, and print their id in a text file named <removed_gene.txt>.



### 2.2 Edit the Configuration File

All paths and parameters are defined in `config/config.yaml`. Key fields:

- `BRBseq_data`: Directory containing raw FASTQ files.
- `annotation_dir`: Directory with genome annotation files.
- `genome_dir`: Directory with genome FASTA files.
- `library_table`: Path to your BRB-seq library table.
- `species_table`: Path to your species table.
- `whitelist`: Path to the barcode whitelist `.txt` file.
- `snakemake_dir_path`: Root directory of the workflow.

> [!IMPORTANT]
>- The barcode whitelist (`whitelist`) is *mandatory* and is typically provided by the sequencing center.
>- Place all FASTQ files under `config["BRBseq_data"]`. Files must follow this naming pattern:
    - `[libraryname]_R1.fastq.gz`
    - `[libraryname]_R2.fastq.gz`


## 2.3 Run Snakemake

Launch the pipeline with:

```bash
snakemake --cores <num_cores> --use-conda
```

---

### Pipeline Overview

| Step | Description |
|:-----|:------------|
| **FastQC** | Quality control for R1 and R2 FASTQ files. |
| **Genome Indexing** | STAR builds species-specific genome indexes. |
| **STARsolo Alignment** | Align reads with STARsolo using genome indexes and perform barcode correction. |
| **BAM Concatenation & Sorting** | BAM files per species are concatenated, sorted, and indexed. |
| **Gene Annotation Extension** | GeneExt extends annotations based on mapped reads. |
| **STARsolo Re-alignment** | Reads are realigned using the updated annotation from GeneExt. |
| **Count Matrix Generation** | An R script converts STARsolo outputs into final count matrices. |



