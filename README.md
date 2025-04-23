# BRB‑seq Workflow

## 1. Demultiplex BRBseq library coming from mutiple species

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

