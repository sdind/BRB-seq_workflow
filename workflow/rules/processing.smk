import pandas as pd
import os

# Import table and create libraries config with all librairies 
library_table = pd.read_csv(config["library_table"])
library_table["speciesId"] = library_table["speciesId"].astype(str)
libraries = set(library_table["library"].astype(str).tolist())
config["libraries"] = libraries

# Create a dict library - species 
lib2species = dict(set(zip(library_table["library"].astype(str), library_table["speciesId"].astype(str))))
config["lib2species"] = lib2species

print(config["lib2species"])
# Create a dic species to libs, easier to concatenate bams later
species_to_libs = library_table.groupby("speciesId")["library"].apply(list).to_dict()
config["species_to_libs"] = species_to_libs
print(species_to_libs)


# Import species info table (speciesID - genome - annotation)
species_df = pd.read_csv(config["species_table"])
species_df["speciesId"] = species_df["speciesId"].astype(str)

# Create a dict to map speciesId -> (genome_id, annotation_id)
species_dict = {}
for row in species_df.itertuples(index=False):
    species_dict[row.speciesId] = (row.genome_id, row.annotation_id)

# Create a dictionary mapping speciesId to annotation_id (easier)
species_annotation = {row.speciesId: row.annotation_id for row in species_df.itertuples(index=False)}
print(species_annotation)

rule fastqc:
    """Run FastQC on both R1 and R2 fastq files."""
    input:
        r1 = lambda wildcards: f"{config['BRBseq_data']}/{wildcards.library}_R1.fastq.gz",
        r2 = lambda wildcards: f"{config['BRBseq_data']}/{wildcards.library}_R2.fastq.gz"
    output:
        html_r1 = f"{config['snakemake_dir_path']}/results/QC/{{library}}_R1_fastqc.html",
        html_r2 = f"{config['snakemake_dir_path']}/results/QC/{{library}}_R2_fastqc.html"
    log:
        os.path.join(config['snakemake_dir_path'], "logs/QC/{library}.log")
    conda:
        "../envs/QC.yaml"
    threads: 12
    resources:
        mem_mb = 10000
    params:
        runtime = "30:00:00",
    shell:
        """
        mkdir -p {config[snakemake_dir_path]}/results/QC
        fastqc -t {threads} --outdir {config[snakemake_dir_path]}/results/QC {input.r1}
        fastqc -t {threads} --outdir {config[snakemake_dir_path]}/results/QC {input.r2}
        """
# Note: Report for R1 fastq file may contain some “red flags” because it contains barcodes/UMIs. Still, it can provide useful information on the sequencing quality of the barcodes/UMIs.



rule star_index:
    """Generate STAR genome index for each speciesId."""
    # We only have one wildcard: speciesId
    input:
        fasta = lambda wc: f"{config['genome_dir']}/{species_dict[wc.speciesId][0]}",
        gtf   = lambda wc: f"{config['annotation_dir']}/{species_dict[wc.speciesId][1]}"
    output:
        # The index folder is named "STAR_index_{speciesId}"
        index = directory(f"{config['snakemake_dir_path']}/results/index/STAR_index_{{speciesId}}")
    threads: 12
    resources:
        mem_mb = 90000
    params:
        runtime = "30:00:00"
    conda:
        "../envs/STAR.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
             --genomeDir {output.index} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --runThreadN {threads} \
             --outTmpDir {config[snakemake_dir_path]}/STAR_tmp_{wildcards.speciesId}
        """




rule star_solo:
    """Align reads using STARsolo and generate count matrices for each library."""
    input:
        # Note: the order is important: R2 (genomic information) comes before R1 (barcode/UMI)
        r2 = lambda wc: f"{config['BRBseq_data']}/{wc.library}_R2.fastq.gz",
        r1 = lambda wc: f"{config['BRBseq_data']}/{wc.library}_R1.fastq.gz",
        # Retrieve the STAR index folder corresponding to this library’s species.
        index = lambda wc: f"{config['snakemake_dir_path']}/results/index/STAR_index_{config['lib2species'][wc.library]}"
    output:
        bams = f"{config['snakemake_dir_path']}/results/star_solo/{{library}}/Aligned.sortedByCoord.out.bam",
        bam_dir = directory(f"{config['snakemake_dir_path']}/results/star_solo/{{library}}/")
    log:
        os.path.join(config['snakemake_dir_path'], "logs/star_solo/{library}.log")
    threads: 2
    resources:
        mem_mb = 35000
    conda:
        "../envs/STAR.yaml"
    params:
        # Define a unique BAM output directory for each library. (doesn't work directly with ouptut directory)
        bam_prefix = lambda wc: f"{config['snakemake_dir_path']}/results/star_solo/{wc.library}/",
        runtime = "60:00:00"
    shell:
        """
        mkdir -p {params.bam_prefix}
        STAR --runMode alignReads \
             --outSAMmapqUnique 60 \
             --runThreadN {threads} \
             --outSAMunmapped Within \
             --soloStrand Forward \
             --quantMode GeneCounts \
             --outBAMsortingThreadN {threads} \
             --genomeDir {input.index} \
             --soloType CB_UMI_Simple \
             --soloCBstart 1 \
             --soloCBlen 14 \
             --soloUMIstart 15 \
             --soloUMIlen 14 \
             --soloUMIdedup NoDedup 1MM_All \
             --soloCellFilter None \
             --soloCBwhitelist {config[whitelist]} \
             --soloBarcodeReadLength 0 \
             --soloFeatures Gene \
             --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
             --outFilterMultimapNmax 1 \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {params.bam_prefix} \
             --readFilesIn {input.r2} {input.r1}
        """

# Note: the '--outFilterMultimapNmax 1' option is recommended for removing multiple mapping reads from the output BAM
# BRB-seq: barcode length: 14, UMI length 14 (from alithea)
# from BRBseq paper: oligo-dT sequences were separated by five random non-T nucleotides (“V”), thus increasing the total UMI length to 15 nt (10 “N” + 5 “V”). This proved to be sufficient to reduce the overrepresentation of “T”-containing UMIs 
# https://25700946.fs1.hubspotusercontent-eu1.net/hubfs/25700946/MERCURIUS-BRBseq-kit-PN-10812-10813-10814-11013-User-Guide-v0.3.2.pdf
#STARsolo performs barcode error correction automatically when you supply a valid whitelist via the --soloCBwhitelist parameter. 
#--soloCBmatchWLtype Exact forces exact matches.
#--soloCBmatchWLtype 1MM allows one mismatch (if only one valid candidate exists).
#--soloCBmatchWLtype 1MM_multi allows one mismatch but marks ambiguous cases. (default)



rule matrix_to_counts:
    """Convert STARsolo output matrix to a count matrix using an R script."""
    input:
        bam_dir = lambda wc: f"{config['snakemake_dir_path']}/results/star_solo/{wc.library}/"
    output:
        counts = f"{config['snakemake_dir_path']}/results/count_matrix/{{library}}/umi.counts.txt"
    conda:
        "../envs/R.yaml"
    log:
        os.path.join(config['snakemake_dir_path'], "logs/count_matrix/{library}.log")
    threads: 1
    resources:
        mem_mb = 5000
    params:
        runtime = "03:00:00",
        R_matrix_count = "workflow/scripts/matrix_to_counts.R"
    shell:
        """
        Rscript {params.R_matrix_count} {input.bam_dir} {output.counts}
        """


rule concat_bam:
    """Concatenate all STARsolo BAM files for libraries belonging to the same species (for geneext)."""
    input:
        bams = lambda wc: expand(
            f"{config['snakemake_dir_path']}/results/star_solo/{{library}}/Aligned.sortedByCoord.out.bam",
            library=config["species_to_libs"][wc.sp])
    output:
        unsorted_bam = f"{config['snakemake_dir_path']}/results/star_solo/concat_bam/{{sp}}_all.unsorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        f"{config['snakemake_dir_path']}/logs/star_solo/concat_bam/{{sp}}.log"
    threads: 1
    resources:
        mem_mb = 10000
    params:
        runtime = "05:00:00"
    shell:
        """
        samtools cat -o {output.unsorted_bam} {input.bams}
        """



rule sort_index_bam:
    """Sort and index the concatenated BAM file for a species."""
    input:
        unsorted_bam = f"{config['snakemake_dir_path']}/results/star_solo/concat_bam/{{sp}}_all.unsorted.bam"
    output:
        sorted_bam = f"{config['snakemake_dir_path']}/results/star_solo/concat_bam/{{sp}}_all.bam",
        sorted_bam_index = f"{config['snakemake_dir_path']}/results/star_solo/concat_bam/{{sp}}_all.bam.bai"
    conda:
        "../envs/samtools.yaml"
    log:
        f"{config['snakemake_dir_path']}/logs/star_solo/concat_bam/{{sp}}.sorted.log"
    threads: 1
    resources:
        mem_mb = 10000
    params:
        runtime = "05:00:00"
    shell:
        """
        samtools sort -o {output.sorted_bam} {input.unsorted_bam}
        samtools index {output.sorted_bam}
        """



rule geneext:
    """Extend genome annotation for a species using geneext.py."""
    input:
        annotation = lambda wc: f"{config['annotation_dir']}/{species_dict[wc.sp][1]}",
        bam = f"{config['snakemake_dir_path']}/results/star_solo/concat_bam/{{sp}}_all.bam"
    output:
        extended_annotation = f"{config['snakemake_dir_path']}/results/extended_annotation/{{sp}}_extended_annotation.gtf"
    conda:
        "../envs/environment.yaml"
    log:
        os.path.join(config['snakemake_dir_path'], "logs/results/extended_annotation/{sp}.log")
    threads: 1
    resources:
        mem_mb = 10000
    params:
        runtime = "05:00:00",
        script_geneext = config['geneext_script']
    shell:
        """
        python {params.script_geneext} -g {input.annotation} -b {input.bam} -o {output.extended_annotation} --peak_perc 0 -j {threads} -inf gtf -force
        """




use rule star_index as star_index_2 with:
    input:
        fasta = lambda wc: f"{config['genome_dir']}/{species_dict[wc.speciesId][0]}",
        gtf = f"{config['snakemake_dir_path']}/results/extended_annotation/{{speciesId}}_extended_annotation.gtf"
    output:
        index = directory(f"{config['snakemake_dir_path']}/results/extended_annotation/index/STAR_index_{{speciesId}}")




use rule star_solo as star_solo_2 with:
    input:
        r2 = lambda wc: f"{config['BRBseq_data']}/{wc.library}_R2.fastq.gz",
        r1 = lambda wc: f"{config['BRBseq_data']}/{wc.library}_R1.fastq.gz",
        index = lambda wc: f"{config['snakemake_dir_path']}/results/extended_annotation/index/STAR_index_{config['lib2species'][wc.library]}"
    output:
        bams = f"{config['snakemake_dir_path']}/results/extended_annotation/star_solo/{{library}}/Aligned.sortedByCoord.out.bam",
        bam_dir = directory(f"{config['snakemake_dir_path']}/results/extended_annotation/star_solo/{{library}}/")
    log:
        os.path.join(config['snakemake_dir_path'], "logs/extended_annotation/star_solo/{library}.log")
    params:
        bam_prefix = lambda wc: f"{config['snakemake_dir_path']}/results/extended_annotation/star_solo/{wc.library}/",
        runtime = "40:00:00"

use rule matrix_to_counts as matrix_to_counts_2 with:
    input:
        bam_dir = lambda wc: f"{config['snakemake_dir_path']}/results/extended_annotation/star_solo/{wc.library}/"
    output:
        counts = f"{config['snakemake_dir_path']}/results/extended_annotation/count_matrix/{{library}}/umi.counts.txt"
    log:
        os.path.join(config['snakemake_dir_path'], "logs/extended_annotation/count_matrix/{library}.log")
