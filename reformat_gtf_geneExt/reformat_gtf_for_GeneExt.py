#!/usr/bin/env python3
import sys
import re

# Regular expression to extract gene_id
gene_re = re.compile(r'gene_id "([^"]+)"')

def get_gene_id(attr_str):
    """Extract gene_id from the attribute string. Returns None if not found."""
    match = gene_re.search(attr_str)
    if match:
        return match.group(1)
    return None

def collect_genes_with_exons(input_file):
    """First pass: Return a set of gene IDs that have at least one exon."""
    genes_with_exon = set()
    with open(input_file, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2].lower()
            if feature == "exon":
                gene_id = get_gene_id(fields[8])
                if gene_id:
                    genes_with_exon.add(gene_id)
    return genes_with_exon

def reorder_attributes(attr_str):
    """Reorder attributes so that gene_id comes first, then transcript_id, then others."""
    attrs = [x.strip() for x in attr_str.split(";") if x.strip()]
    gene_attr = None
    transcript_attr = None
    others = []
    for a in attrs:
        if a.startswith("gene_id"):
            gene_attr = a
        elif a.startswith("transcript_id"):
            transcript_attr = a
        else:
            others.append(a)
    new_attrs = []
    if gene_attr:
        new_attrs.append(gene_attr)
    if transcript_attr:
        new_attrs.append(transcript_attr)
    new_attrs.extend(others)
    return "; ".join(new_attrs) + ";"

def fix_gtf_line(line):
    """Fix a single GTF line (gene, transcript, or other) according to minimal requirements.
       For non-gene/transcript features, reorder attributes if both gene_id and transcript_id are present."""
    if line.startswith("#"):
        return line
    fields = line.strip().split("\t")
    if len(fields) < 9:
        return line
    feature = fields[2].lower()
    attributes = fields[8].strip()
    
    if feature == "gene":
        if "gene_id" not in attributes:
            gene_id = attributes.strip()
            fields[8] = f'gene_id "{gene_id}";'
    elif feature == "transcript":
        if "gene_id" not in attributes or "transcript_id" not in attributes:
            transcript_id = attributes.strip()
            if ".t" in transcript_id:
                gene_id = transcript_id.split(".t")[0]
            else:
                gene_id = transcript_id
            fields[8] = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
    else:
        # For features like exon, CDS, intron, start_codon, etc.
        if "gene_id" in attributes and "transcript_id" in attributes:
            fields[8] = reorder_attributes(attributes)
    return "\t".join(fields) + "\n"

def fix_gtf_file(input_file, output_file, output_removed_genes):
    genes_with_exon = collect_genes_with_exons(input_file)
    removed_genes = set()
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                fout.write(line)
                continue
            gene_id = get_gene_id(fields[8])
            # If gene_id exists and is not in the set of genes with exons, skip it.
            if gene_id and gene_id not in genes_with_exon:
                removed_genes.add(gene_id)
                continue
            fout.write(fix_gtf_line(line))

    # Write the removed gene IDs to "removed_genes.txt"
    with open(output_removed_genes, "w") as rem_f:
        for gene in sorted(removed_genes):
            rem_f.write(gene + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: reformat_gtf_for_GeneExt.py <input.gtf> <output.gtf> <removed_gene.txt>")
    input_gtf = sys.argv[1]
    output_gtf = sys.argv[2]
    output_removed_genes = sys.argv[3]
    fix_gtf_file(input_gtf, output_gtf, output_removed_genes)
