The reformat_gtf_for_GeneExt.py script is designed to correct common formatting issues in GTF files so that they are compatible with GeneExt. 

In particular, it:
- Ensures that gene features include a gene_id attribute.
- Ensures that transcript features include both gene_id and transcript_id attributes.
- Reorders attributes in exon, CDS, intron (and similar) lines so that gene_id appears before transcript_id.
This minimal formatting is required by GeneExt to correctly parse gene, transcript, and exon information.

It also removes genes without any exon, and print their id in a text file named <removed_gene.txt>.

usage: 
./reformat_gtf_for_GeneExt.py <input.gtf> <output.gtf> <removed_gene.txt>
