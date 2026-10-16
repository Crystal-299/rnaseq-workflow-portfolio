# Sanitized historical command outline / 脱敏历史命令轮廓

The earlier local RNA-seq document contains the following command structure. Paths and sample identifiers are replaced by placeholders; HTML escaping and capitalization are normalized for readability. This is **source evidence**, not a runnable script or a verified log. Use the corrected Snakemake implementation in this repository for a new run.

```bash
gffread <annotation.gff3> -T -o <annotation.gtf>
STAR --runMode genomeGenerate --genomeDir <index_dir> \
  --genomeFastaFiles <reference.fa> --sjdbGTFfile <annotation.gtf> \
  --sjdbOverhang 149 --genomeSAindexNbases 12
STAR --runMode alignReads --readFilesCommand zcat \
  --quantMode TranscriptomeSAM GeneCounts --runThreadN 30 \
  --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None \
  --genomeDir <index_dir> --readFilesIn <sample_R1.fq.gz> <sample_R2.fq.gz> \
  --outFileNamePrefix <sample_prefix>
featureCounts -T 30 -t transcript -g gene_id -p \
  -a <annotation.gtf> -o <counts.txt> <aligned.bam>
```

The original R notes mention DESeq2/ggplot2 and PCA but do not supply a clean package-versioned executable here. The public `workflow/scripts/deseq2.R` is part of the pre-existing `4job/rnaseq_workflow` reconstruction, not a verbatim historical file. In particular, `Unsorted` and `-t transcript` above are preserved **only** as provenance; the new workflow uses sorted BAM and exon/gene paired-fragment counts. Do not compare those count matrices as though they were produced by identical parameters.
