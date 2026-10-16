# Limitations / 局限

- The Snakemake workflow is a reconstructed implementation from the local `4job` material, not a demonstrated historical run. Its full STAR/featureCounts/DESeq2 chain was not executed here; synthetic FASTQ generation, configuration validation and a 10-job upstream DAG dry-run were checked.
- A real reference genome, matching annotation, FASTQs, biological replicates, and a Linux environment with STAR/gffread/Subread/DESeq2 are required for a real run. The toy data do not support biological conclusions.
- The historical script did not show a separate FASTQ QC/trimming command. Add and document one before using the workflow on production reads.
- The corrected exon/paired-fragment featureCounts settings differ from the historical transcript setting; gene counts should not be treated as numerically identical.

- Synthetic counts only. No source samples, raw sequencing files, annotations, figures, or findings are included.
- This checks count-table integrity and replicate structure; it does not normalize counts, detect outliers, fit a differential expression model, or validate experimental provenance.
- The original RNA-seq run was completed, but this repository's R code was newly written for the portfolio and has only been validated on its toy fixture.
