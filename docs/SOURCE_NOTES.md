# Source and reconstruction notes / 来源与重构说明

The historical source inspected locally is `4job/RNA-seq实际使用脚本/上游分析流程.docx` (via its text extraction) plus the pre-existing `4job/rnaseq_workflow` files. The document records commands for gffread, STAR index/alignment, featureCounts, and downstream R analysis. It is **not** a complete execution log. The Snakemake files were already present in `4job` as an explicitly **unrun reconstruction**; this repository corrects their toy configuration and validation and does not claim they were used on the historical study.

Representative sanitized command lines are in [HISTORICAL_COMMANDS.md](HISTORICAL_COMMANDS.md); they are evidence for design choices, not instructions to execute unchanged.

历史脚本可确认的参数包括 STAR `--sjdbOverhang 149`、`--genomeSAindexNbases 12`、`--twopassMode Basic`、`--quantMode TranscriptomeSAM GeneCounts`，以及 featureCounts `-t transcript -g gene_id -p`。公开版用虚构组名和参考文件名代替原项目标识；没有复制原 DOCX、FASTQ、比对结果或基因级结果。

Changes made for this reproducible edition / 本次修正：

- `SortedByCoordinate` replaces original unsorted BAM output.
- `-t exon -g gene_id -p --countReadPairs` replaces `-t transcript -g gene_id -p`; feature-level counting and read-pair counting therefore differ from the historical command. These are **new decisions**, not a claim about the original run.
- The synthetic fixture has its own complete `tests/config.toy.yaml` and sample sheet instead of ineffective dotted-key overrides.
- Real reference build and sample metadata still need the user's own review; the original document used apparently inconsistent reference file versions. Never combine unmatched FASTA/GFF builds.
- The earlier separate `scripts/audit.R` remains a small synthetic design/count QC demonstration; it is not a substitute for STAR/featureCounts/DESeq2.

No raw-read QC or trimming step was evidenced in the historical document. This repository does not silently attribute one to it; users should add and validate QC appropriate to their data before production analysis.
