# Bulk RNA-seq workflow | RNA-seq 工作流

## 中文

这是依据本地历史 STAR → featureCounts → R 分析脚本、以及 `4job/rnaseq_workflow` 中已有 Snakemake 重构稿整理的**可运行工作流作品集**。历史命令有记录，但这个 Snakemake 版本**不是当年实跑流程**；未上传真实 FASTQ、参考序列、表达矩阵或未发表结果。[来源与改动](docs/SOURCE_NOTES.md)。

流程：GFF3 → GTF → STAR 索引 → 双端比对/基因计数 → featureCounts 矩阵；DESeq2 差异表达与 PCA 为可选目标。保留原脚本中的 STAR 关键参数，但将 BAM 改为坐标排序、featureCounts 改为 exon/gene 的双端片段计数；这些是本次重构的修正，不冒称原流程已使用。

**下载后验证（Linux/macOS，需 Conda/Mamba）：**

```bash
conda env create -f envs/rnaseq.yaml
conda activate rnaseq
bash tests/run_toy.sh
```

`run_toy.sh` 生成虚构的 1 Mb 参考、6 个基因和 6 份双端 FASTQ，随后校验配置、干跑 DAG、执行 STAR/featureCounts 并检查 BAM、GeneCounts 和计数矩阵。可选运行 `snakemake deseq2 --configfile tests/config.toy.yaml --cores 2`。真实分析请填写 `config/config.yaml`、`config/samples.tsv` 和同一版本的 FASTA/GFF；命令为 `snakemake --configfile config/config.yaml --cores 8`。输入样本表列：`sample, condition, replicate, R1, R2`。所有真实输入必须由使用者提供。

**已验证边界：**本机生成了 6 份玩具 FASTQ，通过样本/配置静态检查和 Python 语法检查；WSL 中用 Snakemake 9.27 对上游流程干跑，解析出 10 个任务。因缺少 STAR、gffread、featureCounts 和 DESeq2，实际工具链执行仍待验证。没有把玩具结果当作真实生物学结果。[局限](LIMITATIONS.md)。仓库另保留 `scripts/audit.R`，可对 `examples/` 中的虚构 count/metadata 做独立设计与计数 QC。

## English

This portfolio repository combines commands evidenced in earlier local STAR/featureCounts/R notes with a pre-existing, **previously unrun** Snakemake reconstruction from `4job/rnaseq_workflow`. The public edition fixes its toy configuration and makes the provenance explicit; no real reads, reference, counts, or unpublished results are included.

Pipeline: GFF3 → GTF → STAR index/alignment → featureCounts; DESeq2 and PCA are optional. Coordinate-sorted BAM and exon/gene paired-fragment counting are **reconstruction changes**, not claims about the historical run. See [source notes](docs/SOURCE_NOTES.md).

On Linux/macOS, install the environment and run the synthetic check with the three commands above. Expected outputs include six sorted BAMs, six STAR GeneCounts files, and `tests/toy/results/counts/counts.txt`; optional DESeq2 is an additional target. Real input uses the five-column sample sheet and matched FASTA/GFF in `config/`. Synthetic generation, validation and a 10-job Snakemake upstream dry-run passed; the external-tool chain was **not** executed here. See [limitations](LIMITATIONS.md).
