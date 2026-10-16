# ============================================================
# bulk RNA-seq 最小工作流 —— 主入口（Snakefile）
# ------------------------------------------------------------
# 范围：bulk RNA-seq（CASE vs CTRL），双端、Illumina、plant reference
# 用法（在本目录下运行）：
#   snakemake --configfile tests/config.toy.yaml -n --cores 2
#   snakemake --configfile config/config.yaml --cores 8
#   snakemake deseq2 --configfile config/config.yaml --cores 8
# 配置：config/config.yaml；样本：config/samples.tsv
# ============================================================

import pandas as pd

# Select an explicit configuration: --configfile config/config.yaml or tests/config.toy.yaml.

# ---------- 样本表 ----------
_samples = pd.read_csv(config["samples"], sep="\t")
SAMPLES = _samples.to_dict("records")
SAMPLE_BY_NAME = {str(s["sample"]): s for s in SAMPLES}
ALL_SAMPLES = [str(s["sample"]) for s in SAMPLES]

REF = config["reference"]
OUT = config["outdirs"]
STAR_DIR = OUT["star"]
DESEQ2_DIR = OUT["deseq2"]

# ---------- 工具参数展平 ----------
# STAR 索引参数（拆自顶层键，便于玩具数据用 --config 覆盖 nbases）
STAR_INDEX_EXTRA = (
    f"--sjdbOverhang {config['star_sjdb_overhang']} "
    f"--genomeSAindexNbases {config['star_genome_sa_index_nbases']}"
)
STAR_ALIGN_EXTRA = config["star"]["align_extra"]
FC_EXTRA = config["featurecounts"]["extra"]
DESEQ2_DESIGN = config["deseq2"]["design"]
DESEQ2_CONTRAST = ",".join(config["deseq2"]["contrast"])

# ---------- 通配符 -> 路径 ----------
def sample_r1(wc):
    return SAMPLE_BY_NAME[wc.sample]["R1"]

def sample_r2(wc):
    return SAMPLE_BY_NAME[wc.sample]["R2"]

# ---------- 最终目标（DESeq2 是可选额外目标，不在此） ----------
rule all:
    input:
        expand(f"{OUT['star']}/{{s}}_Aligned.sortedByCoord.out.bam", s=ALL_SAMPLES),
        expand(f"{OUT['star']}/{{s}}_ReadsPerGene.out.tab", s=ALL_SAMPLES),
        f"{OUT['counts']}/counts.txt",

include: "workflow/rules/pipeline.smk"
