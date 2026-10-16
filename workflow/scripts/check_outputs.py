#!/usr/bin/env python3
"""检查 rnaseq_workflow 的关键输出是否存在且非空（运行后自检）。

用法：python workflow/scripts/check_outputs.py config/config.yaml config/samples.tsv

退出码：0 = 所有期望输出就位；非 0 = 存在缺失或空文件。
注意：只检查「文件存在 + 大小 > 0」，不验证生物学正确性（那是下游 QC 的事）。
"""
import sys
import os


def load_yaml(path):
    import yaml
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_samples(path):
    import pandas as pd
    return pd.read_csv(path, sep="\t", dtype=str)


def main():
    if len(sys.argv) != 3:
        print("用法：check_outputs.py <config.yaml> <samples.tsv>")
        sys.exit(2)

    cfg = load_yaml(sys.argv[1])
    samples = load_samples(sys.argv[2])
    out = cfg["outdirs"]

    expected = []
    # STAR 每个样本：坐标排序 BAM + ReadsPerGene.tab
    for s in samples["sample"]:
        expected.append(f"{out['star']}/{s}_Aligned.sortedByCoord.out.bam")
        expected.append(f"{out['star']}/{s}_ReadsPerGene.out.tab")
    # featureCounts 计数矩阵 + 统计摘要
    expected.append(f"{out['counts']}/counts.txt")
    expected.append(f"{out['counts']}/counts.txt.summary")

    missing, empty = [], []
    for p in expected:
        if not os.path.exists(p):
            missing.append(p)
        elif os.path.getsize(p) == 0:
            empty.append(p)

    if missing:
        print("[FAIL] 缺失输出：")
        for p in missing:
            print(f"  - {p}")
    if empty:
        print("[FAIL] 空文件：")
        for p in empty:
            print(f"  - {p}")

    if missing or empty:
        print(f"\n共 {len(missing)} 缺失、{len(empty)} 空文件。")
        sys.exit(1)

    print(f"[OK] 全部 {len(expected)} 个期望输出均存在且非空。")
    # 额外提示：DESeq2 是可选目标，不强制检查
    if not os.path.exists(f"{out['deseq2']}/deseq2_results.csv"):
        print("      提示：deseq2 结果未生成（可选目标，运行 `snakemake deseq2` 才会产生）。")


if __name__ == "__main__":
    main()
