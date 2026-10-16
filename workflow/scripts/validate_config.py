#!/usr/bin/env python3
"""校验 rnaseq_workflow 的 config.yaml 与 samples.tsv 是否自洽、可运行。

用法：python workflow/scripts/validate_config.py config/config.yaml config/samples.tsv

退出码：0 = 通过（含警告）；非 0 = 致命错误。
仅做「静态」校验（字段齐全、样本唯一、条件成对、版本/参考位点），
不检查软件是否已安装、文件是否真实存在——那属于 check_outputs.py 的职责。
"""
import sys
import os

# 期望样本表列（双端 RNA-seq：需要 R1/R2）
REQUIRED_SAMPLE_COLUMNS = ["sample", "condition", "replicate", "R1", "R2"]

# 参考位点（缺了会导致规则无法解析输入/输出）
REQUIRED_REF_KEYS = ["fasta", "gff", "gtf", "star_index"]


def load_yaml(path):
    try:
        import yaml
    except ImportError:
        print("[ERROR] 未安装 PyYAML：pip install pyyaml（或使用 envs/rnaseq.yaml）")
        sys.exit(1)
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_samples(path):
    import pandas as pd
    if not os.path.exists(path):
        print(f"[ERROR] 样本表不存在：{path}")
        sys.exit(1)
    return pd.read_csv(path, sep="\t", dtype=str)


def main():
    if len(sys.argv) != 3:
        print("用法：validate_config.py <config.yaml> <samples.tsv>")
        sys.exit(2)

    cfg_path, samples_path = sys.argv[1], sys.argv[2]
    cfg = load_yaml(cfg_path)
    samples = load_samples(samples_path)

    errors, warnings = [], []

    # ---- 1. samples.tsv 结构 ----
    missing_cols = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in samples.columns]
    if missing_cols:
        errors.append(f"样本表缺少列：{missing_cols}")
        # 缺列则后续检查无意义，直接返回
        print("\n".join(errors))
        sys.exit(1)

    if samples["sample"].duplicated().any():
        dup = samples["sample"][samples["sample"].duplicated()].tolist()
        errors.append(f"样本名重复：{dup}")

    if samples.duplicated(subset=["condition", "replicate"]).any():
        errors.append("同一 condition 内 replicate 编号重复")

    # 每个 condition 至少 2 个重复（DESeq2 需要生物学重复才能估计离散度）
    for cond, grp in samples.groupby("condition"):
        if len(grp) < 2:
            errors.append(f"condition={cond} 只有 {len(grp)} 个重复，DESeq2 需要 ≥2")

    # 比对条件必须成对出现（有处理组就要有对照组）
    conditions = set(samples["condition"].dropna())
    if len(conditions) < 2:
        errors.append(f"condition 种类少于 2：{conditions}（差异分析需要对照）")

    # ---- 2. config.yaml 结构 ----
    for key in ["reference", "versions", "threads", "mem_mb", "outdirs", "star"]:
        if key not in cfg:
            errors.append(f"config 缺少顶层键：{key}")

    if "reference" in cfg:
        ref = cfg["reference"]
        for rk in REQUIRED_REF_KEYS:
            if rk not in ref:
                errors.append(f"reference 缺少键：{rk}")

    # 版本不允许 latest（固定可复现）
    for tool, ver in cfg.get("versions", {}).items():
        if not ver or ver.lower() in ("latest", "最新"):
            errors.append(f"versions.{tool} 为空或 latest，必须固定版本")

    # ---- 3. STAR 参数合理性（仅静态告警，不作为致命错误） ----
    overhang = cfg.get("star_sjdb_overhang")
    if overhang is None:
        warnings.append("未设置 star_sjdb_overhang（读长-1）；真实读长 150 应为 149")
    nbases = cfg.get("star_genome_sa_index_nbases")
    if nbases is not None and int(nbases) < 10:
        warnings.append(
            f"star_genome_sa_index_nbases={nbases} 较小（真实 135Mb 基因组用 12；"
            "小基因组/玩具数据才需调小，见 README）"
        )

    # ---- 4. 参考位点是否为占位符（resources/... 需用户填充真实文件） ----
    if "reference" in cfg:
        ref = cfg["reference"]
        for rk in ["fasta", "gff"]:
            p = ref.get(rk, "")
            if p.startswith("resources/") and not os.path.exists(p):
                warnings.append(f"参考文件为占位符（尚未填充真实文件）：{rk} = {p}")

    # ---- 输出 ----
    for w in warnings:
        print(f"[WARN] {w}")
    if errors:
        for e in errors:
            print(f"[ERROR] {e}")
        print(f"\n共 {len(errors)} 个致命错误，{len(warnings)} 个警告。")
        sys.exit(1)

    print(f"[OK] config 与 samples 校验通过（{len(warnings)} 个警告）。")
    print(f"     样本数={len(samples)}，condition={sorted(conditions)}")


if __name__ == "__main__":
    main()
