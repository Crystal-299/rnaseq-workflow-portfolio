#!/usr/bin/env python3
"""生成 RNA-seq 工作流的确定性玩具数据（自测用，不用于真实分析）。

产物（写入 tests/toy/）：
  genome.fa             —— 3 条 contig 约 1Mb（STAR 可索引的最小规模）
  genome.gff3           —— 6 个单外显子基因 g1..g6
  raw/{样本}_R{1,2}.fastq.gz —— CASE_rep1..3 / CTRL_rep1..3 双端 150bp 读段

设计：
  - seed=42 固定，可复现
  - 差异表达：CASE 高表达 g1..g3、低表达 g4..g6；CTRL 相反
  - 读段直接从基因组切出（R2 为反向互补），确保 STAR 能唯一比对

重要限制：STAR 比 bowtie2 重得多。1Mb 是能索引的底限，但索引仍需几百 MB 内存；
真实运行时用 --config star_genome_sa_index_nbases=8（小基因组必须调小，否则报错）。
"""
import gzip
import os
import random

SEED = 42
READ_LEN = 150
FRAG_LEN = 250        # 插入片段长度
HIGH = 400            # 每个高表达基因的读段对数
LOW = 40              # 每个低表达基因的读段对数

# 6 个单外显子基因坐标（放在不同的 contig 位置，保证互不重叠）
GENES = [
    # (gene_id, contig, start, end)  1-based, 闭区间
    ("g1", "chr1", 100_000, 101_000),
    ("g2", "chr1", 200_000, 201_000),
    ("g3", "chr1", 300_000, 301_000),
    ("g4", "chr2", 100_000, 101_000),
    ("g5", "chr2", 200_000, 201_000),
    ("g6", "chr3", 100_000, 101_000),
]

CONTIGS = [("chr1", 400_000), ("chr2", 300_000), ("chr3", 300_000)]

# 条件 -> 高表达基因集合（另一组为低表达）
EXPR = {"CASE": {"g1", "g2", "g3"}, "CTRL": {"g4", "g5", "g6"}}

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def make_genome(rng):
    seqs = {}
    for name, length in CONTIGS:
        seqs[name] = "".join(rng.choice("ACGT") for _ in range(length))
    return seqs


def write_fasta(seqs, path):
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def write_gff(genes, path):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        for gene_id, contig, start, end in genes:
            fh.write(f"{contig}\t.\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id}\n")
            # 单外显子基因：一个 exon，与 gene 同区间；gene_id 用于 featureCounts 计数
            fh.write(f"{contig}\t.\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene_id}.t1;Parent={gene_id}\n")
            fh.write(f"{contig}\t.\texon\t{start}\t{end}\t.\t+\t.\tID={gene_id}.t1.exon1;Parent={gene_id}.t1\n")


def sample_reads(seqs, gene, count, rng):
    """从某个基因区间切出 count 对读段。返回 [(R1, R2), ...]"""
    _, contig, start, end = gene  # start/end 为 1-based 闭区间
    seq = seqs[contig]
    # 有效起始位点（0-based，Python 切片）：保证片段落在基因内
    usable = end - start + 1 - FRAG_LEN
    if usable <= 0:
        raise ValueError(f"基因 {gene[0]} 长度不足以容纳 {FRAG_LEN}bp 片段")
    reads = []
    for _ in range(count):
        pos = start - 1 + rng.randint(0, usable - 1)   # 0-based 起点
        r1 = seq[pos:pos + READ_LEN]
        r2 = revcomp(seq[pos + FRAG_LEN - READ_LEN:pos + FRAG_LEN])
        reads.append((r1, r2))
    return reads


def write_fastq(reads, r1_path, r2_path):
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i, (r1, r2) in enumerate(reads):
            name = f"@read{i}"
            f1.write(f"{name}/1\n{r1}\n+\n{'I' * READ_LEN}\n")
            f2.write(f"{name}/2\n{r2}\n+\n{'I' * READ_LEN}\n")


def main():
    outdir = os.path.join(os.path.dirname(__file__), "toy")
    rawdir = os.path.join(outdir, "raw")
    os.makedirs(rawdir, exist_ok=True)

    rng = random.Random(SEED)
    seqs = make_genome(rng)
    write_fasta(seqs, os.path.join(outdir, "genome.fa"))
    write_gff(GENES, os.path.join(outdir, "genome.gff3"))

    for cond, hi in EXPR.items():
        for rep in (1, 2, 3):
            sample = f"{cond}_rep{rep}"
            reads = []
            for gene_id, contig, start, end in GENES:
                gene = (gene_id, contig, start, end)
                n = HIGH if gene_id in hi else LOW
                reads.extend(sample_reads(seqs, gene, n, rng))
            rng.shuffle(reads)
            write_fastq(reads,
                        os.path.join(rawdir, f"{sample}_R1.fastq.gz"),
                        os.path.join(rawdir, f"{sample}_R2.fastq.gz"))
            print(f"[OK] {sample}: {len(reads)} 对读段")

    print(f"\n玩具数据已生成到 {outdir}")
    print("运行示例：snakemake --configfile tests/config.toy.yaml -n --cores 2")


if __name__ == "__main__":
    main()
