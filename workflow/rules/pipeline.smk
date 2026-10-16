# ============================================================
# 规则集：gff2gtf -> star_index -> star_align -> featurecounts -> (deseq2 可选)
# 每条规则均含 input/output/log/threads/benchmark/resources
# ============================================================

# ---- 1. GFF -> GTF（gffread） ----
rule gff2gtf:
    input:
        REF["gff"],
    output:
        REF["gtf"],
    log:
        "logs/gff2gtf/gff2gtf.log",
    threads: config["threads"]["gff2gtf"]
    resources:
        mem_mb = config["mem_mb"]["gff2gtf"],
    benchmark:
        "benchmark/gff2gtf/gff2gtf.tsv",
    shell:
        "gffread {input} -T -o {output} 2> {log}"

# ---- 2. STAR 索引（genomeGenerate） ----
rule star_index:
    input:
        fasta = REF["fasta"],
        gtf   = REF["gtf"],
    output:
        idx = directory(REF["star_index"]),
    log:
        "logs/star_index/star_index.log",
    threads: config["threads"]["star_index"]
    resources:
        mem_mb = config["mem_mb"]["star_index"],
    benchmark:
        "benchmark/star_index/star_index.tsv",
    shell:
        "mkdir -p {output.idx} && "
        "STAR --runMode genomeGenerate --genomeDir {output.idx} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "{STAR_INDEX_EXTRA} --runThreadN {threads} 2> {log}"

# ---- 3. STAR 比对（alignReads） ----
rule star_align:
    input:
        r1  = sample_r1,
        r2  = sample_r2,
        idx = REF["star_index"],
    output:
        bam    = f"{OUT['star']}/{{sample}}_Aligned.sortedByCoord.out.bam",
        counts = f"{OUT['star']}/{{sample}}_ReadsPerGene.out.tab",
    log:
        f"logs/star_align/{{sample}}.log",
    threads: config["threads"]["star_align"]
    resources:
        mem_mb = config["mem_mb"]["star_align"],
    benchmark:
        f"benchmark/star_align/{{sample}}.tsv",
    shell:
        "STAR --runMode alignReads --readFilesCommand zcat "
        "{STAR_ALIGN_EXTRA} --runThreadN {threads} "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.r1} {input.r2} "
        "--outFileNamePrefix {STAR_DIR}/{wildcards.sample}_ 2> {log}"

# ---- 4. featureCounts（单个计数矩阵，一次吃所有 BAM） ----
rule featurecounts:
    input:
        bams = expand(f"{OUT['star']}/{{s}}_Aligned.sortedByCoord.out.bam", s=ALL_SAMPLES),
        gtf  = REF["gtf"],
    output:
        counts  = f"{OUT['counts']}/counts.txt",
        summary = f"{OUT['counts']}/counts.txt.summary",
    log:
        "logs/featurecounts/featurecounts.log",
    threads: config["threads"]["featurecounts"]
    resources:
        mem_mb = config["mem_mb"]["featurecounts"],
    benchmark:
        "benchmark/featurecounts/featurecounts.tsv",
    shell:
        "featureCounts -T {threads} {FC_EXTRA} "
        "-a {input.gtf} -o {output.counts} {input.bams} 2> {log}"

# ---- 5. DESeq2（可选，手动 snakemake deseq2） ----
rule deseq2:
    input:
        counts  = f"{OUT['counts']}/counts.txt",
        samples = config["samples"],
    output:
        results = f"{OUT['deseq2']}/deseq2_results.csv",
        pca     = f"{OUT['deseq2']}/pca.pdf",
    log:
        "logs/deseq2/deseq2.log",
    threads: 1
    benchmark:
        "benchmark/deseq2/deseq2.tsv",
    shell:
        "Rscript workflow/scripts/deseq2.R "
        "--counts {input.counts} --samples {input.samples} "
        "--design '{DESEQ2_DESIGN}' --contrast {DESEQ2_CONTRAST} "
        "--outdir {DESEQ2_DIR} 2> {log}"
