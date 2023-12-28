import os
import glob
from pathlib import Path

##############################环境变量修改#################################

# 定义基本目录和文件
INDEX_DIR = Path("index")

# 获取环境变量中的GATK路径或使用默认路径
GATK_PATH = os.environ.get('GATK_PATH', '/work/apps/tools/gatk/gatk-4.4.0.0/gatk')
GENOME = INDEX_DIR / "GCF_000002315.6_GRCg6a_genomic.fna"

fastq_files = sorted(glob.glob("rawdata/*_1.fq.gz"))
SAMPLES = [Path(f).stem.rsplit("_", 1)[0] for f in fastq_files]
print(SAMPLES)

##############################环境变量修改#################################

rule all:
    input:
        expand("gatk/{sample}.HC.g.vcf.gz", sample=SAMPLES)

#构建索引
rule index_genome:
    input:
        fa=GENOME
    output:
        bwa_files = [str(GENOME) + x for x in [".amb", ".ann", ".bwt", ".pac", ".sa"]],
        fai = str(GENOME) + ".fai",
        dict = Path(GENOME).with_suffix(".dict")
    log:
        "logs/index.log"
    params:
        prefix = GENOME,
        partition = config['partitions'].get('default')
    shell:
       """
       bwa index -p {params.prefix} {input.fa}
       samtools faidx {input.fa}  
       {GATK_PATH} CreateSequenceDictionary -R {input.fa} -O {output.dict} > {log} 2>&1
       """
# 使用 fastp 进行质量控制和过滤
rule fastp_qc:
    input:
        r1="rawdata/{sample}_1.fq.gz",
        r2="rawdata/{sample}_2.fq.gz"
    output:
        r1_out=temp("fastp/{sample}_1_clean.fq.gz"),
        r2_out=temp("fastp/{sample}_2_clean.fq.gz"),
        jsn="fastp/{sample}.json",
        hml="fastp/{sample}.html"
    threads: 2
    log:
        "logs/fastp_qc/{sample}.log"
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        mkdir -p fastp
        fastp -i {input.r1} -I {input.r2} -o {output.r1_out} -O {output.r2_out} -j {output.jsn} -h {output.hml} -w {threads} --length_required=50 --n_base_limit=6 --compression=6 > {log} 2>&1
        """

# 使用 BWA 进行序列比对
rule bwa_mem:
    input:
        ref = GENOME,
        r1 = "fastp/{sample}_1_clean.fq.gz",
        r2 = "fastp/{sample}_2_clean.fq.gz",
        amb = str(GENOME) + ".amb",
        ann = str(GENOME) + ".ann",
        bwt = str(GENOME) + ".bwt",
        pac = str(GENOME) + ".pac",
        sa = str(GENOME) + ".sa"
    output:
        bam =temp("bwa/{sample}.bam") 
    threads: 4
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        partition = config['partitions']['bwa_mem'],
    shell:
        """
        mkdir -p bwa
        (time bwa mem -t {threads} -M -Y -R "@RG\\tID:foo_lane\\tPL:ILLUMINA\\tLB:library\\tSM:{wildcards.sample}" {input.ref} {input.r1} {input.r2} | samtools view -Sb - > {output.bam}) 2>&1 | tee {log}
        echo "** BWA MEM done **" 2>&1 | tee -a {log}
        """

# 使用 samtools 对 BAM 文件进行排序
rule sort_bam:
    input:
        bam="bwa/{sample}.bam",
        fai= str(GENOME) + ".fai"
    output:
        sorted_bam=temp("bwa/{sample}.sorted.bam")
    threads: 4
    log:
        "logs/sort_bam/{sample}.log"
    resources:
        mem_mb=8000
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        (time samtools sort -@ {threads} -m {resources.mem_mb}M -O bam -o {output.sorted_bam} {input.bam} 2>&1; echo "** sorted raw bamfile done **") | tee {log}
        """

# 使用 GATK 进行重复序列标记
rule mark_duplicates:
    input:
        bam="bwa/{sample}.sorted.bam"
    output:
        markdup_bam=temp("bwa/{sample}.sorted.markdup.bam"),
        metrics=temp("bwa/{sample}.markdup_metrics.txt")
    threads: 4
    log:
        "logs/mark_duplicates/{sample}.log"
    params:
        validation_stringency="LENIENT",
        max_file_handles=1000,
        partition = config['partitions'].get('default')
    shell:
        """
        (time {GATK_PATH} MarkDuplicates -I {input.bam} -M {output.metrics} -O {output.markdup_bam} --VALIDATION_STRINGENCY {params.validation_stringency} --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP {params.max_file_handles} 2>&1; echo "** {wildcards.sample}.markdup.bam done **") | tee {log}
        """

# 对标记了重复序列的 BAM 文件进行索引
rule samtools_index:
    input:
        bam="bwa/{sample}.sorted.markdup.bam"
    output:
        bai=temp("bwa/{sample}.sorted.markdup.bam.bai")
    threads: 4
    log:
        "logs/samtools_index/{sample}.log"
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        (time samtools index {input.bam} 2>&1; echo "** {wildcards.sample}.sorted.markdup.bam index done **") | tee {log}
        """

# 使用 GATK 进行变异检测，并生成 GVCF 文件
rule haplotype_caller:
    input:
        bam="bwa/{sample}.sorted.markdup.bam",
        bai="bwa/{sample}.sorted.markdup.bam.bai",
        fai= str(GENOME) + ".fai",
        gatk_dict= Path(GENOME).with_suffix(".dict"),
        ref= GENOME
    output:
        "gatk/{sample}.HC.g.vcf.gz"
    threads:1
    log:
        "logs/haplotypecaller/{sample}.log"
    resources:
        mem_mb=16000
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        mkdir -p gatk
        (time {GATK_PATH} --java-options "-Xmx16G -Djava.io.tmpdir=./" HaplotypeCaller -R {input.ref} -I {input.bam} --emit-ref-confidence GVCF --native-pair-hmm-threads {threads} -O {output} 2>&1; echo "** {wildcards.sample}.HC.g.vcf.gz done **") | tee {log}
        """
