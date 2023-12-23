import os
import glob
from pathlib import Path

##############################环境变量修改#################################

# 定义基本目录和文件
INDEX_DIR = Path("index")

# 获取环境变量中的GATK路径或使用默认路径
GATK_PATH = os.environ.get('GATK_PATH', '/work/apps/tools/gatk/gatk-4.4.0.0/gatk')
GENOME = INDEX_DIR / "GCF_000002315.6_GRCg6a_genomic.fna"
FINAL_OUTPUT_VCF = "final/combine_SNP.filtered.vcf"

fastq_files = sorted(glob.glob("fastp/*.fq.gz"))
SAMPLES = sorted(set(os.path.basename(f).split('.', 1)[0] for f in fastq_files))

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
    params:
        prefix = GENOME,
        partition = config['partitions'].get('default')
    shell:
       """
       bwa index -p {params.prefix} {input.fa}
       samtools faidx {input.fa}  
       {GATK_PATH} CreateSequenceDictionary -R {input.fa} -O {output.dict}
       """

# 使用 BWA 进行序列比对
rule bwa_mem:
    input:
        ref = GENOME,
        r1 = "fastp/{sample}.clean.R1.fq.gz",
        r2 = "fastp/{sample}.clean.R2.fq.gz",
        amb = str(GENOME) + ".amb",
        ann = str(GENOME) + ".ann",
        bwt = str(GENOME) + ".bwt",
        pac = str(GENOME) + ".pac",
        sa = str(GENOME) + ".sa"
    output:
        bam = "bwa/{sample}.bam"
    threads: 4
    params:
        partition = config['partitions']['bwa_mem'],
    shell:
        """
        mkdir -p bwa
        time bwa mem -t {threads} -M -Y -R "@RG\\tID:foo_lane\\tPL:ILLUMINA\\tLB:library\\tSM:{wildcards.sample}" {input.ref} {input.r1} {input.r2} | samtools view -Sb - > {output.bam}
        echo "** BWA MEM done **"
        """

# 使用 samtools 对 BAM 文件进行排序
rule sort_bam:
    input:
        bam="bwa/{sample}.bam",
        fai= str(GENOME) + ".fai"
    output:
        sorted_bam="bwa/{sample}.sorted.bam"
    threads: 4
    resources:
        mem_mb=24576
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        time samtools sort -@ {threads} -m {resources.mem_mb}M -O bam -o {output.sorted_bam} {input.bam}
        echo "** sorted raw bamfile done **"
        """

# 使用 GATK 进行重复序列标记
rule mark_duplicates:
    input:
        bam="bwa/{sample}.sorted.bam"
    output:
        markdup_bam="bwa/{sample}.sorted.markdup.bam",
        metrics="bwa/{sample}.markdup_metrics.txt"
    threads: 4
    params:
        validation_stringency="LENIENT",
        max_file_handles=1000,
        partition = config['partitions'].get('default')
    shell:
        """
        time {GATK_PATH} MarkDuplicates -I {input.bam} -M {output.metrics} -O {output.markdup_bam} --VALIDATION_STRINGENCY {params.validation_stringency} --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP {params.max_file_handles}
        echo "** {wildcards.sample}.markdup.bam done **"
        """

# 对标记了重复序列的 BAM 文件进行索引
rule samtools_index:
    input:
        bam="bwa/{sample}.sorted.markdup.bam"
    output:
        bai="bwa/{sample}.sorted.markdup.bam.bai"
    threads: 4
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        time samtools index {input.bam}
        echo "** {wildcards.sample}.sorted.markdup.bam index done **"
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
    threads:2
    resources:
        mem_mb=16384
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        mkdir -p gatk
        time {GATK_PATH} --java-options "-Xmx16G -Djava.io.tmpdir=./" HaplotypeCaller -R {input.ref} -I {input.bam} --emit-ref-confidence GVCF --native-pair-hmm-threads {threads} -O {output}
        echo "** {wildcards.sample}.HC.g.vcf.gz done **"
        """

