#!/bin/bash

# 激活 Conda 环境
source ~/miniconda3/etc/profile.d/conda.sh
echo "conda ok"

conda activate wga
echo "activate wga"

# 启动最新JAVA
module load tools/java/v20.0.1
echo "Java-v20.0.1 ok"

export GATK_PATH="/work/apps/tools/gatk/gatk-4.4.0.0/gatk"

# 运行 Snakemake 并使用 SLURM 提交任务
mkdir -p log
snakemake -s fastp_gvcf_v2.smk --configfile config.yaml --use-conda --cluster "sbatch --partition={params.partition} --cpus-per-task={threads} --mem={resources.mem_mb}M --output=log/slurm_{rule}.out --error=log/slurm_{rule}.err" -j 200

echo "脚本执行完成" 
curl https://api.day.app/EFoykKEWiQMTNRq7wscPQZ/CallSNP_Completed