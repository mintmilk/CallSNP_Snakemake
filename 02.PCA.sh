#!/bin/bash
#SBATCH --job-name=PCA     # 你的作业名称
#SBATCH --partition=YongDingHe      # 使用的分区名称
#SBATCH --nodes=1                  # 需要的节点数
#SBATCH --ntasks-per-node=1        # 每个节点的任务数
#SBATCH --cpus-per-task=4         # 每个任务的CPU核数
#SBATCH --output=PCA.out     # 输出文件
#SBATCH --mail-user=2392593414@qq.com 
#SBATCH --mail-type=ALL

source /work/home/zhgroup02/miniconda3/bin/activate

conda activate wga
echo "conda enter wga"

# 设置vcf文件路径
VCF_FILE="feisha_IDnew_SNP_filtered_highqual.vcf"

# 使用plink进行PCA
plink \
  --vcf $VCF_FILE \
  --allow-extra-chr \
  --pca 3\
  --out PCA_output

echo "PCA analysis finished!"

Rscript PCA_ggplot2.R PCA_output.eigenvec