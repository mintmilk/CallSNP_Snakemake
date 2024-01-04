#!/bin/bash
#SBATCH --job-name=tree     # 你的作业名称
#SBATCH --partition=XiaoYueHe      # 使用的分区名称
#SBATCH --nodes=1                  # 需要的节点数
#SBATCH --ntasks-per-node=1        # 每个节点的任务数
#SBATCH --cpus-per-task=24          # 每个任务的CPU核数
#SBATCH --mem=128G                  # 所需的内存
#SBATCH --output=phylogenetic_tree.out     # 输出文件
#SBATCH --mail-user=2392593414@qq.com 
#SBATCH --mail-type=ALL

source /work/apps/tools/conda/minconda3/20230202/bin/activate

conda activate bioinfo
echo "conda enter bioinfo"

# Define input, intermediate, and output paths
INPUT_VCF="combine_SNP.filtered.vcf"
FILTERED_VCF="combine_SNP.filtered.highqual.vcf"
INTERMEDIATE_PHY="converted.phy"
OUTPUT_TREE="final_tree.newick"

# Filter the VCF for high quality variants
bcftools filter -i 'QUAL > 30 && INFO/DP > 10' $INPUT_VCF -o $FILTERED_VCF
echo "VCF filtering completed."

# Convert VCF to FASTA
python3 vcf2phylip.py -i $FILTERED_VCF --output-prefix $INTERMEDIATE_PHY
echo "VCF to FASTA conversion completed."

mv converted.phy.min4.phy $INTERMEDIATE_PHY
echo "phy renamed"

raxmlHPC -s $INTERMEDIATE_PHY -n tree_result -m GTRGAMMA -p 12345

echo "Tree building completed."
