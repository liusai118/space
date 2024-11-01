# 定义路径变量
input=/path/to/your/input/files  # 修改为实际输入文件路径
output=/home/liusai/24/space/result/workplace/oral

# 激活 Conda 环境
source /home/liusai/miniconda3/etc/profile.d/conda.sh
conda activate rna

# 创建通用输出目录
mkdir -p $output/fastp
mkdir -p $output/sam
mkdir -p $output/bam

# 循环处理样本
cat oral_id.txt | while read id
do
    id=${id%%.*}  # 去除文件扩展名

    # 循环处理类型
    cat indexoral.txt | while read type
    do
        # 创建每个类型和样本的输出目录
        mkdir -p $output/FC/${type}
        mkdir -p $output/sam/${type}/${id}

        # Bowtie2 比对
        bowtie2 -t -p 60 -x /home/liusai/24/space/tempindex/oral/${type}/${type} \
        -1 $input/${id}.1.fq  \
        -2 $input/${id}.2.fq  \
        -S $output/sam/${type}/${id}/${id}.sam

        # SAM 转 BAM
        samtools view -@10 -bS $output/sam/${type}/${id}/${id}.sam > $output/FC/${type}/${id}_bowtie2.bam
    done
done

# 使用 featureCounts 计数
cat index.txt | while read type
do
    featureCounts -T 10 -a /home/liusai/24/space/tempindex/oral/${type}/${type}.gtf \
                  -o $output/${type}_read.count -p -t exon -g gene_id \
                  $output/FC/${type}/*.bam
done
