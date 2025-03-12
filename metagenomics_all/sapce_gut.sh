path_fastq=/home/liusai/24/aging/RNA-seq
path_sam=/home/liusai/24/aging/sam
path_kraken=/home/liusai/24/aging/kraken2
path_bracken=/home/liusai/24/aging/bracken
path_megahit=/home/liusai/24/aging/MegaHit
path_prodigal=/home/liusai/24/aging/prodigal
path_humann=/home/liusai/24/aging/humann


mkdir -p ${path_sam} ${path_kraken} ${path_bracken} ${path_megahit} ${path_prodigal} ${path_humann}


source /home/liusai/miniconda3/etc/profile.d/conda.sh
conda activate rna


cat space.txt | while read id
do
    echo "开始处理样本: ${id} 时间: $(date +"%Y-%m-%d %T")"
    id=$(echo "$id" | tr -d '\r')
    bowtie2 -p 32 -x /home/liusai/index/bowtie2/human/ncbi/p14/human \
            -1 ${path_fastq}/GLDS-599_GMetagenomics_${id}_R1_filtered.fastq.gz \
            -2 ${path_fastq}/GLDS-599_GMetagenomics_${id}_R2_filtered.fastq.gz \
            -S ${path_sam}/${id}.sam --un-conc ${path_sam}/${id}.fq --very-sensitive
    rm ${path_sam}/${id}.sam 
done


cat space.txt | while read id
do
    id=$(echo "$id" | tr -d '\r')
    kraken2 --threads 16 \
            --db /home/liusai/index/kraken2/db \
            --confidence 0.05 \
            --output ${path_kraken}/${id}.output \
            --report ${path_kraken}/${id}.kreport \
            --paired \
            ${path_sam}/${id}.1.fq \
            ${path_sam}/${id}.2.fq
done


cat space.txt | while read id
do
    id=$(echo "$id" | tr -d '\r')
    bracken -d /home/liusai/index/kraken2/db \
            -i ${path_kraken}/${id}.kreport \
            -o ${path_bracken}/${id}.bracken \
            -w ${path_bracken}/${id}.bracken.report \
            -r 150 -l S
done
