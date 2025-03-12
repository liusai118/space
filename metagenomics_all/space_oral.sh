path_fastq=/home/liusai/24/space/osd572/metagenomics
path_result=/home/liusai/24/space/result/oral/metagenomics
path_sam=${path_result}/sam
path_kraken=${path_result}/kraken2
path_bracken=${path_result}/bracken
path_megahit=${path_result}/MegaHit
path_prodigal=${path_result}/prodigal
path_humann=${path_result}/humann


mkdir -p ${path_sam} ${path_kraken} ${path_bracken} ${path_megahit} ${path_prodigal} ${path_humann}


source /home/liusai/miniconda3/etc/profile.d/conda.sh
conda activate rna


cat space_oral.txt | while read id
do
    echo "开始处理样本: ${id} 时间: $(date +"%Y-%m-%d %T")"
    id=$(echo "$id" | tr -d '\r')
    bowtie2 -p 32 -x /home/liusai/index/bowtie2/human/ncbi/p14/human \
            -1 ${path_fastq}/${id}_R1_HRremoved_raw.fastq.gz \
            -2 ${path_fastq}/${id}_R2_HRremoved_raw.fastq.gz \
            -S ${path_sam}/${id}.sam --un-conc ${path_sam}/${id}.fq --very-sensitive
    rm ${path_sam}/${id}.sam 
done


cat space_oral.txt | while read id
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


cat space_oral.txt | while read id
do
    id=$(echo "$id" | tr -d '\r')
    bracken -d /home/liusai/index/kraken2/db \
            -i ${path_kraken}/${id}.kreport \
            -o ${path_bracken}/${id}.bracken \
            -w ${path_bracken}/${id}.bracken.report \
            -r 150 -l S
done
