
path_fastq=/home/liusai/24/old/oral/RAW
path_result=/home/liusai/24/old/oral/result
path_clean=${path_result}/clean
path_sam=${path_result}/sam
path_kraken=${path_result}/kraken2
path_bracken=${path_result}/bracken
path_megahit=${path_result}/MegaHit
path_prodigal=${path_result}/prodigal
path_humann=${path_result}/humann
path_report=${path_result}/report


mkdir -p ${path_report} ${path_clean} ${path_sam} ${path_kraken} ${path_bracken} ${path_megahit} ${path_prodigal} ${path_humann}


source /home/liusai/miniconda3/etc/profile.d/conda.sh
conda activate rna


cat oral.txt | while read id
do
  id=$(echo "$id" | tr -d '\r')
  fastp -i ${path_fastq}/${id}.1.fq.gz -o ${path_clean}/${id}_filter_R1.fq.gz \
        -I ${path_fastq}/${id}.2.fq.gz -O ${path_clean}/${id}_filter_R2.fq.gz \
        -h ${path_report}/${id}_report.html
done


cat oral.txt | while read id
do
    echo "开始: $(date +"%Y-%m-%d %T")"
    id=$(echo "$id" | tr -d '\r')
    bowtie2 -p 32 -x /home/liusai/index/bowtie2/human/ncbi/p14/human \
            -1 ${path_clean}/${id}_filter_R1.fq.gz \
            -2 ${path_clean}/${id}_filter_R2.fq.gz \
            -S ${path_sam}/${id}.sam --un-conc ${path_sam}/${id}.fq --very-sensitive
    rm ${path_sam}/${id}.sam  
done


cat oral.txt | while read id
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


cat oral.txt | while read id
do
    id=$(echo "$id" | tr -d '\r')
    bracken -d /home/liusai/index/kraken2/db \
            -i ${path_kraken}/${id}.kreport \
            -o ${path_bracken}/${id}.bracken \
            -w ${path_bracken}/${id}.bracken.report \
            -r 150 -l S
done
