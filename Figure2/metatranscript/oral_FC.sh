
input=/home/liusai/24/space/result/oral/metatranscript/sam 
output=/home/liusai/24/space/result/workplace/oral
source /home/liusai/miniconda3/etc/profile.d/conda.sh
conda activate rna
mkdir -p $output/fastp
mkdir -p $output/sam
mkdir -p $output/bam
cat oral_id.txt | while read id
do
    id=${id%%.*} 

 
    cat indexoral.txt | while read type
    do
      
        mkdir -p $output/FC/${type}
        mkdir -p $output/sam/${type}/${id}

      
        bowtie2 -t -p 60 -x /home/liusai/24/space/tempindex/oral/${type}/${type} \
        -1 $input/${id}.1.fq  \
        -2 $input/${id}.2.fq  \
        -S $output/sam/${type}/${id}/${id}.sam

      
        samtools view -@10 -bS $output/sam/${type}/${id}/${id}.sam > $output/FC/${type}/${id}_bowtie2.bam
    done
done
cat index.txt | while read type
do
    featureCounts -T 10 -a /home/liusai/24/space/tempindex/oral/${type}/${type}.gtf \
                  -o $output/${type}_read.count -p -t exon -g gene_id \
                  $output/FC/${type}/*.bam
done
