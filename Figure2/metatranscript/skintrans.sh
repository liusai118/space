
path_fastq=/home/liusai/24/space/osd572/metatranscriptomics
path_result=/home/liusai/24/space/result/skin/metatranscript
path_sam=${path_result}/sam
path_kraken=${path_result}/kraken2
path_bracken=${path_result}/bracken
path_mega=${path_result}/MegaHit
path_msalmon=${path_result}/msalmon
path_prodigal=${path_result}/prodigal
path_humann=${path_result}/humann
path_output=${path_prodigal}
path_input=${path_mega}
error_log="${path_output}/prodigal.txt"
eggno_data_dir=/home/liusai/index/eggno
NR=${path_result}/NR

mkdir -p ${path_result} ${path_sam} ${path_kraken} ${path_bracken} ${path_prodigal}
mkdir -p ${path_output}/gene_fa ${path_output}/gene_gff ${path_output}/fullgene_fa ${NR}
mkdir -p ${path_mega} ${path_result}/summ_table

> ${error_log}


source /home/liusai/miniconda3/etc/profile.d/conda.sh
conda activate rna

cat skin_id.txt | while read id; do
    echo ": $(date +"%Y-%m-%d %T")"
    id=$(echo "$id" | tr -d '\r')
    bowtie2 -p 32 -x /home/liusai/index/bowtie2/human/ncbi/p14/human \
    -1 ${path_fastq}/GLDS-564_metatranscriptomics_metaT_${id}_R1_HRremoved_raw.fastq.gz \
    -2 ${path_fastq}/GLDS-564_metatranscriptomics_metaT_${id}_R2_HRremoved_raw.fastq.gz \
    -S ${path_sam}/${id}.sam --un-conc ${path_sam}/${id}.fq --very-sensitive
    rm ${path_sam}/${id}.sam
done


cat skin_id.txt | while read id; do
    id=$(echo "$id" | tr -d '\r')
    kraken2 --threads 16 \
    --db /home/liusai/index/kraken2/db \
    --confidence 0.05 \
    --output ${path_kraken}/${id}.output \
    --report ${path_kraken}/${id}.kreport \
    --paired ${path_sam}/${id}.1.fq ${path_sam}/${id}.2.fq
done


cat skin_id.txt | while read id; do
    id=$(echo "$id" | tr -d '\r')
    bracken -d /home/liusai/index/kraken2/db \
    -i ${path_kraken}/${id}.kreport \
    -o ${path_bracken}/${id}.bracken \
    -w ${path_bracken}/${id}.bracken.report \
    -r 150 -l S
done


cat skin_id.txt | while read id; do
    echo "Start processing ${id}: $(date +"%Y-%m-%d %T")"
    id=$(echo "$id" | tr -d '\r')
    output_dir="${path_mega}/${id}"
    
    if [ -d "${output_dir}" ]; then
        echo "Output directory ${output_dir} already exists. Skipping ${id}."
        continue
    fi

    megahit -t 16 -1 ${path_sam}/${id}.1.fq -2 ${path_sam}/${id}.2.fq -o ${output_dir} --out-prefix ${id} --min-contig-len 135

    if [ $? -ne 0 ]; then
        echo "Error processing ${id}: $(date +"%Y-%m-%d %T")" >> ${error_log}
    else
        echo "Finished processing ${id}: $(date +"%Y-%m-%d %T")"
    fi
done

cat skin_id.txt | while IFS= read -r id || [ -n "$id" ]; do
    echo "Start processing ${id}: $(date +"%Y-%m-%d %T")"
    id=$(echo "$id" | tr -d '\r')

    input_file="${path_input}/${id}/${id}.contigs.fa"
    gene_fa_file="${path_output}/gene_fa/${id}.gene.fa"
    gene_gff_file="${path_output}/gene_gff/${id}.gene.gff"
    fullid_file="${path_output}/gene_fa/${id}.fullid"
    fullgene_fa_file="${path_output}/fullgene_fa/${id}.gene.fa"

    if [ -f "${gene_fa_file}" ]; then
        echo "Gene fasta file ${gene_fa_file} already exists, skipping ${id}"
        continue
    fi

    if [ ! -f "${input_file}" ]; then
        echo "Input file ${input_file} does not exist, skipping ${id}" | tee -a ${error_log}
        continue
    fi

    prodigal -i "${input_file}" -d "${gene_fa_file}" -o "${gene_gff_file}" -p meta -f gff || {
        echo "Prodigal failed for ${id}" | tee -a ${error_log}
        continue
    }

    if [ -f "${gene_fa_file}" ]; then
        grep 'partial=00' "${gene_fa_file}" | cut -f1 -d ' ' | sed 's/>//' > "${fullid_file}"
        seqkit grep -f "${fullid_file}" "${gene_fa_file}" > "${fullgene_fa_file}" || {
            echo "Seqkit grep failed for ${id}" | tee -a ${error_log}
        }
    else
        echo "Gene fasta file ${gene_fa_file} does not exist, skipping ${id}" | tee -a ${error_log}
    fi
done


cat ${path_prodigal}/gene_fa/*.gene.fa > ${path_prodigal}/all.gene.fa
cat ${path_prodigal}/fullgene_fa/*.gene.fa > ${path_prodigal}/all.fullgene.fa


mkdir -p ${NR}
mmseqs easy-linclust ${path_prodigal}/all.gene.fa ${NR}/nucleotide mmseqs_tmp --min-seq-id 0.9 -c 0.9 --cov-mode 1 --threads 8
mv ${NR}/nucleotide_rep_seq.fasta ${NR}/nucleotide.fa
seqkit translate ${NR}/nucleotide.fa > ${NR}/protein.fa
sed 's/\*//g' ${NR}/protein.fa > ${NR}/protein_no_end.fa


salmon index -t ${NR}/nucleotide.fa -p 4 -i ${path_msalmon}/index
salmon quantmerge --quants ${path_msalmon}/quant/*.quant -o ${path_msalmon}/gene.TPM
salmon quantmerge --quants ${path_msalmon}/quant/*.quant --column NumReads -o ${path_msalmon}/gene.count


emapper.py --no_annot --no_file_comments --override --data_dir ${eggno_data_dir} -i ${NR}/protein_no_end.fa --cpu 16 -m diamond -o ${NR}/protein
emapper.py --annotate_hits_table ${NR}/protein.emapper.seed_orthologs --data_dir ${eggno_data_dir} --cpu 8 --no_file_comments --override -o ${NR}/output

sed '1 i Name\tortholog\tevalue\tscore\ttaxonomic\tprotein\tGO\tEC\tKO\tPathway\tModule\tReaction\trclass\tBRITE\tTC\tCAZy\tBiGG\ttax_scope\tOG\tbestOG\tCOG\tdescription' \
    ${NR}/output.emapper.annotations > ${NR}/eggnog_anno_output


python /home/liusai/index/kraken2/script/summarizeAbundance.py -i ${path_msalmon}/gene.count -m ${NR}/eggnog_anno_output -c '9,16,21' -s ',+,+*' -n raw -o ${path_result}/summ_table/eggnog
sed -i 's/^ko://' ${path_result}/summ_table/eggnog.KO.raw.txt
