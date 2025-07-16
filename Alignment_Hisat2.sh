
ref_genome="hzao2"
index="/home/paul/rna_seq/genome_index/sinensis_v2"
genome="/home/paul/rna_seq/genome_index/sinensis_v2"
trimmed_sequences="/home/paul/rna_seq/data/my_data/trimmed_sequences"
alignment="/home/paul/rna_seq/data/my_data/alignment"


#Récuper la liste de tous les reads
for file in $trimmed_sequences/*trimmed*
do
        sample_name=$(basename  "$file" | cut -c1-6)
        if [[ ! ${read_prefix[@]} =~ $sample_name ]]
        then
                read_prefix+=($sample_name)
                echo $sample_name
        fi
done

echo $read_prefix

#Align with Hisat2
for element in ${read_prefix[@]}
do
        # be sure the name of the index is correct. I usually use genome_index_ + name of the genome
        echo "alignment : $element"
        hisat2 --dta -x $index \
        -1 $trimmed_sequences/${element}_1_trimmed_PE.fastq.gz \
        -2 $trimmed_sequences/${element}_2_trimmed_PE.fastq.gz \
        --fr \
        -S $alignment/${element}.sam \
        --summary-file $alignment/log/${element}_alignment_summary.txt
        samtools sort \
        -o $alignment/${element}.bam \
        $alignment/${element}.sam
        #remove the .bam files to save space
        rm $alignment/${element}.sam
done

# Reads assembly 
for element in $alignment/*.bam
do
        filename=$(basename -- $element .bam)
        stringtie $alignment/${filename}.bam \
        -l $filename \
        -G $index/*.gtf \
        -o $assembly/${filename}.gtf 
done

stringtie --merge \
-G $index/genome.gtf \
-o $assembly/merged/merged_assembly.gtf \
$assembly/merged/mergelist.txt