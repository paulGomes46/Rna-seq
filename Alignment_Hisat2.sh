
ref_genome="hzao2"
index="/mnt/usb-TOSHIBA_EXTERNAL_USB_20230522004738F-0:0-part1/genome/${ref_genome}/"
genome="/mnt/usb-TOSHIBA_EXTERNAL_USB_20230522004738F-0:0-part1/genome/${ref_genome}/"
trimmed_sequences="/mnt/usb-TOSHIBA_EXTERNAL_USB_20230522004738F-0:0-part1/data/mydata/trimmed"
alignment="/mnt/usb-TOSHIBA_EXTERNAL_USB_20230522004738F-0:0-part1/data/mydata/alignment"


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

#Align with Hisat2
for element in ${read_prefix[@]}
do
        # be sure the name of the index is correct. I usually use genome_index_ + name of the genome
        hisat2 --dta -x $index/genome_index_${ref_genome} \ 
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
        -G $index/genome.gtf \
        -o $assembly/${filename}.gtf 
done

stringtie --merge \
-G $index/genome.gtf \
-o $assembly/merged/merged_assembly.gtf \
$assembly/merged/mergelist.txt