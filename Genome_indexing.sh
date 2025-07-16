##### 1. SET VARIABLES #####


####### Preparation of annotation files ########

#Convertion gff to gtf
agat_convert_sp_gff2gtf.pl --gff citrus_sinensis.gff -o citrus_sinensis.gtf


#Extract informations from gtf

hisat2_extract_splice_sites.py \
/home/paul/Documents/bioinfo/training/ncbi_dataset/data/GTF/genome.gtf \
> genome.ss

hisat2_extract_exons.py \
/home/paul/Documents/bioinfo/training/ncbi_dataset/data/GTF/genome.gtf \
> genome.exons


# Hisat2

hisat2-build -p 24 \
--ss $index/${ref_genome}.ss \
--exon $index/${ref_genome}.exons \
$genome/genome_${ref_genome}* \
genome_index_${ref_genome}



# Kallisto

kallisto index \
-i sinensis_V2.idx \ 
reference_transcriptome.fasta