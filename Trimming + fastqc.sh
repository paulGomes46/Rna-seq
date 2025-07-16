

sequences="/home/paul/rna_seq/data/ramsey/sequences"
trimmed_sequences="/home/paul/rna_seq/data/ramsey/trimmed_sequences"
fastqc="/home/paul/rna_seq/data/ramsey/fastqc"

# preliminary fastqc
for file in $sequences/*
do
        fastqc $file -o $fastqc/before_trimming
done



#Récuper la liste de tout les reads
for file in $sequences/*
do
        sample_name=$(basename  "$file" | cut -c1-10)
        if [[ ! ${read_prefix[@]} =~ $sample_name ]]
        then
                read_prefix+=($sample_name)
        fi
done

#loop sur tout les elements 
for file in ${read_prefix[@]}
do
        trimmomatic \
        PE \
        -phred33 \
        $sequences/${file}_1.fastq \
        $sequences/${file}_2.fastq \
        $trimmed_sequences/${file}_1_trimmed_PE.fastq.gz \
        $trimmed_sequences/${file}_1_trimmed_SE.fastq.gz \
        $trimmed_sequences/${file}_2_trimmed_PE.fastq.gz \
        $trimmed_sequences/${file}_2_trimmed_SE.fastq.gz \
        ILLUMINACLIP:/home/paul/miniconda3/envs/kallisto/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
        #-trimlog $trimmed_sequences/${file}_trimmomatic_log
done

# delete SE files
for file in $trimmed_sequences/*_SE.fastq.gz
do
	rm $file
done


# after-trimming fastqc
for file in $trimmed_sequences/*trimmed*
do
        fastqc $file -o $fastqc/after_trimming
done