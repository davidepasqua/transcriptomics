####### set input directory 

INP_DIR=""



################################## FASTP ##########################################

# Loop through all the forward files in the directory

for forward_file in $INP_DIR/*_1.fastq.gz; do
    # Extract the corresponding reverse file name
    reverse_file=${forward_file/*_1_*/*_2_*}
    # Run fastp on the paired-end reads
    echo "Running fastp analysis for files $forward_file and $reverse_file"
    fastp -i $forward_file -I $reverse_file --detect_adapter_for_pe
    # Remove fastp output files
    rm fastp.json fastp.html
done


###################################### CUTADAPT #################################################

# Create output folders for processed forward and reverse reads after cutadapt

mkdir -p $INP_DIR/processed_data/after_cutadapt

# Loop through all the forward files in the directory

for file_forward in $INP_DIR/*_1.fastq.gz; do
    # Extract the corresponding reverse file name
    file_reverse=${file_forward/_1./_2.}
    # Set output file names
    basename_forward=$(basename $file_forward)
    basename_reverse=$(basename $file_reverse)
    out_for=$INP_DIR/processed_data/after_cutadapt/$basename_forward
    out_rev=$INP_DIR/processed_data/after_cutadapt/$basename_reverse
    # Run cutadapt to remove adapter sequences from reads
    echo "Running cutadapt for $basename_forward and $basename_reverse"
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o $out_for -p $out_rev \
        --cores 10 \
        $file_forward $file_reverse
done



########### TRIMMOMATIC ###########

# Create output directory for Trimmomatic output

mkdir -p $INP_DIR/processed_data/Trimmomatic_output

# Set output directory variable

out_dir="$INP_DIR/processed_data/Trimmomatic_output"

# Loop through each forward read file in the directory

for file_forward in `ls $INP_DIR/processed_data/after_cutadapt/*_1.fastq.gz`; do

    # Create name for reverse read file by replacing "_1." with "_2." in the forward read file name
    file_reverse=${file_forward/_1./_2.}

    # Set output file names
    basename_f=$(basename $file_forward)
    basename_r=$(basename $file_reverse)
    for_paired=$out_dir/${basename_f/_1/_1_PAIRED}
    rev_paired=$out_dir/${basename_r/_2/_2_PAIRED}
    for_unpaired=$out_dir/${basename_f/_1/_1_UNPAIRED}
    rev_unpaired=$out_dir/${basename_r/_2/_2_UNPAIRED}

    # Print output file names to the console
    echo $for_paired

    # Run Trimmomatic with specified parameters

    trimmomatic PE -phred33 -threads 10 $file_forward $file_reverse $for_paired $for_unpaired $rev_paired $rev_unpaired TRAILING:30 LEADING:30 HEADCROP:20 MINLEN:75 AVGQUAL:30
done

########### HISAT2 ###########

# Create output directory for HISAT2 output

mkdir -p $INP_DIR/processed_data/Hisat2_output

# Set output directory variable

out_dir="$INP_DIR/processed_data/Hisat2_output"

# Loop through each paired forward read file in the directory
for file_forward in `ls $INP_DIR/processed_data/Trimmomatic_output/*_1_PAI.fastq.gz`; do

    # Remove "_1_PAI.fastq.gz" from forward read file name to get file ID
    file_id=${file_forward/_1_PAI.fastq.gz/}
    file_id_basename=$(basename $file_id)

    # Set output file name for SAM file
    out_file_sam=$out_dir/$file_id_basename".sam"

    # Set input file names for HISAT2 alignment
    for_pai=$file_forward
    rev_pai=${file_forward/_1_PAI.fastq.gz/_2_PAI.fastq.gz}
    for_unp=${file_forward/_1_PAI.fastq.gz/_1_UNP.fastq.gz}
    rev_unp=${file_forward/_1_PAI.fastq.gz/_2_UNP.fastq.gz}

    # Run HISAT2 alignment with specified parameters
    hisat2 -x /g100/home/userexternal/dpasquan/human_genome/genome -1 $for_pai -2 $rev_pai -U $for_unp,$rev_unp -S $out_file_sam
done

########### FEATURECOUNTS ###########

# Set path to the GTF file containing gene information

genes=/g100/home/userexternal/dpasquan/human_genome/genome.gtf

# Create output directory for featureCounts output

mkdir -p $INP_DIR/processed_data/FeatureCoun
ts_output

# featureCounts for creating count matrix 
featureCounts -T 5 -p -B -t exon -g gene_id -a $genes -o /g100/home/userexternal
/dpasquan/PRJNA867309/processed_data/FeatureCounts_output/count_matrix /g100/hom
e/userexternal/dpasquan/PRJNA867309/processed_data/Hisat2_output/*.sam
