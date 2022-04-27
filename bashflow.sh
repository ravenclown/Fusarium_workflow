mkdir -p aligned aligned_files
mkdir -p bam_files
mkdir -p index

echo "Building index"

#bowtie2-build GCF_000149955.1_ASM14995v2_genomic.fna index/index

for i in 1 3 4 5 7 8 9 11 12 14 34 35 36 37 38;
do
path1=Fusarium/${i}/clean_data/${i}_read1.fq.gz
path2=Fusarium/${i}/clean_data/${i}_read2.fq.gz
echo "Path1 is: "$path1
echo "Path2 is: "$path2
echo "Aligning with bowtie2.."
bowtie2 -x index/index -1 ${path1} -2 ${path2} -S aligned_files/${i}.sam
echo "Changing to bam and sorting..."
samtools view --bam aligned_files/${i}.sam | samtools sort - > bam_files/${i}_sorted.bam
echo "Sorting done, indexing..."
samtools index bam_files/${i}_sorted.bam
done
