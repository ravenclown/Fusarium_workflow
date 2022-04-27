configfile:"config.yaml"

rule all:
    input:
        expand("{specimen}_1.fastq.gz", specimen= config["specimen"]),
        expand("{specimen}_2.fastq.gz", specimen= config["specimen"]),
        expand("trimmed_{specimen}_1.fastq.gz", specimen= config["specimen"]),
        expand("trimmed_{specimen}_1.fastq.gz", specimen= config["specimen"]),
        expand("aligned_files/{specimen}_sorted.bam", specimen= config["specimen"]),
        expand("aligned_files/{specimen}_sorted.bam.bai", specimen= config["specimen"])

rule fastqc:
    input:
        read1="{specimen}_1.fastq.gz",
        read2="{specimen}_2.fastq.gz"
    output:
        "fastqc_reports/{specimen}_1.html",
        "fastqc_reports/{specimen}_2.hmtl"
    conda:
        "env.yml"
    shell:
        "mkdir -p fastqc_reports && fastqc {input.read1} {input.read2} -o fastqc_reports"

rule fastP:
    input:
        read1="{specimen}_1.fastq.gz",
        read2="{specimen}_2.fastq.gz"
    output:
        out1="trimmed_reads/trimmed_{specimen}_1.fastq.gz",
        out2="trimemd_reads/trimmed_{specimen}_2.fastq.gz"
    conda:
        "env.yml"
    shell:
        "mkdir -p trimmed_reads && \
        fastp -i {input.read1} -o {output.out1}\
        -I {input.read2} -O {output.out2}"

rule bowtie2_index:
    input:
        "reference.fasta"
    params:
        "bowtie_index/reference_index"
    output:
        output1="bowtie_index/reference_index.1.bt2",
        output2="bowtie_index/reference_index.2.bt2",
        output3="bowtie_index/reference_index.3.bt2",
        output4="bowtie_index/reference_index.4.bt2",
        outputrev1="bowtie_index/reference_index.rev1.bt2",
        outputrev2="bowtie_index/reference_index.rev2.bt2"
    shell:
        "mkdir -p bowtie_index && \
        bowtie2-build {input} {params}"

rule bowtie2:
    input:
        read1="{specimen}_1.fastq.gz",
        read2="{specimen}_2.fastq.gz"
    params:
        "bowtie_index/reference_index"
    output:
        sam="aligned_files/{specimen}.sam"
    conda:
        "sam-bam.env"
    shell:
        "mkdir -p aligned_files && \
        bowtie2 -x {params} -1 {input.read1} -2 {input.read2} -S {output}"

rule samtools:
    input:
        sam="aligned_files/{specimen}.sam"
    output:
        bam="aligned_files/{specimen}.bam",
        sorted_bam="aligned_files/{specimen}_sorted.bam",
        bai="aligned_files/{specimen}_sorted.bam.bai"
    conda:
        "sam-bam.env"
    shell:
        "samtools view --bam {input.sam} > {output.bam} \
        samtools sort {output.bam} > {output.sorted_bam} \
        samtools index {output.sorted_bam}"
