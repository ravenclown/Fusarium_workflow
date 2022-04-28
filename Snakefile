configfile:"config.yaml"

rule all:
    input:
        "reference.fasta",
        expand("fastqc_reports/{specimen}_1.html",specimen= config["specimen"]),
        expand("fastqc_reports/{specimen}_2.html",specimen= config["specimen"]),
        expand("trimmed_{specimen}_read1.fq.gz", specimen= config["specimen"]),
        expand("trimmed_{specimen}_read2.fq.gz", specimen= config["specimen"]),
        "index/index.1.bt2",
        "index/index.2.bt2",
        "index/index.3.bt2",
        "index/index.4.bt2",
        "index/index.rev.1.bt2",
        "index/index.rev.2.bt2",
        expand("aligned_files/{specimen}_sorted.bam", specimen= config["specimen"]),
        expand("aligned_files/{specimen}_sorted.bam.bai", specimen= config["specimen"])

rule fastqc:
    input:
        read1="reads/{specimen}_read1.fq.gz",
        read2="reads/{specimen}_read2.fq.gz"
    output:
        "fastqc_reports/{specimen}_1.html",
        "fastqc_reports/{specimen}_2.hmtl"
    conda:
        "env.yml"
    shell:
        "mkdir -p fastqc_reports && fastqc {input.read1} {input.read2} -o fastqc_reports"

rule fastP:
    input:
        read1="reads/{specimen}_read1.fq.gz",
        read2="reads/{specimen}_read2.fq.gz"
    output:
        out1="trimmed_reads/trimmed_{specimen}_read1.fq.gz",
        out2="trimemd_reads/trimmed_{specimen}_read2.fq.gz"
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
        "index/index"
    output:
        output1="index/index.1.bt2",
        output2="index/index.2.bt2",
        output3="index/index.3.bt2",
        output4="index/index.4.bt2",
        outputrev1="index/index.rev1.bt2",
        outputrev2="index/index.rev2.bt2"
    conda:
        "sam-bam.yml"
    shell:
        "mkdir -p bowtie_index && \
        bowtie2-build {input} {params}"

rule bowtie2:
    input:
        "index/index.1.bt2",
        "index/index.2.bt2",
        "index/index.3.bt2",
        "index/index.4.bt2",
        "index/index.rev.1.bt2",
        "index/index.rev.2.bt2",
        read1="reads/{specimen}_read1.fq.gz",
        read2="reads/{specimen}_read2.fq.gz"

    params:
        "index/index"
    output:
        sam="aligned_files/{specimen}.sam"
    conda:
        "sam-bam.yml"
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
        "sam-bam.yml"
    shell:
        "samtools view --bam {input.sam} |  samtools sort - > {output.sorted_bam} && \
        samtools index {output.sorted_bam} && \
        rm aligned_files/{wildcards.specimen}.sam"
