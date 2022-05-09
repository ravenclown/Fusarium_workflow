configfile:"config.yaml"
rule all:
    input:
        "reference.fasta",
        expand("fastqc_reports/{sample}_1.html",sample= config["sample"]),
        expand("fastqc_reports/{sample}_2.html",sample= config["sample"]),
        expand("fastqc_report/{sample}_trimmed_1.html",sample= config["sample"]),
        expand("fastqc_report/{sample}_trimmed_2.html",sample= config["sample"]),
        "index/index.1.bt2",
        "index/index.2.bt2",
        "index/index.3.bt2",
        "index/index.4.bt2",
        "index/index.rev.1.bt2",
        "index/index.rev.2.bt2",
        expand("aligned_files/{sample}_sorted.bam", sample= config["sample"]),
        expand("aligned_files/{sample}_sorted.bam.bai", sample= config["sample"]),
        expand("aligned_files/{sample}_gatk.vcf", sample= config["sample"]),
        expand("aligned_files/{sample}_gatk_snp.vcf", sample= config["sample"])


rule fastqc:
    input:
        read1="reads/{sample}_read1.fq.gz",
        read2="reads/{sample}_read2.fq.gz",
        trim1="trimmed_reads/trimmed_{sample}_read1.fq.gz",
        trim2="trimmed_reads/trimmed_{sample}_read2.fq.gz"
    output:
        "fastqc_reports/{sample}_1.html",
        "fastqc_reports/{sample}_2.hmtl",
        "fastqc_reports/{sample}_trimmed_1.html",
        "fastqc_reports/{sample}_trimmed_2.html"
    conda:
        "env.yml"
    shell:
        "mkdir -p fastqc_reports && fastqc {input.read1} {input.read2} {input.trim1} {input.trim2} -o fastqc_reports"

rule fastP:
    input:
        read1="reads/{sample}_read1.fq.gz",
        read2="reads/{sample}_read2.fq.gz"
    output:
        out1=temp("trimmed_reads/trimmed_{sample}_read1.fq.gz"),
        out2=temp("trimemd_reads/trimmed_{sample}_read2.fq.gz")
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
        read1="reads/{sample}_read1.fq.gz",
        read2="reads/{sample}_read2.fq.gz"
    params:
        "index/index"
    output:
        sam=temp("aligned_files/{sample}.sam")
    conda:
        "sam-bam.yml"
    shell:
        "mkdir -p aligned_files && \
        bowtie2 -x {params} -1 {input.read1} -2 {input.read2} -S {output}"

rule samtools:
    input:
        sam="aligned_files/{sample}.sam"
    output:
        sorted_bam="aligned_files/{sample}_sorted.bam",
        bai="aligned_files/{sample}_sorted.bam.bai"
    conda:
        "sam-bam.yml"
    shell:
        "samtools view --bam {input.sam} |  samtools sort - > {output.sorted_bam} && \
        samtools index {output.sorted_bam}"

rule PrepareReference:
    input:
        "reference.fasta"
    output:
        dict="reference.dict",
        fai="reference.fasta.fai"
    shell:
        "./gatk-4.2.6.1/gatk CreateSequenceDictionary -R {input} && \
        samtools faidx {input}"
        
rule MarkDuplicates:
    input:
        "{specimen}_sorted.bam"
    output:
        md=temp("{specimen}_MD.bam"),
        metric="{specimen}_duplicate_metrics.txt"
    conda:
        "gatk.env"
    shell:
        "picard MarkDuplicates -I {input} -M {output.metric} -O {output.md}"

rule SortSam:
    input:
        "{specimen}_MD.bam"
    output:
        bam="{specimen}_clean_sorted.bam"
    shell:
        "picard SortSam -I {input} -O {output.bam} -SO coordinate"

rule HaplotypeCaller:
    input:
        bam="{specimen}_clean_sorted.bam",
        ref="reference.fasta",
        dict="reference.dict",
        fai="reference.fasta.fai"
    output:
        "{specimen}_gatk.vcf"
    conda:
        "gatk.env"
    shell:
        "./gatk-4.2.6.1/gatk HaplotypeCaller -I {input.bam} -R {input.ref} -O {output} -ploidy 1"

rule SelectVariants:
    input:
        vcf="{specimen}_gatk.vcf",
        ref="reference.fasta",
        dict="reference.dict",
        fai="reference.fasta.fai"
    output:
        "{sample}_gatk_snp.vcf"
    conda:
        "gatk.env"
    shell:
        "./gatk-4.2.6.1/gatk SelectVariants \
         -R {input.ref} \
         -V {input.vcf} \
         --select-type-to-include SNP \
         -O {output}"
