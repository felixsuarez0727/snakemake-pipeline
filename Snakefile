rule bwa_map:
    input:
        "ref/hg19_chr8.fa",
        "fastq/father_R1.fq.gz",
        "fastq/father_R2.fq.gz"
    output:
        "mapped_reads/father.bam"
    shell:
        "bwa mem {input} | samtools view -Shb -o {output} -"
