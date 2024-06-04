rule sort_align_bam:
    input:
        rules.align.output['bam']
    output:
        bam=temp("results/{sample}/{conf}/Aligned.out.sorted.bam")
    log:
        "logs/{sample}/{conf}/sort_align_bam.log"
    params:
        extra = lambda wildcards: "--reference " + star_genome(wildcards.sample),
    threads:
        workflow.cores
    wrapper:
        "file:resources/snakemake-wrappers/bio/samtools/sort"

rule index_align_bam:
    input:
        rules.sort_align_bam.output['bam']
    output:
        bai=temp("results/{sample}/{conf}/Aligned.out.sorted.bam.bai")
    log:
        "logs/{sample}/{conf}/index_align_bam.log"
    params:
        "-@ {workflow.cores}"
    wrapper:
        "file:resources/snakemake-wrappers/bio/samtools/index"

rule subset_bam:
    input:
        bam=rules.sort_align_bam.output['bam'],
        bai=rules.index_align_bam.output['bai']
    output:
        bam=temp("results/{sample}/{conf}/{conf}.subset.bam")
    log:
        "logs/{sample}/{conf}/subset_bam.log"
    params:
        extra = lambda wildcards: "-O BAM --reference " + star_genome(wildcards.sample),
	regions = regions_of_interest(),
    threads:
        workflow.cores
    shell:
        "samtools view -@ {threads} -b -o {output.bam} {params.extra} "
	"{input.bam} {params.regions} > {log} 2>&1"

rule sort_bam:
    input:
        rules.subset_bam.output['bam']
    output:
        "results/{sample}/{conf}/{conf}.sorted.subset.bam"
    log:
        "logs/{sample}/{conf}/sort_bam.log"
    params:
        extra = lambda wildcards: "--reference " + star_genome(wildcards.sample),
    threads:
        workflow.cores
    wrapper:
        "file:resources/snakemake-wrappers/bio/samtools/sort"

rule index_bam:
    input:
        rules.sort_bam.output
    output:
        "results/{sample}/{conf}/{conf}.sorted.subset.bam.bai"
    log:
        "logs/{sample}/{conf}/index_bam.log"
    params:
        "-@ {workflow.cores}"
    wrapper:
        "file:resources/snakemake-wrappers/bio/samtools/index"
