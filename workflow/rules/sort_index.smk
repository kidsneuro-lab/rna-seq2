rule sort_cram:
    input:
        rules.align.output['bam']
    output:
        "results/{sample}/{conf}/{conf}.sorted.cram"
    params:
        extra = lambda wildcards: "-O CRAM --reference " + star_genome(wildcards.sample),
        tmp_dir = directory("/tmp/")
    threads:
        workflow.cores
    wrapper:
        "file:resources/snakemake-wrappers/bio/samtools/sort"

rule index_cram:
    input:
        rules.sort_cram.output
    output:
        "results/{sample}/{conf}/{conf}.sorted.cram.crai"
    params:
        "-@ {workflow.cores}"
    wrapper:
        "file:resources/snakemake-wrappers/bio/samtools/index"
