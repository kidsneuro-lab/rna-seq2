import os
from snakemake.shell import shell

rule fastqc_pre_trim:
    input:
        r1=lambda wildcards: os.path.join(config['data']['fastq'], config['samples'][wildcards.sample]['reads']['r1']),
        r2=lambda wildcards: os.path.join(config['data']['fastq'], config['samples'][wildcards.sample]['reads']['r2'])
    output:
        fastqc=directory("results/{sample}/fastqc_pre_trim"),
        done=touch("results/{sample}/fastqc_pre_trim.done")
    log:
        "logs/{sample}/fastqc_pre_trim.log"
    threads:
        workflow.cores
    shell:
        "mkdir -p {output.fastqc} && fastqc {params} -t {threads} "
        "--outdir {output.fastqc} {input.r1} {input.r2} > {log} 2>&1"

rule trim:
    input:
        fastqc_pre_trim_done=rules.fastqc_pre_trim.output['done'],
        r1=lambda wildcards: os.path.join(config['data']['fastq'], config['samples'][wildcards.sample]['reads']['r1']),
        r2=lambda wildcards: os.path.join(config['data']['fastq'], config['samples'][wildcards.sample]['reads']['r2'])
    output:
        r1=os.path.join(config['data']['trimmed'], "{sample}/{sample}.1.paired.fastq.gz"),
        r2=os.path.join(config['data']['trimmed'], "{sample}/{sample}.2.paired.fastq.gz"),
        # reads where trimming entirely removed the mate
        r1_unpaired=os.path.join(config['data']['trimmed'], "{sample}/{sample}.1.unpaired.fastq.gz"),
        r2_unpaired=os.path.join(config['data']['trimmed'], "{sample}/{sample}.2.unpaired.fastq.gz")
    log:
        "logs/{sample}/trimmomatic.log"
    params:
        # list of trimmers (see manual)
        trimmer="ILLUMINACLIP:config/nextera.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
        # trimmer="ILLUMINACLIP:/usr/local/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10",
        # trimmer="ILLUMINACLIP::2:30:10",
        extra=""
        # optional parameters
        # extra="-phred33"
    threads:
        workflow.cores
        # optional specification of memory usage of the JVM that snakemake will respect with global
        # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
        # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
        # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    shell:
        "trimmomatic PE -threads {threads} {params.extra} "
        "{input.r1} {input.r2} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer} > {log} 2>&1"

rule fastqc_post_trim:
    input:
        r1=rules.trim.output['r1'],
        r2=rules.trim.output['r2']
    output:
        done=touch("results/{sample}/fastqc_post_trim.done"),
        fastqc=directory("results/{sample}/fastqc_post_trim")
    log:
        "logs/{sample}/fastqc_post_trim.log"
    threads:
        workflow.cores
    shell:
        "mkdir -p {output.fastqc} && fastqc {params} -t {threads} "
        "--outdir {output.fastqc} {input.r1} {input.r2} > {log} 2>&1"


