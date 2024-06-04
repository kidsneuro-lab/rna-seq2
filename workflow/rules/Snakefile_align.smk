configfile: "config/config.yaml"

include: "../../workflow/rules/common.smk"
include: "../../workflow/rules/trim.smk"
include: "../../workflow/rules/align.smk"
include: "../../workflow/rules/sort_index.smk"
include: "../../workflow/rules/subset_sort_index.smk"

rule all:
    input:
        [expand("results/{sample}/{conf}/{conf}.sorted.cram.crai", sample=key, conf=value) for key, value in get_sample_variations()],
        [expand("results/{sample}/{conf}/{conf}.sorted.subset.bam.bai", sample=key, conf=value) for key, value in get_sample_variations()]
