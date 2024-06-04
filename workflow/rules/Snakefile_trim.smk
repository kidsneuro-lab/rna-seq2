configfile: "config/config.yaml"

include: "../../workflow/rules/common.smk"
include: "../../workflow/rules/trim.smk"

rule all:
    input:
        [expand("results/{sample}/fastqc_post_trim", sample=key, conf=value) for key, value in get_sample_variations()]
