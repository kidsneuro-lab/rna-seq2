configfile: "config/config.yaml"

include: "../../workflow/rules/common.smk"
include: "../../workflow/rules/star_index.smk"

rule all:
    input:
        [expand("{genomepath}/{genome}_index", genome=key, genomepath=value) for key, value in get_genomes()]