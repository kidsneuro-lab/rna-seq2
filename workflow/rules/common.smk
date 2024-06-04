import os

def get_sample_configs() -> list:
    sample_variations = []
    for sample in list(config['samples'].keys()):
        for variation in list(config['samples'][sample]['variations'].keys()):
            sample_variations.append(f"{sample}_{variation}")

    return sample_variations


def get_star_config(sample: str, conf: str) -> str:
    STAR_parameters = []
    sample_variation = re.sub(f"^{sample}_", "", conf) # e.g. 'default' from 430_RERE_default

    for param in config['samples'][sample]['variations'][sample_variation]:
        if param in config['config']['star']['alignment'].keys():
            STAR_parameters.extend(config['config']['star']['alignment'][param])
        else:
            STAR_parameters.append(param)    
    
    return ' '.join(STAR_parameters)
    
def get_sample_variations():
    for sample in config['samples'].keys():
        for variation in config['samples'][sample]['variations']:
            yield sample, f"{sample}_{variation}"

def regions_of_interest() -> str:
    region_of_interest = []

    for region in config['region_of_interest']:
        region_of_interest.append(region)

    return ' '.join(region_of_interest)

def get_genomes():
    for genome in config['config']['star']['genomes'].keys():
        yield genome, config['config']['star']['genomes'][genome]['data']

def star_index_path(sample: str) -> str:
    genome=config['samples'][sample]['genome']
    genome_path=config['config']['star']['genomes'][genome]['data']

    return os.path.join(genome_path, f"{genome}_index")

def star_genome(sample: str) -> str:
    genome=config['samples'][sample]['genome']
    fasta_path=os.path.join(config['config']['star']['genomes'][genome]['data'], config['config']['star']['genomes'][genome]['fasta'])

    return fasta_path