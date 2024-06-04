# RNASeq alignment pipeline
A pipeline to run multiple alignments for RNASeq using the STAR alignment. The intent here is to test multiple STAR alignment configurations including insertion of custom splice junctions.

# Pre-requities
## Software
- bedtools
- fastqc
- htslib
- python 3.10
- samtools
- star
- tabix
- trimmomatic

**Note:** The workflow has been tested with these specific versions. Updated versions may work
```
module load bedtools/2.29.2
module load fastqc/0.11.8
module load htslib/1.10.2
module load python/3.9.15
module load samtools/1.9
module load star/2.7.8a
module load tabix/0.2.6
module load trimmomatic/0.39
```

## Genome reference
- hg38 fasta
- hg38 annotation (e.g. Gencode)

## Runtime environment
- HPC using slurm scheduler (e.g. qsub)
- It is possible to extract the commands from the `*.qsub` files and run them using a different workflow engine

# Set up
## Create a virtual environment
```shell
python -m venv venv_rnaseq2
. venv_rnaseq2/bin/activate
pip3 install -U setuptools wheel pip
```

## Install requirements
```shell
pip3 install -r requirements.txt
```
# Running the workflow
## Configure the workflow
- Update `config/config.yaml`
- Refer to the configuration guide below

## Generating STAR index (This is a one off task)
`qsub run_star_index.pbs`

## Generating STAR index
`qsub run_star_index.pbs`

## Trim the FASTQ
- Check adapter sequence
- Update `config/nextera.fasta`
- Run `qsub run_pe_trim.pbs`

## Run alignment
qsub `run_pe_align.pbs`

# Configuration guide
## `config/config.yaml`
```yaml
star:
    genome: /<location>/Homo_sapiens_assembly38.fasta
    index: /<location>/star_index_gencode38
    gtf: /<location>/gencode.v38.annotation.gtf
    sjdbOverhang: 150

data:
    fastq: /<location>/raw
    trimmed: /<location>/trimmed

region_of_interest: # This is for creation of the subset BAM file only containing your region of interest
    - "chr1:8,557,893-8,561,603" # chr1:8,350,404-8,819,640

samples:
     "SAMPLE":
         "reads":
             "r1": "SAMPLE_R1.fastq.gz"
             "r2": "SAMPLE_R2.fastq.gz"
         "variations":
             "default":
                 - "default"

config:
    star:
        default:
            - "--outSAMtype BAM Unsorted"
            - "--outSAMattributes NH HI AS NM MD"
            - "--outFilterMultimapNmax 20"
            - "--alignSJoverhangMin 6"
            - "--alignSJDBoverhangMin 6"
            - "--outFilterMismatchNmax 10"      
            - "--outFilterMismatchNoverLmax 0.3"
            - "--alignIntronMin 21"
            - "--alignIntronMax 1000000"
            - "--alignMatesGapMax 0"
            - "--outFilterType BySJout"
            - "--outFilterScoreMinOverLread 0.66"
            - "--outFilterMatchNminOverLread 0.66"
            - "--outSAMstrandField intronMotif"
            - "--outFilterIntronMotifs None"
            - "--alignSoftClipAtReferenceEnds Yes"
            - "--quantMode GeneCounts"
            - "--outSAMunmapped Within"
            - "--chimSegmentMin 12"
            - "--chimJunctionOverhangMin 12"
            - "--chimOutType WithinBAM"
            - "--chimMainSegmentMultNmax 10"
            - "--outSJfilterCountUniqueMin 3 1 1 1"
```
