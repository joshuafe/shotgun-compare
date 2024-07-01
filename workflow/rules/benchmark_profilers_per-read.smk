import os
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


envvars:
    "TMPDIR",
    "SNAKEMAKE_PROFILE",


SHARDS = make_shard_names(config["nshards"])


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)
    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,


kraken_outputs = f"kraken2/{config['sample']}_kraken2.report"
kraken_unclassified_outputs = expand(
    "kraken2/{sample}_kraken2_unclassified_{readdir}.fastq.gz",
    sample=config["sample"],
    readdir=[1, 2],
)
phanta_db_names = ["phanta_unmasked_db", "phanta_masked_db"]
phanta_outputs = expand(
    f'phanta_{config["sample"]}/{{db}}/results/final_merged_outputs/counts.txt',
    db=phanta_db_names,
)

all_outputs = [kraken_outputs, kraken_unclassified_outputs]
if not config["skip_phanta"]:
    all_outputs.extend(phanta_outputs)
    all_outputs.extend(phanta_extra_outputs)


rule all:
    input:
        all_outputs,


#
# Utils Module
if len(config["R1"]) == 1:
    input_R1 = config["R1"]
    input_R2 = config["R2"]
else:
    input_R1 = f"concatenated/{config['sample']}_R1.fastq.gz"
    input_R2 = f"concatenated/{config['sample']}_R2.fastq.gz"


module utils:
    snakefile:
        "utils.smk"
    config:
        config
    skip_validation:
        True


use rule concat_R1_R2 from utils as utils_concat_R1_R2 with:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        R1=temp("concatenated/{sample}_R1.fastq.gz"),
        R2=temp("concatenated/{sample}_R2.fastq.gz"),
    log:
        e="logs/concat_r1_r2_{sample}.e",


rule kraken_standard_run:
    """ Profile the microbiome with Kraken2.

    We run this on the raw reads rather than the host-depleted trimmed reads because
    - its nice to have validation of the host detection percentages
    - we want uniform read lengths for bracken
    - its what the authors seem to recommend, although
      Nick has asked for confirmation: https://github.com/DerrickWood/kraken2/issues/646

    We set the confidence threshold after reading
    https://www.biorxiv.org/content/10.1101/2022.04.27.489753v1.full

    We wanted to
      - minimize the L1 distance to the known community composition
      - maximize the alpha diversity similarity to the actual diversity

    We don't care as much about reads assigned so long as the assignments we
    do get are high quality and representative
    """
    input:
        R1=input_R1,
        R2=input_R2,
        db=config["kraken2_db"],
    output:
        out="kraken2/{sample}_kraken2.out",
        unclass_R1=temp("kraken2/{sample}_kraken2_unclassified_1.fastq"),
        unclass_R2=temp("kraken2/{sample}_kraken2_unclassified_2.fastq"),
        report="kraken2/{sample}_kraken2.report",
    params:
        inpstr=lambda wc, input: (
            f"--paired {input.R1} {input.R2}" if hasattr(input, "R2") else input.R1
        ),
        unclass_template=lambda wildcards, output: output["unclass_R1"].replace(
            "_1.fastq", "#.fastq"
        ),
    container:
        config["docker_kraken"]
    conda:
        "../envs/kraken2.yaml"
    log:
        e="logs/kraken_{sample}.log",
    resources:
        # use 8GB mem if using a minikraken db to make testing easier
        mem_mb=lambda wildcards, attempt: (
            8 * 1024 if "mini" in config["kraken2_db"] else (64 * 1024 * attempt)
        ),
    threads: 16
    shell:
        """
        kraken2 \
            --threads {threads} \
            --use-names \
            --confidence 0.2 \
            --unclassified-out {params.unclass_template} \
            --db {input.db} \
            --report {output.report} \
            {params.inpstr} \
            > {output.out} 2> {log.e}

        # some (mock) datasets are perfect but we still need these files
        if [ ! -f "{output.unclass_R1}" ]
        then
            touch {output.unclass_R1}
            touch {output.unclass_R2}
        fi
        """



rule kraken_merge_shards:
    """ This is only needed when running the whole workflow's snakefile;
    otherwise we do not assume sample is split. with the --only-combined
    option, krakentools will merge individual reports into a single report
    """
    input:
        reports=[
            f"kraken2/{config['sample']}_shard{shard}_kraken2.report"
            for shard in SHARDS
        ],
    output:
        out="kraken2/{sample}_kraken2_merged.report",
    conda:
        "../envs/kraken2.yaml"
    resources:
        mem_mb=2000,
    container:
        config["docker_kraken"]
    threads: 1
    log:
        e="logs/kraken_merge_{sample}.e",
        o="logs/kraken_merge_{sample}.o",
    shell:
        """
        combine_kreports.py \
            --only-combined \
            --no-headers \
            -r {input.reports} \
            -o {output.out} \
            > {log.o} 2> {log.e}
        """


rule make_phanta_manifest:
    """ {params.thisdir} is needed to list the absolute path when using concatenated libraries.
    PHANTA needs the full path, and Input_R1 is relative the relative paths if we are using the concatenated lanes
    """
    input:
        R1=input_R1,
        R2=input_R2,
    output:
        manifest=temp("phanta_inputs.tsv"),
    params:
        thisdir=os.getcwd() if isinstance(input_R1, str) else "",
        sample=config["sample"],
    shell:
        """echo -e "{params.sample}\t{params.thisdir}/{input.R1}\t{params.thisdir}/{input.R2}" > {output.manifest}"""


def get_singularity_args(wildcards):
    with open(os.path.join(os.environ["SNAKEMAKE_PROFILE"], "config.yaml"), "r") as f:
        return yaml.safe_load(f)["singularity-args"]


rule phanta:
    # we need the results dir because if we try to initiate two snakemake workflows from the same --directory we get LockErrors
    input:
        R1=input_R1,
        R2=input_R2,
        manifest=rules.make_phanta_manifest.output.manifest,
        readlen_file=parse_read_lengths,
    output:
        counts="phanta_{sample}/{db}/results/final_merged_outputs/counts.txt",
        bracken_failed="phanta_{sample}/{db}/results/classification/samples_that_failed_bracken.txt",
        bracken_report="phanta_{sample}/{db}/results/classification/{sample}.krak.report_bracken_species.filtered",
        bracken_species_abundance="phanta_{sample}/{db}/results/classification/{sample}.krak.report.filtered.bracken.scaled",
        relative_read_abundance="phanta_{sample}/{db}/results/final_merged_outputs/relative_read_abundance.txt",
        relative_taxonomic_abundance="phanta_{sample}/{db}/results/final_merged_outputs/relative_taxonomic_abundance.txt",
        total_reads="phanta_{sample}/{db}/results/final_merged_outputs/total_reads.tsv",
        single_end=temp(
            directory("phanta_{sample}/{db}/results/classification/single_end")
        ),
    container:
        config["docker_phanta"]
    params:
        dbpath=lambda wc: config[wc.db],
        sing_args=lambda wc: get_singularity_args(wc),
        cov_thresh_viral=config["cov_thresh_viral"],
        cov_thresh_arc=config["cov_thresh_arc"],
        cov_thresh_bacterial=config["cov_thresh_bacterial"],
        cov_thresh_euk=config["cov_thresh_euk"],
        minimizer_thresh_viral=config["minimizer_thresh_viral"],
        minimizer_thresh_arc=config["minimizer_thresh_arc"],
        minimizer_thresh_bacterial=config["minimizer_thresh_bacterial"],
        minimizer_thresh_euk=config["minimizer_thresh_euk"],
        single_end_krak=config["single_end_krak"],
    resources:
        mem_mb=lambda wildcards, attempt: 58 * 1024 * attempt,
    threads: 16
    # we have to do the dummy profile to keep the job from inheriting SNAKEMAKE_PROFILE;
    #   we want this job to be submitted locally.  Otherwise, we have this unpleasant
    # situation where we need the docker container for its snakefile, but also
    # need the workflow itself to execute within the same container.  Cue mount point conflicts,
    # the inability to test, etc.
    # The workaround is running it without a --cluster directive so it runs 'locally',
    # which is actually on the node this rule got submitted to. While there is a slight
    # inefficiency, phanta doesn't benefit enough from distributing its workflow rules to
    # warrant going the alternative route of adding it to this repo as a submodule,
    # importing the rules, etc.
    shell:
        """
        READLEN=$(basename {input.readlen_file} | sed 's|len||' | sed  's|approx||')
        echo "cores: {threads}" > config.yaml && \
        export SNAKEMAKE_PROFILE="" && \
        snakemake  --profile $PWD --notemp \
        --singularity-args "{params.sing_args}" \
        --snakefile /home/mambauser/phanta/Snakefile \
        --configfile /home/mambauser/phanta/config.yaml \
        --directory phanta_{wildcards.sample}/{wildcards.db}/ \
        --config \
        outdir="results/" \
        read_length=$READLEN \
        cov_thresh_viral={params.cov_thresh_viral} \
        cov_thresh_bacterial={params.cov_thresh_bacterial} \
        cov_thresh_euk={params.cov_thresh_euk} \
        cov_thresh_arc={params.cov_thresh_arc} \
        minimizer_thresh_viral={params.minimizer_thresh_viral} \
        minimizer_thresh_bacterial={params.minimizer_thresh_bacterial} \
        minimizer_thresh_euk={params.minimizer_thresh_euk} \
        minimizer_thresh_arc={params.minimizer_thresh_arc} \
        single_end_krak={params.single_end_krak} \
        class_mem_mb={resources.mem_mb} \
        database={params.dbpath} \
        pipeline_directory=/home/mambauser/phanta/ \
        sample_file=$PWD/{input.manifest}
        """
