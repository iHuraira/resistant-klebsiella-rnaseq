configfile: "config/config.yaml"

with open(config["sra_path"], "r") as sra_text_file:
    samples = [line.strip() for line in sra_text_file if line.strip()]

include: "rules/sra_tools.smk"

rule all:
    input:
        expand("results/sra/{sample}/{sample}.sra", sample=samples)
