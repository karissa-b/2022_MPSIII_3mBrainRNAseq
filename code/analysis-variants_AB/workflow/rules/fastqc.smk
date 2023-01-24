rule fastqc_raw:
    input:
        os.path.join("results", rawData_dir, "fastq", "{SAMPLE}{MERGETAG, .*}{PAIRTAG, .*}" + config["fastq_ext"]),
    output:
        multiext(os.path.join("results", rawData_dir, "FastQC", "{SAMPLE}{MERGETAG, .*}{PAIRTAG, .*}_fastqc"), ".zip", ".html"),
    params:
        outDir = os.path.join("results", rawData_dir, "FastQC"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00",
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule fastqc_trim:
    input:
        os.path.join("results", trim_dir, "fastq", "{SAMPLE}{MERGETAG, .*}{PAIRTAG, .*}" + config["fastq_ext"]),
    output:
        multiext(os.path.join("results", trim_dir, "FastQC", "{SAMPLE}{MERGETAG, .*}{PAIRTAG, .*}_fastqc"), ".zip", ".html"),
    params:
        outDir = os.path.join("results", trim_dir, "FastQC"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 2000,
        time = "00-01:00:00",
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}
        """

rule fastqc_align:
    input:
        os.path.join("results", align_dir, "bam", "{SAMPLE}{MERGETAG, .*}.bam"),
    output:
        multiext(os.path.join("results", align_dir, "FastQC", "{SAMPLE}{MERGETAG, .*}_fastqc"), ".zip", ".html"),
    params:
        outDir = os.path.join("results", align_dir, "FastQC"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00",
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}
        """