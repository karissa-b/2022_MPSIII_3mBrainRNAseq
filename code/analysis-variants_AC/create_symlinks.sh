#!/bin/bash

## Setup snakemake pipeline by creating symlinks to raw data for specified samples
## Rename symlinks as required by pipeline while leaving raw data filenames in tact

declare -a SAMPS=(
    "22-01839_S1"
    "22-01840_S2"
    "22-01841_S3"
    "22-01842_S4"
    "22-01843_S5"
    "22-01844_S6"
    "22-01845_S7"
    "22-01846_S8"
    "22-01847_S9"
    "22-01848_S10"
    "22-01849_S11"
    "22-01850_S12"
    "22-01851_S13"
    "22-01852_S14"
    "22-01853_S15"
    "22-01854_S16"
    "22-01855_S17"
    "22-01856_S18"
    "22-01857_S19"
    "22-01858_S20"
)

RAW_DIR=/hpcfs/users/a1647910/2022_MPSIII_3mBrainRNAseq/raw_data
DEST_DIR=/hpcfs/users/a1647910/2022_MPSIII_3mBrainRNAseq/code/analysis-variants_AC/results/00_rawData/fastq/

for i in "${SAMPS[@]}"
do
    ## == Run 1 ==
    ## Loop through fastq files in run1 raw data
    for file in ${RAW_DIR}/SAGCQA0413_run1/fastq/L0[1-4]/${i}*.fastq.gz
    do
        # Create the symbolic link
        ln -s ${file} ${DEST_DIR}
        # Rename the symlink so that it contains the run number in its filename
        # This will be used for BQSR and when merging files in the snakemake workflow
        SYMLINK=${DEST_DIR}$(basename ${file})
        # Here we use sed grouping, defined using parenthesis, to retain the lane number when renaming
        mv ${SYMLINK} $(echo ${SYMLINK} | sed -E 's/_(L0[1-4])_R/_\1_run1_R/g')
    done
    ## == Run 2 ==
    ## Do the same for run 2
    for file in ${RAW_DIR}/SAGCQA0413_run2/fastq/L0[1-4]/${i}*.fastq.gz
    do
        ln -s ${file} ${DEST_DIR}
        SYMLINK=${DEST_DIR}$(basename ${file})
        mv ${SYMLINK} $(echo ${SYMLINK} | sed -E 's/_(L0[1-4])_R/_\1_run2_R/g')
    done
done