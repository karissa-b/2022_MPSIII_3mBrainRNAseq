#!/bin/bash

## Setup snakemake pipeline by creating symlinks to raw data for specified samples
## Rename symlinks as required by pipeline while leaving raw data filenames in tact

declare -a SAMPS=(
    "22-01859_S21"
    "22-01860_S22"
    "22-01861_S23"
    "22-01862_S24"
    "22-01863_S25"
    "22-01864_S26"
    "22-01865_S27"
    "22-01866_S28"
    "22-01867_S29"
    "22-01868_S30"
    "22-01869_S31"
    "22-01870_S32"
    "22-01871_S33"
    "22-01872_S34"
    "22-01873_S35"
    "22-01874_S36"
    "22-01875_S37"
    "22-01876_S38"
    "22-01877_S39"
    "22-01879_S40"
    "22-01881_S41"
    "22-01882_S42"
    "22-01883_S43"
    "22-01884_S44"
    "22-01885_S45"
)

RAW_DIR=/hpcfs/users/a1647910/2022_MPSIII_3mBrainRNAseq/raw_data
DEST_DIR=/hpcfs/users/a1647910/2022_MPSIII_3mBrainRNAseq/code/analysis-variants_BC/results/00_rawData/fastq/

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