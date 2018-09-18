#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name= index_1
#SBATCH --output=/home/hklein3/outputs/index_1
#SBATCH --error=/home/hklein3/errors/index_1
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=helenasfklein@gmail.com

module load python3
module load racs-eb
module load matplotlib/2.1.1-intel-2017b-Python-3.6.3

python3 mean_per_base.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
