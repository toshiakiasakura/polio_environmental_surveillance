#!/bin/bash
#SBATCH --job-name=polio_stochastic1 # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lsh2200567@lshtm.ac.uk # Where to send mail
#SBATCH --ntasks=1                    # Run on a single core
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=120:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test3_%j.log # Standard output and error log
pwd; hostname; date;
echo "time started";
date;
module load singularity;
singularity exec --bind /home/lsh2200567/polio_environmental_surveillance:/workdir \
    --env JULIA_DEPOT_PATH="/workdir/.julia::/opt/julia" \
    /home/lsh2200567/polio_environmental_surveillance/image-2024-02-23.simg \
    /bin/sh compute3.sh;

echo "time ended";
date;
module load singularity
