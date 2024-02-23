#!/bin/bash
#SBATCH --job-name=polio_stochastic1 # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lsh2200567@lshtm.ac.uk # Where to send mail
#SBATCH --ntasks=1                    # Run on a single core
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test1.log # Standard output and error log
pwd; hostname; date;
echo "time started";
date;
module load singularity;
singularity exec --bind /home/lsh2200567/polio_environmental_surveillance:/workdir \
    --env JULIA_DEPOT_PATH="/workdir/.julia::/opt/julia" \
    polio_stochastic-2023-12-14-bf6ead4075b7.simg \
    /bin/sh cd src && julia 4j_develop_rel_waste.jl;

echo "time ended";
date;
module load singularity
