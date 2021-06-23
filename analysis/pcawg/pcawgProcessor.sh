#!/bin/bash



wd="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj"
job_dir="${wd}/slurm_out"
job_name="pcawg-lung-subsetting"
dt=$(date +%y%m%d_%T | sed -e "s/:/-/g")



job_script=${job_dir}/jobs/${dt}.${job_name}.jobscript.sh

if [[ ! -p ${job_script} ]]; then
        touch ${job_script};
fi


echo "#!/bin/bash

#SBATCH --job-name=${job_name}
#SBATCH --output=${job_dir}/outs/${dt}.${job_name}.output.txt
#SBATCH --error=${job_dir}/errs/${dt}.${job_name}.error.txt
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=a.movasati@uu.nl

guixr load-profile /hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/.guix-profile-berner-proj --<<EOF

Rscript ${wd}/analysis/pcawg/explorePcawgLungCohort.R

EOF" > ${job_script}


sbatch ${job_script}


