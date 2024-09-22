#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --output=projects/jnorhos/dolsen/NS_to_HIC_Bash/Merge_Job_%j/output_merge_%j.txt
#SBATCH --error=projects/jnorhos/dolsen/NS_to_HIC_Bash/Merge_Job_%j/error_merge_%j.txt
#SBATCH --account=qgp
#SBATCH --partition=qgp
#SBATCH --mail-user=dolse0058@ung.edu
#SBATCH --mail-type=BEGIN,END

# Load necessary modules
module load python/3.9.16
module load anaconda/2022-May/3

# Create directories for output and error logs
merge_dir=projects/jnorhos/dolsen/NS_to_HIC_Bash/Merge_Job_${SLURM_JOB_ID}
mkdir -p ${merge_dir}/Outputs
mkdir -p ${merge_dir}/Errors
mkdir -p projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_${SLURM_ARRAY_JOB_ID}/Merged_DB

# Run Python code to merge databases
output_path=${merge_dir}/Outputs/output_merge_${SLURM_JOB_ID}.txt
error_path=${merge_dir}/Errors/error_merge_${SLURM_JOB_ID}.txt

echo "Starting merge process at $(date)" >> $output_path
python projects/jnorhos/dolsen/NS_to_HIC_Bash/merge_databases.py projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_${SLURM_ARRAY_JOB_ID}/Databases/ projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_${SLURM_ARRAY_JOB_ID}/Merged_DB >> $output_path 2>> $error_path
echo "Finished merge process at $(date)" >> $output_path
