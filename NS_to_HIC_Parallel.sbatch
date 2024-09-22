#!/bin/bash
#SBATCH --time=05:00:01
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-49  # Running 50 tasks
#SBATCH --output=projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_%A/Outputs/output_%A_%a.txt
#SBATCH --error=projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_%A/Errors/error_%A_%a.txt
#SBATCH --account=qgp
#SBATCH --partition=qgp
#SBATCH --mail-user=dolse0058@ung.edu
#SBATCH --mail-type=BEGIN,END

# Number of EOS files to process
N=500

# Load necessary modules
module load python/3.9.16
module load anaconda/2022-May/3

# Create directories for output and error logs only once
if [ $SLURM_ARRAY_TASK_ID -eq 0 ]; then
    batch_dir=projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_${SLURM_ARRAY_JOB_ID}
    mkdir -p ${batch_dir}/Outputs
    mkdir -p ${batch_dir}/Errors
    mkdir -p ${batch_dir}/Databases/Checks
    mkdir -p ${batch_dir}/Databases/Pressure
    mkdir -p ${batch_dir}/Databases/CS2
    mkdir -p ${batch_dir}/Databases/Eden
    mkdir -p ${batch_dir}/Databases/nB

    # Create a batch job info file with the creation details
    timestamp=$(date +"Batch Job ${SLURM_ARRAY_JOB_ID} was submitted at %d %B %Y at %I:%M:%S %p")
    echo $timestamp > ${batch_dir}/Batch_Job_${SLURM_ARRAY_JOB_ID}_Info.txt
fi

# Ensure the directories exist for all jobs
batch_dir=projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_${SLURM_ARRAY_JOB_ID}
while [ ! -d ${batch_dir}/Outputs ]; do
    sleep 1
done
while [ ! -d ${batch_dir}/Errors ]; do
    sleep 1
done
while [ ! -d ${batch_dir}/Databases/Checks ]; do
    sleep 1
done
while [ ! -d ${batch_dir}/Databases/Pressure ]; do
    sleep 1
done
while [ ! -d ${batch_dir}/Databases/CS2 ]; do
    sleep 1
done
while [ ! -d ${batch_dir}/Databases/Eden ]; do
    sleep 1
done
while [ ! -d ${batch_dir}/Databases/nB ]; do
    sleep 1
done

# Calculate the start and end index for N EOS files
start_index=$(( SLURM_ARRAY_TASK_ID * N ))
end_index=$(( start_index + N - 1 ))

# Run Python code
output_path=${batch_dir}/Outputs/output_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
error_path=${batch_dir}/Errors/error_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

echo "Starting job $SLURM_ARRAY_TASK_ID at $(date)" >> $output_path
python projects/jnorhos/dolsen/NS_to_HIC_Bash/convert_ns_to_hic.py $start_index $end_index $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID ${batch_dir}/Databases/Checks ${batch_dir}/Databases/Pressure ${batch_dir}/Databases/CS2 ${batch_dir}/Databases/Eden ${batch_dir}/Databases/nB >> $output_path 2>> $error_path

echo "Finished job $SLURM_ARRAY_TASK_ID at $(date)" >> $output_path
