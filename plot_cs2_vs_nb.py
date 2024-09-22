import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
import numpy as np

# Set the agg.path.chunksize to avoid the OverflowError
plt.rcParams['agg.path.chunksize'] = 10000

# Define the batch job number
batch_job_number = 12479328

# Read the data for the first N rows of EOS n
N = 100
n = 22

# Connect to the databases
db_dir = f"/projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_number}/Databases/"
cs2_db_path = os.path.join(db_dir, f'eos_cs2_ALL_{batch_job_number}.sqlite')
nb_db_path = os.path.join(db_dir, f'eos_nB_ALL_{batch_job_number}.sqlite')

conn_cs2 = sqlite3.connect(cs2_db_path)
conn_nb = sqlite3.connect(nb_db_path)

# Initialize lists to store the data
all_cs2_data = []
all_nb_data = []

for row_num in range(N):
    query_cs2 = f"SELECT cs2 FROM eos_cs2 WHERE eos_number = {n} AND row_number = {row_num}"
    query_nb = f"SELECT nb FROM eos_nb WHERE eos_number = {n} AND row_number = {row_num}"

    df_cs2 = pd.read_sql(query_cs2, conn_cs2).dropna().reset_index(drop=True)
    df_nb = pd.read_sql(query_nb, conn_nb).dropna().reset_index(drop=True)

    # Filter out rows where nb == 0
    df_cs2 = df_cs2[df_nb['nb'] != 0]
    df_nb = df_nb[df_nb['nb'] != 0]

    # Ensure the lengths of the arrays are equal
    min_length = min(len(df_cs2), len(df_nb))
    df_cs2 = df_cs2.iloc[:min_length]
    df_nb = df_nb.iloc[:min_length]

    all_cs2_data.append(df_cs2)
    all_nb_data.append(df_nb)

# Close the connections
conn_cs2.close()
conn_nb.close()

# Load the data from the .dat file
dat_file_path = f"/projects/jnorhos/emilyad/Binary_Love/complete_eos/EOS_{n}.dat"
dat_data = np.loadtxt(dat_file_path, skiprows=1)
dat_nb = dat_data[:, 1]  # Second column for nB
dat_cs2 = dat_data[:, 3]  # Fourth column for cs2

# Filter out rows where nb == 0
valid_idx = dat_nb != 0
dat_nb = dat_nb[valid_idx]
dat_cs2 = dat_cs2[valid_idx]

# Plot the data
plt.figure(figsize=(10, 6))
for cs2_data, nb_data in zip(all_cs2_data, all_nb_data):
    plt.plot(nb_data['nb'], cs2_data['cs2'], '-', color='gray', alpha=0.5)  # SQL data
plt.plot(dat_nb, dat_cs2, 'k-', linewidth=4, label='NS EOS')  # .dat file data
plt.xlabel(r'$n_B \slash n_{sat}$')
plt.ylabel(r'$c_s^2$')
plt.title(f'NS EOS {n} with {N} HIC EOS')

# Save the plot
output_path = f"/projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_number}/Databases/cs2_vs_nb_EOS_{n}_N={N}.png"
plt.savefig(output_path)

# Show the plot
plt.show()
