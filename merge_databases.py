import os
import sqlite3
import glob

# Define batch job number from running previous BASH script
batch_job_ID = 12479328

# Define directories and corresponding output file names
directories = {
#    "checks": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/Checks",
#    "pressure": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/Pressure",
#    "cs2": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/CS2",
    "eden": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/Eden",
    "nb": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/nB"
}

output_files = {
#    "checks": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/eos_checks_ALL_{batch_job_ID}.sqlite",
#    "pressure": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/eos_pressure_ALL_{batch_job_ID}.sqlite",
#    "cs2": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/eos_cs2_ALL_{batch_job_ID}.sqlite",
    "eden": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/eos_eden_ALL_{batch_job_ID}.sqlite",
    "nb": f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_ID}/Databases/eos_nB_ALL_{batch_job_ID}.sqlite"
}

for dir_name, dir_path in directories.items():
    output_db_path = output_files[dir_name]
    print(f"Creating combined database: {output_db_path}")

    # Create a new SQLite database to combine tables into
    conn_output = sqlite3.connect(output_db_path)
    cur_output = conn_output.cursor()

    merged_databases = 0

    # Loop through each file in the directory
    for db_file in glob.glob(os.path.join(dir_path, "*.sqlite")):
        print(f"Processing file: {db_file}")

        # Attach the input database file
        attach_name = f"attached_db_{merged_databases}"
        cur_output.execute(f"ATTACH DATABASE '{db_file}' AS {attach_name}")

        # Get the list of tables in the input database
        cur_output.execute(f"SELECT name FROM {attach_name}.sqlite_master WHERE type='table';")
        tables = cur_output.fetchall()

        # Loop through each table and copy its content to the output database
        for table in tables:
            table_name = table[0]

            if table_name == "sqlite_sequence":
                continue  # Skip the internal sqlite_sequence table

            # Create the table if it does not exist in the output database
            cur_output.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}';")
            if not cur_output.fetchall():
                cur_output.execute(f"CREATE TABLE {table_name} AS SELECT * FROM {attach_name}.{table_name} WHERE 1=0")

            # Copy data from the attached table to the output table
            cur_output.execute(f"INSERT OR IGNORE INTO {table_name} SELECT * FROM {attach_name}.{table_name}")

        # Commit any pending transactions before detaching
        conn_output.commit()

        # Detach the input database
        cur_output.execute(f"DETACH DATABASE {attach_name}")

        merged_databases += 1

    # Commit and close the connection to the output database
    conn_output.commit()
    conn_output.close()

    print(f"All tables from {dir_name} have been combined successfully into {output_db_path}.")
    print(f"Merged {merged_databases} databases from {dir_name}.")
