import sqlite3
import pandas as pd

def read_sql_table(database_path):
    try:
        # Connect to the SQLite database
        conn = sqlite3.connect(database_path)

        # Get the table names in the database
        query = "SELECT name FROM sqlite_master WHERE type='table';"
        table_names = pd.read_sql_query(query, conn)

        if table_names.empty:
            print("No tables found in the database.")
            return

        # Display the names of the tables
        print("Tables in the database:")
        print(table_names)

        # Loop through the tables and display their contents
        for table in table_names['name']:
            print(f"\nContents of table '{table}':")
            df = pd.read_sql_query(f"SELECT * FROM {table};", conn)
            print(df)

        # Close the connection
        conn.close()

    except sqlite3.Error as e:
        print(f"Error reading SQLite database: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    batch_job_id = 12473748
    database_path = f'projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_id}/Databases/Checks/eos_checks_{batch_job_id}_0.sqlite'
    read_sql_table(database_path)
