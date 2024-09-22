These codes are designed to take in 100,000 NS EOS and convert them to HIC. The procedure is this:



  1) Initialize the NS_to_HIC_Parallel.sbatch by changing the output directories in the SLURM headings and in the body. For 100,000 NS EOS, create 200 tasks where each task processes 500 NS EOS (200 x 500 = 100,000).

  2) NS_to_HIC_Parallel.sbatch will create a directory for the Batch job you submit (Batch_Job_#). In this directory will three more directories: Databases, Outputs, and Errors. 

  3) NS_to_HIC_Parallel.sbatch will fill the Databases directory with four EOS directories (CS2, Eden, nB, Pressure) and one directory for the symmetry energy expansion coefficients and stability/causality checks.

  4) For 200 tasks, there will be 200 sqlite files in each directory. Each sqlite file will have information from 500 NS EOS converted to HIC EOS using convert_ns_to_hic.py.

  5) merge_sql_databases.sbatch will merge 200 sqlite files into one sqlite file ending with "ALL.sqlite."

  6) Once merged, you can run plot_cs2_vs_nb. Here, you specify the Batch Job ID, the NS EOS you are plotting (n), and how many HIC EOS you are plotting with it (N). You can also get cornerplots with cornerplots.py

  7) Because SQL cannot be opened easily, read_sql.py allows you to enter a sqlite table to ensure data is stored in the database correctly.



This is the flow of things for this project. Of course, changing the directories and keeping up with the Batch Job ID are the most important things to keep track of. I used Cyberduck (https://cyberduck.io/) to easily navigate the directories, it was extremely helpful.

