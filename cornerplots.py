import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import sqlite3

batch_job_id = 12292109

sqlite_db_path = f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_id}/Databases/eos_checks_ALL_{batch_job_id}.sqlite"

conn = sqlite3.connect(sqlite_db_path)
df = pd.read_sql_query("SELECT * FROM eos_checks", conn)
conn.close()

# Drop all columns except eos_number, esym, lsym, ksym, and jsym
df = df[['eos_number', 'esym', 'lsym', 'ksym', 'jsym']]

# EOS number
eos_num = 11

df = df[df['eos_number'] == eos_num]
csv_file_path = f"projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_id}/Databases/eos_cs2_{eos_num}_{batch_job_id}.csv"
df = df[['esym', 'lsym', 'ksym', 'jsym']]
df.to_csv(csv_file_path, index=False, header=False)
data_full = pd.read_csv(csv_file_path, header=None)
data_full.columns = ['Esym', 'Lsym', 'Ksym', 'Jsym']

# Create an instance of the PairGrid class.
grid = sns.PairGrid(data = data_full,height = 4,corner=True,diag_sharey = False)

grid.tight_layout(pad=3,w_pad=8)
grid = grid.map_upper(plt.scatter, color = 'paleturquoise')

col_list = ['Esym', 'Lsym', 'Ksym', 'Jsym']
cols = iter(col_list)

bins = {'Esym' : 5, 'Lsym' : 4,
        'Ksym' : 6, 'Jsym' : 6}

#bins = {'Esym' : 5, 'Lsym' : 4,
#        'Ksym' : 6, 'Jsym' : 6}

# Map a histogram to the diagonal
def myhist(*args, **kwargs):
    b = bins[next(cols)]
    n,nbins,patches=plt.hist(*args, bins=b, density=True, stacked=True, histtype = 'bar', edgecolor = 'k', **kwargs)
  #  sns.histplot(*args,bins=b,stat='count', **kwargs)
   # print(patches)
    cm = plt.cm.get_cmap('BuPu')
    n /= max(n)
    n = n-0.3

   # print(n)
    for c, p in zip(n, patches):
        plt.setp(p, 'facecolor', cm(0.5))
       # print(p.get_height()/sum(n))
     #   p.set_height(p.get_height()/sum(n))
        #plt.bar_label(p)


#grid = grid.map_diag(plt.hist,bins = 15,edgecolor = 'k')
grid = grid.map_diag(myhist)
#grid = grid.map_diag(sns.histplot)


#def mykdeplot(*args, **kwargs):
    #sns.kdeplot(shade=True, cbar=True)
# Map a density plot to the lower triangle
grid = grid.map_lower(sns.kdeplot,cmap = 'BuPu', levels=9, thresh=0.01,bw_adjust=2.7)

#[0][0]
#Turn on axis for corner figure
#grid.axes[0][0].set_ylabel(r'$E_{sym}[MeV]$',fontsize = 30);
grid.axes[0][0].tick_params(axis='both', which='major', labelsize=23)
#grid.axes[0][0].get_yaxis().set_visible(True)
grid.axes[0][0].get_xaxis().set_visible(True)
#grid.axes[0][0].tick_params(axis='y', which='major', labelsize=23)

#grid.axes[0][0].spines[['top', 'right']].set_visible(True)

#[1][0]
grid.axes[1][0].set_ylabel(r'$L_{sym}[MeV]$',fontsize = 30);
grid.axes[1][0].tick_params(axis='both', which='major', labelsize=23)
grid.axes[1][0].get_xaxis().set_visible(True)
#grid.axes[1][0].spines[:].set_visible(True)

#[1][1]
#grid.axes[1][1].spines.set_visible(True)
grid.axes[1][1].get_yaxis().set_visible(True)

#[2][0]
grid.axes[2][0].set_ylabel(r'$K_{sym}[MeV]$',fontsize = 30);
grid.axes[2][0].tick_params(axis='both', which='major', labelsize=23)
#grid.axes[2][0].spines[:].set_visible(True)
#[2][1]
#grid.axes[2][1].spines[:].set_visible(True)

#[2][2]
#grid.axes[2][2].spines[:].set_visible(True)
grid.axes[2][2].get_yaxis().set_visible(True)
#[3][0]
grid.axes[3][0].set_ylabel(r'$J_{sym}[MeV]$',fontsize = 30);
grid.axes[3][0].set_xlabel(r'$E_{sym}[MeV]$',fontsize = 30);
grid.axes[3][0].tick_params(axis='both', which='major', labelsize=23)
#grid.axes[3][0].set_xlim(25, 40)
#ind=[25,30,35,40]
#grid.axes[3][0].set_xticks(ind)
#grid.axes[3][0].spines[:].set_visible(True)

#[3][1]
grid.axes[3][1].set_xlabel(r'$L_{sym}[MeV]$',fontsize = 30);
grid.axes[3][1].tick_params(axis='both', which='major', labelsize=23)
grid.axes[3][1].set_xlim(-20, 140)
ind=[-20,60,140]
grid.axes[3][1].set_xticks(ind)
#grid.axes[3][1].spines[:].set_visible(True)
#ind=[25,75,125]
#grid.axes[3][1].set_xticks(ind)
#label=['25','75', '125']
#[3][2]
grid.axes[3][2].set_xlabel(r'$K_{sym}[MeV]$',fontsize = 30);
grid.axes[3][2].tick_params(axis='both', which='major', labelsize=21)
grid.axes[3][2].set_xlim(-350, 400)
ind=[-350,-0,400]
grid.axes[3][2].set_xticks(ind)
#label=['-200','0', '200']
#grid.axes[3][2].set_xticklabels(label)
#grid.axes[3][2].spines[:].set_visible(True)

#[3][3]
grid.axes[3][3].set_xlabel(r'$J_{sym}[MeV]$',fontsize = 30);
grid.axes[3][3].tick_params(axis='both', which='major', labelsize=23)
grid.axes[3][3].set_xlim(-400, 1000)
ind=[-400,300,1000]
grid.axes[3][3].set_xticks(ind)
grid.axes[3][3].get_yaxis().set_visible(True)
#grid.axes[3][3].spines[:].set_visible(True)

grid.savefig(f'projects/jnorhos/dolsen/NS_to_HIC_Bash/Batch_Job_{batch_job_id}/Databases/eos_{eos_num}_cornerplot.png')
