# Copyright 2024 Nanxi Yao, University of Illinois, Urbana-Champaign
#modules
import h5py
import numpy as np
from numpy import array
from scipy import interpolate
from scipy.interpolate import CubicSpline
from matplotlib import pyplot as plt
import pandas as pd
import math
import sys
from itertools import chain
from scipy.interpolate import LinearNDInterpolator
import os
import json
import jsonschema
import logging
import time
import sqlite3


logger = logging.getLogger('Convert NS to HIC')
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

#constants
n_sat = 0.16 # 1/fm^3
hc = 197.3
hbar=6.582119569*pow(10,-22) #MeV s
PI  =3.1415926
c=2.998*pow(10,23) # in fm
melectron=0.51099895/c/c #Mass of an electron in natural units
mnucleon = 938.919 # Mass of nucleon in MeV


####################################################################################
################################Calculation helpers#################################
####################################################################################

def tanh_match(x, low,high,matching_point,smoothness):#tanh matching function

    if(math.isnan(high)):
        return low

    low_part = 1/2*low*(1-np.tanh((x-matching_point)/smoothness))
    high_part = 1/2*high*(1+np.tanh((x-matching_point)/smoothness))

    return low_part + high_part





def getFiniteDer(x,y,i): #Get 1st order derivative by finite difference method

    if(i==0):
        return (y[i+1]-y[i])/(x[i+1]-x[i])
    if(i==len(x)-1):
        return (y[i]-y[i-1])/(x[i]-x[i-1])
    else:
        return (y[i+1]-y[i-1])/(x[i+1]-x[i-1])





def interpolation(x,y,n):#interpolates and finds values at nth derivative
    try:
        func = interpolate.splrep(x, y)# B spline interpolation
        ynew = interpolate.splev(x, func, der=n)#find the value of the nth deriv. of the interpolation at the x pts
    except Exception as e:
        print(f"Error in interpolation function: {e}")
        raise
    return ynew





def interpolation_single(x,y,atx):#find values of interpolation (of x,y) at a set of pts (atx)
    tck = interpolate.splrep(x, y, s=0)#interpolate B spline
    xnew = atx
    ynew = interpolate.splev(atx, tck, der=0)#find the values of the spline at the atx set of pts
    return ynew.item()





def interpolation_single_der(x,y,atx):#interpolates and finds values at points (atx) for 1st deriv
    # Ensure no invalid values (NaN or inf) are in the input arrays
    valid_indices = np.isfinite(x) & np.isfinite(y) & (x != 0) & (y != 0)
    x_valid = x[valid_indices]
    y_valid = y[valid_indices]

    if len(x_valid) == 0 or len(y_valid) == 0:
        print(f"Empty valid input arrays for x: {x} and y: {y}")
        return np.nan

    try:
        tck = interpolate.splrep(x_valid, y_valid, s=0)
        ynew = interpolate.splev(atx, tck, der=1)
    except Exception as e:
        print(f"Error in interpolation_single_der: {e}")
        print(f"Inputs x: {x}, y: {y}, atx: {atx}")
        return np.nan




    return ynew.item()




####################################################################################
############################Check saturation properties#############################
####################################################################################
def check_sat(nB,e,pressure):
    sat_pos = 0
    eovernB=np.divide(e,nB)
    for i in range(len(nB)):
        if(pressure[i]<0):
          #  print(nB[i],pressure[i],e[i], eovernB[i]-mnucleon)
            sat_pos = i
 #   print(nB[sat_pos])

    dpdnB = getFiniteDer(nB,pressure,sat_pos)
    if(0.87*n_sat > nB[sat_pos] or nB[sat_pos]>1.12*n_sat or dpdnB<22.22 ):
        return False

    elif((eovernB[sat_pos]-mnucleon)>-14 or (eovernB[sat_pos]-mnucleon)<-18):
        return False

def check_stability_causality(nB,cs2):
    for i in range(len(nB)):
        if(nB[i]>0.9*n_sat and (cs2[i]<0 or cs2[i]>1)):
            return False
    return True

def check_Jorgebound(nB,cs2):
    for i in range(len(nB)):
        if (cs2[i]>0.781):
            return False
    return True
####################################################################################
#############################Symmetry energy expansion##############################
####################################################################################

def getYQ(nB, E, L, K, J,previous_val):#get ye for B-eq
    x=(nB-n_sat)/3.0/n_sat
    esym=(E+L*x+K*x*x/2.0+J*x*x*x/6.0)/hc
   # print(E,L,K,J)
    term2=288*pow(esym,12.0)*pow(nB,2.0)+PI*PI*pow(esym,9.0)*pow(nB,3.0)
    term1=-24*pow(esym,6.0)*nB+pow(2,1.0/2.0)*pow(term2,1.0/2.0)
    term3=pow(term1,1.0/3.0)
    #print(nB,E,L,K,J,term1,term2,term3)
    yq_high=1.0/16.0*(8-pow(PI,4.0/3.0)*nB/term3/pow(2,1.0/3.0)+pow(PI/2,2.0/3.0)*term3/pow(esym,3.0))
    if(math.isnan(yq_high)): #take the last valu�| e of yQ
        return previous_val
    return yq_high

def getFermiEnergy(mass,ni): #Calculate lepton contributions to energy by Fermi energy
    kF=pow(ni/2.0*6*PI*PI,1.0/3.0)*hbar
    epsilon=pow(mass,4.0)*pow(c,5.0)/PI/PI/hbar/hbar/hbar
    x=kF/mass/c
    eden=epsilon/8.0*((2*pow(x,3.0)+x)*pow((1+x*x),1.0/2.0)-math.asinh(x));
    return eden

def getFermiPres(mass,ni): #Calculate lepton contributions to pressure by Fermi pressure
    kF=pow(ni/2.0*6*PI*PI,1.0/3.0)*hbar
    epsilon=pow(mass,4.0)*pow(c,5.0)/PI/PI/hbar/hbar/hbar
    x=kF/mass/c
    pressure=epsilon/24.0*((2*pow(x,3.0)-3*x)*pow((1+x*x),1.0/2.0)+3*math.asinh(x));
    return pressure

def geteHIC(nB,eNS,e0,L,K,J,Y_HIC): #Calcualte energy from symmetry energy expansion at a particular charge fraction and nB
    previous_val=0
    eHIC_arr=[]
    for i in range(len(nB)):
       yq = getYQ(nB[i], e0, L, K, J,previous_val)
       ne=nB[i]*yq
       previous_val=yq
       fermiE=getFermiEnergy(melectron,ne)
       expansion = e0+L*(nB[i]/n_sat-1)/3+K/18*(nB[i]/n_sat-1)*(nB[i]/n_sat-1)+J/162*(nB[i]/n_sat-1)*(nB[i]/n_sat-1)*(nB[i]/n_sat-1)
       converted_epsilon = eNS[i] - fermiE - expansion*4*((Y_HIC-yq)+(yq*yq-Y_HIC*Y_HIC))*nB[i]
       eHIC_arr.append(converted_epsilon)
    return eHIC_arr#calc complete converted energy density

def convert(nB, energy, e0, L, K, J,Y_HIC):

    try:
        eHIC_arr = geteHIC(nB,energy,e0,L,K,J,Y_HIC)
    except Exception as e:
        print(f"Error in geteHIC: {e}")

    try:
        # Prevent division by zero and handle zeros in nB
        eovernB_arr = np.divide(eHIC_arr, nB, out=np.zeros_like(eHIC_arr), where=nB!=0)

        # Replace NaN values in eovernB_arr
        eovernB_arr = np.nan_to_num(eovernB_arr, nan=0.0)
    except Exception as e:
        print(f"Error in calculating eovernB_arr: {e}")
        return None

    # Debugging output
    #print("nB:", nB)
    #print("eovernB_arr:", eovernB_arr)

    # Validate arrays before interpolation
    validate_array("nB", nB)
    validate_array("eovernB_arr", eovernB_arr)

    try:
        # Filter out zeros in nB before interpolation
        valid_indices = np.nonzero(nB)
        nB_valid = nB[valid_indices]
        eovernB_valid = eovernB_arr[valid_indices]

        pressure_arr = np.zeros_like(nB)
        pressure_valid = np.multiply(interpolation(nB_valid, eovernB_valid, 1), np.multiply(nB_valid, nB_valid))

        # Assign valid pressure values back to the original array
        pressure_arr[valid_indices] = pressure_valid

    except Exception as e:
        print(f"Error during interpolation: {e}")
        return None


    # Debugging output
    #print("pressure_arr:", pressure_arr)

    try:
        nB_sat_pos = 0
        for i in range(len(pressure_arr)):
           if pressure_arr[i] < 0:
              #print(nB[i],pressure_arr[i])
              nB_sat_pos = i


        nB_sat = nB[nB_sat_pos]
        binding_E = eovernB_arr[nB_sat_pos] - mnucleon
        K_0=interpolation_single_der(nB,pressure_arr,nB_sat)*9
        cs2_arr = []
        for i in range(len(nB)):
           cs2 = getFiniteDer(eHIC_arr,pressure_arr,i)
           cs2_arr.append(cs2)

        causality_stability_check = check_stability_causality(nB,cs2_arr)
        sym_param = [e0,L,K,J]
        jorge_check = check_Jorgebound(nB,cs2_arr)

   # converted_dic = {'Sym_Par':sym_param,'Binding_E':binding_E,'K_0':K_0,'nB':nB,'e':eHIC_arr,'p':pressure_arr,'cs2':cs2_arr,'causality_stability':causality_stability_check}
        converted_dic = {'Sym_Par':sym_param,'Binding_E':binding_E,'K_0':K_0,'nsat': nB_sat, 'causality_stability':causality_stability_check,'jorge_bound_check':jorge_check, 'Pressure':pressure_arr, 'CS2': cs2_arr, 'Energy Density':eHIC_arr, 'nB':nB}

        return converted_dic


    except Exception as e:
        print(f"Error in final calculation steps: {e}")
        return None

def validate_array(name, array):
    if np.any(np.isnan(array)):
        print(f"Warning: {name} contains NaN values")
    if np.any(np.isinf(array)):
        print(f"Warning: {name} contains Inf values")
    if len(array) == 0:
        print(f"Warning: {name} is empty")

def print_column_names(df):
    print("Column names:")
    for col in df.columns:
        print(col)







####################################################################################
#############################Get the table!##############################
####################################################################################






def convert_ns_to_hic(filename):

    try:
        crust = np.loadtxt(filename, skiprows=1).T
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    nB = np.multiply(crust[1],n_sat)
    eNS = crust[2]
    cs2NS = crust[3]
    y_HIC=0.5

    validate_array('nB', nB)
    validate_array('eNS', eNS)
    validate_array('cs2NS', cs2NS)

    #Loop parameters

    E_arr = np.arange(27,40,1)
    L_arr = np.arange(20,130,10)
    K_arr = np.arange(-250,300,50)
    J_arr = np.arange(-200,800,100)

    #Arrays to append to

    all_EOS_checks = []
    all_EOS_pressure = []
    all_EOS_cs2 = []
    all_EOS_eden = []
    all_EOS_nb = []

    row_number = 0
    eos_number = int(filename.split('_')[-1].split('.')[0])
    print(f"Processing file: {filename} with EOS number: {eos_number}")

    for esym in E_arr:
       for lsym in L_arr:
          for ksym in K_arr:
             for jsym in J_arr:
                reject = False
                for nB_0 in nB:
                   x = (nB_0-n_sat)/3.0/n_sat
                   esym_0 = (esym+lsym*x+ksym*x*x/2.0+jsym*x*x*x/6.0)/hc
                   if (esym_0 < 0 and nB_0 < 6*n_sat):
                      reject = True
                      #print(esym,lsym,ksym,jsym,esym_0)
                      break
                if reject == True:
                   continue
                else:
                   #print(esym,lsym,ksym,jsym,esym_0)
                   converted_EOS = convert(nB,eNS,esym,lsym,ksym,jsym,y_HIC)
                   if converted_EOS is None:
                       continue

                   cs2_NS_check = check_Jorgebound(nB,cs2NS)
                   converted_EOS['cs2_NS_jorge_bound_check'] = cs2_NS_check
                   sym_param = [esym, lsym, ksym, jsym]
                   all_EOS_checks.append({
                            'eos_number': eos_number,
                            'row_number': row_number,
                            'Sym_Par': sym_param,
                            'Binding_E': converted_EOS['Binding_E'],
                            'nsat': converted_EOS['nsat'],
                            'causality_stability': converted_EOS['causality_stability'],
                            'jorge_bound_check': converted_EOS['jorge_bound_check'],
                            'cs2_NS_jorge_bound_check': converted_EOS['cs2_NS_jorge_bound_check']
                   })
                   # Unpack arrays into individual elements
                   for idx, (p, cs2, eden, nb) in enumerate(zip(converted_EOS['Pressure'], converted_EOS['CS2'], converted_EOS['Energy Density'], converted_EOS['nB'])):
                       all_EOS_pressure.append({'eos_number': eos_number, 'row_number': row_number, 'pressure': p})
                       all_EOS_cs2.append({'eos_number': eos_number, 'row_number': row_number, 'cs2': cs2})
                       all_EOS_eden.append({'eos_number': eos_number, 'row_number': row_number, 'energy_density': eden})
                       all_EOS_nb.append({'eos_number': eos_number, 'row_number': row_number, 'nb': nb / n_sat})
                   row_number += 1


    if not all_EOS_checks:
        print("No valid EOS conversions found.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    df_eos_checks = pd.DataFrame(all_EOS_checks)
    df_eos_checks[['Esym', 'Lsym', 'Ksym', 'Jsym']] = pd.DataFrame(df_eos_checks['Sym_Par'].tolist(), index=df_eos_checks.index)
    df_eos_checks.drop(columns=['Sym_Par'], inplace=True)

    df_eos_checks = df_eos_checks[df_eos_checks['causality_stability'] != 0]

    df_eos_pressure = pd.DataFrame(all_EOS_pressure)
    df_eos_cs2 = pd.DataFrame(all_EOS_cs2)
    df_eos_eden = pd.DataFrame(all_EOS_eden)
    df_eos_nb = pd.DataFrame(all_EOS_nb)

    return df_eos_checks, df_eos_pressure, df_eos_cs2, df_eos_eden, df_eos_nb






def main():
    start_index = int(sys.argv[1])
    end_index = int(sys.argv[2])
    job_id = sys.argv[3]
    task_id = sys.argv[4]
    db_dir_checks = sys.argv[5]
    db_dir_pressure = sys.argv[6]
    db_dir_cs2 = sys.argv[7]
    db_dir_eden = sys.argv[8]
    db_dir_nb = sys.argv[9]

    db_path_checks = os.path.join(db_dir_checks, f'eos_checks_{job_id}_{task_id}.sqlite')
    db_path_pressure = os.path.join(db_dir_pressure, f'eos_pressure_{job_id}_{task_id}.sqlite')
    db_path_cs2 = os.path.join(db_dir_cs2, f'eos_cs2_{job_id}_{task_id}.sqlite')
    db_path_eden = os.path.join(db_dir_eden, f'eos_eden_{job_id}_{task_id}.sqlite')
    db_path_nb = os.path.join(db_dir_nb, f'eos_nb_{job_id}_{task_id}.sqlite')

    conn_checks = sqlite3.connect(db_path_checks)
    cur_checks = conn_checks.cursor()
    cur_checks.execute('''
        CREATE TABLE IF NOT EXISTS eos_checks (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            eos_number INTEGER,
            row_number INTEGER,
            esym REAL,
            lsym REAL,
            ksym REAL,
            jsym REAL,
            Binding_E REAL,
            nsat REAL,
            causality_stability INTEGER,
            jorge_bound_check INTEGER,
            cs2_NS_jorge_bound_check INTEGER
        )
    ''')

    conn_pressure = sqlite3.connect(db_path_pressure)
    cur_pressure = conn_pressure.cursor()
    cur_pressure.execute('''
        CREATE TABLE IF NOT EXISTS eos_pressure (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            eos_number INTEGER,
            row_number INTEGER,
            pressure REAL
        )
    ''')

    conn_cs2 = sqlite3.connect(db_path_cs2)
    cur_cs2 = conn_cs2.cursor()
    cur_cs2.execute('''
        CREATE TABLE IF NOT EXISTS eos_cs2 (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            eos_number INTEGER,
            row_number INTEGER,
            cs2 REAL
        )
    ''')

    conn_eden = sqlite3.connect(db_path_eden)
    cur_eden = conn_eden.cursor()
    cur_eden.execute('''
        CREATE TABLE IF NOT EXISTS eos_eden (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            eos_number INTEGER,
            row_number INTEGER,
            energy_density REAL
        )
    ''')

    conn_nb = sqlite3.connect(db_path_nb)
    cur_nb = conn_nb.cursor()
    cur_nb.execute('''
        CREATE TABLE IF NOT EXISTS eos_nb (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            eos_number INTEGER,
            row_number INTEGER,
            nb REAL
        )
    ''')

    conn_checks.commit()
    conn_pressure.commit()
    conn_cs2.commit()
    conn_eden.commit()
    conn_nb.commit()

    print(f"Start index: {start_index}")
    print(f"End index: {end_index}")

    for i in range(start_index, end_index + 1):
        try:
            filename = f"/projects/jnorhos/emilyad/Binary_Love/complete_eos/EOS_{i}.dat"
            print(f"Processing file: {filename}")


            df_eos_checks, df_eos_pressure, df_eos_cs2, df_eos_eden, df_eos_nb = convert_ns_to_hic(filename)

            if not df_eos_checks.empty:
                df_eos_checks.to_sql('eos_checks', conn_checks, if_exists='append', index=False)

            if not df_eos_pressure.empty:
                df_eos_pressure.to_sql('eos_pressure', conn_pressure, if_exists='append', index=False)

            if not df_eos_cs2.empty:
                df_eos_cs2.to_sql('eos_cs2', conn_cs2, if_exists='append', index=False)

            if not df_eos_eden.empty:
                df_eos_eden.to_sql('eos_eden', conn_eden, if_exists='append', index=False)

            if not df_eos_nb.empty:
                df_eos_nb.to_sql('eos_nb', conn_nb, if_exists='append', index=False)


        except Exception as e:
            print(f"Error processing file {i}: {e}")

    conn_checks.close()
    conn_pressure.close()
    conn_cs2.close()
    conn_eden.close()
    conn_nb.close()

if __name__ == '__main__':
    main()
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ~                                                       
