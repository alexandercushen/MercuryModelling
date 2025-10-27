#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 15:14:36 2025

@author: atcushen
"""

# Import packages

# General Purpose and Data Handling Libraries
import os
import re
import glob
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
from natsort import natsorted
import pickle
from operator import add
import random
import math

# MatPlotlib for Plotting and Visualization
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.tri as tri
import matplotlib as mpl
import matplotlib.colors as mcolors
from matplotlib import cm, ticker
from matplotlib.colors import LogNorm, LightSource, ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import LogFormatter, LogFormatterSciNotation
from matplotlib.ticker import LogLocator, MultipleLocator, NullFormatter
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import Axes3D
from streamtracer import StreamTracer, VectorGrid
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from cmap import Colormap
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec


# Scipy for Scientific Computing and Analysis
from scipy import stats, interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, griddata
from scipy.ndimage import label, gaussian_filter
from scipy.spatial import ConvexHull
from scipy.interpolate import RegularGridInterpolator
from skimage import measure
from shapely.geometry import Polygon
from scipy.ndimage import label
from scipy.ndimage import distance_transform_edt
from scipy.ndimage import binary_fill_holes
from scipy.interpolate import LSQUnivariateSpline

# Image Handling and Processing
from PIL import Image

# Tecplot for Scientific Data Visualization
import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *

# For 3d plotting
from skimage import measure

# FLEKS toolkit
import sys  
sys.path.insert(1, '/Users/atcushen/SWMF/PC/FLEKS/tools/')
import fleks,yt

# From UMGPT
from collections import defaultdict
import statistics

# Define Constants
amu = 1.67e-27
k_b = 1.38e-23
mu_0 = 1.257e-6
R_M = 2440e3 #m
m_p = 1.67e-27 # kg
e = 1.60218e-19 # C
c = 2.99e8 # m/s
eV = 6.242e18 # 1 J in eV
m_e = 9.11e-31 # kg

# Reconnection score minimums
c1_min = 0.01 
c2_min = 2.5e-7

'''User input start'''

# Declare directories for 2d and 3d datasets
folder2D = "/Users/atcushen/Documents/MercuryModelling/runs/nightside_v5_run3/alldata/"   # Should contain the FR data  in dir/FRs/
folder3D = "/Volumes/My Book Duo/runs/nightside_v5_run3/alldata/" 
output = "/Users/atcushen/Documents/MercuryModelling/analysis/figures/"

# Specify earliest data time present in folder
start_time = 120            

# Specify time range to plot (t_bound[0] used when a single time is required)
t_bound = [120.0,200.0]       

# Specify time step size 
dt = 0.05

# Specify grid spacing [R_M]
dx = 1/64   

# Specify ion:electron mass ratio
mi_me = 100

# Choose whether to actually read in data (may take a long time) or just plot from already loaded data
read_data = True

# Set the plotting axes limits
x_region = [-2.5,-1.0]
y_region = [-1.2,1.20]
z_region = [-0.3,0.7]

# Set the location to plot at (as required), or else the plane used for plotting
loc = [-1.6,0.15,0.3]

'''User input end'''

# Utility Functions

def get_files(mydir, key=".*cut_particle_region0_0.*", read_time = False, reduce = True):
    '''LEGACY CODE
    For a directory "dir", return a list of all files which match the regex expression "key"'''
    
    all_files = [f for f in listdir(mydir) if isfile(join(mydir, f))]

    files=[]
    for file in all_files:
        match = re.search(key,file)
        if match != None:
            files.append(file)
    files.sort()
    
    # Now give them the appropriate name for their time
    # If we haven't already named these files with their time, do that now
    named_files = {}
    if read_time == False:
        for i in range(len(files)):
            time = round(i*dt+start_time,3)
            named_files[time] = files[i]
    # Otherwise, read the time right from the (last 6 elements) filename
    else:
        for i in range(len(files)):
            time = str("%.2f"%float(files[i][-6:]))
            named_files[time] = files[i]
    
    # Now cut the list down to files inside t_bound
    if reduce:
        reduced_files = {}
        filtered_keys = [file_time for file_time in list(named_files.keys()) if t_bound[0] <= float(file_time) < t_bound[1]]
        for file_time in filtered_keys:
            reduced_files[file_time] = str(named_files[file_time])
        return reduced_files

    else:
        return named_files
    
def load_df_data(file = "df_data5", mydir = folder2D):
    '''Load the pickle file represeting tracked DFs
    Output is dictionary, where each item is the DF name, and contains an array of its position and computed quantities'''
    
    with open(mydir+file, 'rb') as f:
        df_data = pickle.load(f)
    return df_data
        
def round_to_dt(time):
    ''' Convert input to float, rounded to the nearest timestep, and return a formatted string '''
    try:
        value = float(time)
    except (ValueError, TypeError):
        raise ValueError("Input must be convertible to a number.")

    # Round to nearest 0.05
    rounded_final = round(value / dt) * dt

    # Format to always show at least two decimals
    return f"{rounded_final:.2f}"

def read_3d_data(time, debug = False):
    '''Read in the 3d data at specified time, from the folder_3D directory.
    Data is format as a dictionary of 3d numpy arrays, for each variable, and is stored on disk as a pickle file.
    time can be input as string or float, and will be automatically rounded to nearest timestep.
    Note that the arrays are indexed as [Y,X,Z]'''
    
    # Reformat time to appropriate string
    time = round_to_dt(time)
    
    # Read file
    file3D = str(files3D[time])
    with open(folder3D+file3D, 'rb') as f:
        if debug:
            print("reading 3d data at time = "+time+": ",str(folder3D+file3D))
        return pickle.load(f) 

def read_2d_data(time, debug = False):
    '''Read in the 2d data at specified time, from the folder2D directory.
    Data is format as a dictionary of 2d numpy arrays, for each variable, and is stored on disk as a pickle file.
    time can be input as string or float, and will be automatically rounded to nearest timestep.
    Note that the arrays are indexed as [Y,X]'''
    
    # Reformat time to appropriate string
    time = round_to_dt(time)
    
    # Read file
    file2D = str(files2D[time])
    with open(folder2D+file2D, 'rb') as f:
        if debug:
            print("reading 2d data at time = "+time+": ",str(folder2D+file2D))
        return pickle.load(f)
    
def compute_distance(p1,p2):
    ''' Return the distance between two n-dim points '''
    total = 0
    for i in range(len(p1)):
        total += (p1[i]-p2[i])**2
    return total**0.5

def get_coord_index(p,data):
    '''for a 3d point p = (x,y,z), get the cell coords (yi,xi,zi)'''
    xi = np.where(data["X"][0,:,0]>p[0])[0][0]
    yi = np.where(data["Y"][:,0,0]>p[1])[0][0]
    zi = np.where(data["Z"][0,0,:]>p[2])[0][0]
    return yi,xi,zi

def get_interpolator(X3d,Y3d,Z3d,A3d):
    ''' For a given 3d value array A, and the 3d coordinate arrays X,Y,Z, create an interpolating function '''
    return RegularGridInterpolator((X3d[0,:,0], Y3d[:,0,0], Z3d[0,0,:]), np.swapaxes(A3d,0,1), bounds_error=False, fill_value=None)

def get_2d_interpolator(X2d,Y2d,A2d):
    ''' For a given 2d value array A, and the 2d coordinate arrays X,Y, create an interpolating function '''
    return RegularGridInterpolator((X2d[0,:], Y2d[:,0]), np.swapaxes(A2d,0,1), bounds_error=False, fill_value=None)

def interpolate_onto_2d_mesh(xx,yy,zz,interp_A3d):
    '''For the 2d coordinate meshgrids xx,yy,zz, and the interpolating function interp_A3d,
    return a 2d array of the same shape as the meshgrids representing the interpolated function'''
    
    # Wrap points on the plane
    plane = np.array([xx.flatten(),yy.flatten(),zz.flatten()]).T
    
    interpolated_values_flat = interp_A3d(plane)

    # Reshape the interpolated values back to the 2D plane's shape
    return interpolated_values_flat.reshape(xx.shape)

def compute_var(data3d,var):
    '''Compute a derived quantity eg. alfven speed, for the given data3d:
        B_mag, E_mag, p_i,p _e, T_i, T_e, u_mag,i, u_mag,e, J_x, J_y, J_z, J_para, J_perp, V_A,
        M_A,i, M_A,e, c_s, V_MS, M_MS,i, M_MS,e, d_i, c1, c2, recon_score, S_integrand'''
    
    if var == "B_mag":
        # Field magnitude [nT]
        return np.sqrt(data3d['Bx']**2+data3d["By"]**2+data3d["Bz"]**2)
    
    if var == "E_mag":
        # Electric fuekd magnitude [mV/m]
        return np.sqrt(data3d['Ex']**2+data3d["Ey"]**2+data3d["Ez"]**2) * 1e-3
    
    elif var == "p_i":
        # Ion scalar pressure [nPa]
        return (data3d["pxxS1"]+data3d["pyyS1"]+data3d["pzzS1"])/3
    
    elif var == "p_e":
        # Electron scalar pressure [nPa]
        return (data3d["pxxS0"]+data3d["pyyS0"]+data3d["pzzS0"])/3
    
    elif var == "T_i":
        # Electron temperature [eV]
        return compute_var(data3d, "p_i")*1e-9 / (data3d["rhoS1"]*1e6*k_b) / 11605
    
    elif var == "T_e":
        # Electron temperature [eV]
        return compute_var(data3d, "p_e")*1e-9 / (data3d["rhoS1"]*1e6*k_b) / 11605
    
    elif var == "u_mag,i":
        # Ion bulk speed [km/s]
        return np.sqrt(data3d['uxS1']**2+data3d["uyS1"]**2+data3d["uzS1"]**2)
    
    elif var == "u_mag,e":
        # Electron bulk speed [km/s]
        return np.sqrt(data3d['uxS0']**2+data3d["uyS0"]**2+data3d["uzS0"]**2)
    
    elif var == "J_x":
        # X-component of current density [nA/m^2]
        return e*data3d['rhoS1']*1e6*(data3d['uxS1']-data3d['uxS0'])*1e3 * 1e9
    
    elif var == "J_y":
        # Y-component of current density [nA/m^2]
        return e*data3d['rhoS1']*1e6*(data3d['uyS1']-data3d['uyS0'])*1e3 * 1e9
    
    elif var == "J_z":
        # Z-component of current density [nA/m^2]
        return e*data3d['rhoS1']*1e6*(data3d['uzS1']-data3d['uzS0'])*1e3 * 1e9
    
    elif var == "J_para":
        # Field-aligned current density [nA/m^2]
        return (compute_var(data3d, "J_x") * data3d["Bx"]+compute_var(data3d, "J_y") * data3d["By"]+compute_var(data3d, "J_z") * data3d["Bz"])/compute_var(data3d, "B_mag")
    
    elif var == "J_perp":
        # (magntidue of) perpendicular current density [nA/m^2]
        return np.sqrt(compute_var(data3d, "J_x")**2+compute_var(data3d, "J_y")**2+compute_var(data3d, "J_z")**2) - np.abs(compute_var(data3d, "J_para"))
    
    elif var == "V_A":
        # Alfven speed [km/s]
        return compute_var(data3d, "B_mag")*1e-9 / np.sqrt(mu_0*data3d['rhoS1']*amu*1e6) / 1e3
    
    elif var == "M_A,i":
        # Ion alfven mach number [unitless]
        return compute_var(data3d,"u_mag,i") / compute_var(data3d, "V_A")
    
    elif var == "M_A,e":
        # Electron alfven mach number [unitless]
        return compute_var(data3d,"u_mag,e") / compute_var(data3d, "V_A")
    
    elif var == "c_s":
        # Sound speed [km/s]
        return np.sqrt(5/3 * (compute_var(data3d, "p_i")+compute_var(data3d, "p_e")) * 1e-9 / (data3d["rhoS1"]*amu*1e6)) / 1e3
    
    elif var == "V_MS":
        # Fast magnetosonic speed [km/s]
        return np.sqrt(0.5*(compute_var(data3d, "V_A")**2+compute_var(data3d, "c_s")**2) + 0.5*np.sqrt((compute_var(data3d, "V_A")**2-compute_var(data3d, "c_s")**2)**2))
    
    elif var == "M_MS,i":
        # Ion magnetosonic mach number [unitless]
        return compute_var(data3d,"u_mag,i") / compute_var(data3d, "V_MS")
    
    elif var == "M_MS,e":
        # Electron magnetosonic mach number [unitless]
        return compute_var(data3d,"u_mag,e") / compute_var(data3d, "V_MS")
    
    elif var == "d_i":
        # Ion inertial length / skin depth [R_M]
        return amu/e * np.sqrt(1/(mu_0*data3d["rhoS1"]*amu*1e6)) / R_M
    
    elif var == "c1":
        dX = R_M/64 # cell size [m]
        # First sub-term of reconnection score
        Bx = data3d["Bx"]*1e-9 # [T]
        By = data3d["Bx"]*1e-9 # [T]
        Bz = data3d["Bx"]*1e-9 # [T]
        Jx = compute_var(data3d, "J_x")*1e-9 # [A/m^2]
        Jy = compute_var(data3d, "J_y")*1e-9 # [A/m^2]
        Jz = compute_var(data3d, "J_z")*1e-9 # [A/m^2]
        epsilon = 1 
        return (Jx**2+Jy**2+Jz**2)*dX / (np.sqrt((Jy*Bz-Jz*By)**2 + (Jz*Bx-Jx*Bz)**2 + (Jx*By-Jy*Bx)**2) + np.sqrt(Jx**2+Jy**2+Jz**2)*epsilon)
    
    elif var == "c2":
        # Second sub-term of reconnection score
        dX = R_M/64 # cell size [m]
        Bx = data3d["Bx"]*1e-9 # [T]
        By = data3d["Bx"]*1e-9 # [T]
        Bz = data3d["Bx"]*1e-9 # [T]
        dB_dx = data3d["dB_dx"]*1e-9 # [T/m]
        dB_dy = data3d["dB_dy"]*1e-9 # [T/m]
        dB_dz = data3d["dB_dz"]*1e-9 # [T/m]
        dBx_dx = data3d["dBx_dx"]*1e-9 # [T/m]
        dBy_dx = data3d["dBy_dx"]*1e-9 # [T/m]
        dBz_dx = data3d["dBz_dx"]*1e-9 # [T/m]
        dBx_dy = data3d["dBx_dy"]*1e-9 # [T/m]
        dBy_dy = data3d["dBy_dy"]*1e-9 # [T/m]
        dBz_dy = data3d["dBz_dy"]*1e-9 # [T/m]
        dBx_dz = data3d["dBx_dz"]*1e-9 # [T/m]
        dBy_dz = data3d["dBy_dz"]*1e-9 # [T/m]
        dBz_dz = data3d["dBz_dz"]*1e-9 # [T/m]
        # Compute second criteria, based on curvature divergence. If large enough, it seperates x lines from o-lines or flux ropes
        B_mag = compute_var(data3d, "B_mag")*1e-9# [T]
        bx = Bx/B_mag
        by = By/B_mag
        bz = Bz/B_mag
        # We need partial derivatives of the magnetic field unit vectors.
        # To reduce the number of derivates to take, we compute them using the chain rule and the pre-computed derivatives
        dbx_dx = (1/B_mag) * (dBx_dx - bx*dB_dx)
        dbx_dy = (1/B_mag) * (dBx_dy - bx*dB_dy)
        dbx_dz = (1/B_mag) * (dBx_dz - bx*dB_dz)
        dby_dx = (1/B_mag) * (dBy_dx - by*dB_dx)
        dby_dy = (1/B_mag) * (dBy_dy - by*dB_dy)
        dby_dz = (1/B_mag) * (dBy_dz - by*dB_dz)
        dbz_dx = (1/B_mag) * (dBz_dx - bz*dB_dx)
        dbz_dy = (1/B_mag) * (dBz_dy - bz*dB_dy)
        dbz_dz = (1/B_mag) * (dBz_dz - bz*dB_dz)
        # Precompute the terms of c2
        pre_c2_x = bx*dbx_dx + by*dbx_dy + bz*dbx_dz
        pre_c2_y = bx*dby_dx + by*dby_dy + bz*dby_dz
        pre_c2_z = bx*dbz_dx + by*dbz_dy + bz*dbz_dz
        # Compute the gradient of each term
        c2_x = np.gradient(pre_c2_x,dX*R_M,axis=1)
        c2_y = np.gradient(pre_c2_y,dX*R_M,axis=0)
        c2_z = np.gradient(pre_c2_z,dX*R_M,axis=2)
        return (c2_x+c2_y+c2_z)*(dX)**2
        
    
    elif var == "recon_score":
        # MHD-AEPIC reconnection metric based off magnetic topology, from Wang et al. 2022
        c1 = compute_var(data3d, "c1")
        c2 = compute_var(data3d, "c2")
        # Return =1 where criteria are met
        recon_sites = np.zeros_like(c1)+1
        recon_sites[c1<c1_min] = 0
        recon_sites[c2<c2_min] = 0
        return recon_sites#, c1, c2
    
    elif var == "S_integrand":
        # Entropy integrand, p^gamma/B to be integrated over field line using field line tracing #[(nPa)^5/3 / nT]
        B_mag = compute_var(data3d, "B_mag")
        p = compute_var(data3d, "p_i") + compute_var(data3d, "p_e")
        gamma = 5/3
        return p**gamma / B_mag 
    
    else:
        print("ERROR: VAR\"",var,"\" NOT RECOGNIZED")
        return np.zeros_like(data3d["X"],dtype='bool')

def argmax_in_region(var,data,xlim,ylim,zlim):
    '''Find the max value of "var" in the volume bounded by the provided limits'''
    '''Return the index (iy,ix,iz), and the max value.'''
    
    if var not in data.keys():
        A = compute_var(data, var)
    else:
        A = data[var]
        
    # Get coordinate indices
    yi_min,xi_min,zi_min = get_coord_index([xlim[0],ylim[0],zlim[0]],data)
    yi_max,xi_max,zi_max = get_coord_index([xlim[1],ylim[1],zlim[1]],data)
    
    # Crop the var array
    A_crop = A[yi_min:yi_max,xi_min:xi_max,zi_min:zi_max]
    
    # Find max value of var
    flat_max_index = np.nanargmax(A_crop)
    
    # Reshape into 3d index
    max_index = np.unravel_index(flat_max_index, A_crop.shape)
    
    # Add back the cropped off indices
    max_index = max_index + np.array([yi_min,xi_min,zi_min])
    
    return max_index, A[max_index[0],max_index[1],max_index[2]]

def mean_in_region(var,data,xlim,ylim,zlim):
    '''Find the mean and standard deviation of "var" in the volume bounded by the provided limits'''
    
    if var not in data.keys():
        A = compute_var(data, var)
    else:
        A = data[var]
        
    # Get coordinate indices
    yi_min,xi_min,zi_min = get_coord_index([xlim[0],ylim[0],zlim[0]],data)
    yi_max,xi_max,zi_max = get_coord_index([xlim[1],ylim[1],zlim[1]],data)
    
    # Crop the var array
    A_crop = A[yi_min:yi_max,xi_min:xi_max,zi_min:zi_max]
    
    return np.nanmean(A_crop), np.nanstd(A_crop)

def timeseries_at_loc(t_bound,p_ls):
    '''Extract timeseries of variables in var_ls at p=(x,y,z) in t_bound)
    p_ls should be a list of (x,y,z) tuples.
    eg: timeseries_at_loc([150,151],[[-1.5,0,0.2],[-1.25,0,0.2]])
    or
    p_ls = [[-2,0,0.18],[-1.9,0,0.18],[-1.8,0,0.18],[-1.7,0,0.18],[-1.6,0,0.18],[-1.5,0,0.18],[-1.4,0,0.18],[-1.3,0,0.18],[-1.2,0,0.18]]
    data = timeseries_at_loc([120,180],p_ls)
    '''
    
    var_ls = ["Bx","By","Bz","rhoS1","T_i","T_e","uxS1","uyS1","uzS1","uxS0","uyS0","uzS0"]
    
    t_ls = np.arange(t_bound[0],t_bound[1],0.05)
    
    # Set up dictionary for saving values
    var_dict = {"t":t_ls}
    for var in var_ls:
        var_dict[var] = np.zeros((len(t_ls),len(p_ls))) # rows are times, columns are each point
    
    for itime,time in enumerate(t_ls):
        
        print("Loading time",round(time,3))
        
        data = read_3d_data(time)
        
        for ivar,var in enumerate(var_ls):
            
            if var in data.keys():
                A = data[var]
            else:
                A = compute_var(data, var)
                
            for ip,p in enumerate(p_ls):
                
                index = get_coord_index(p, data)
            
                var_dict[var][itime,ip] = A[index[0],index[1],index[2]]
    
    return var_dict

def dewey_DF_identification(t,Bx,By,Bz,ux,showPlot = False, advFilter = False):
    '''For a magnetic timeseries, apply Dewey DF identification algorithm.
    Return a list of times representing the DF rise'''
    
    dBz_dt_min = 5 # [nT/s]
    
    Delta_Bz_min = 10 # [nT]
    
    Delta_t_min = 0.4 # [s]
    
    dBz_dt = np.diff(Bz) / dt #np.gradient(Bz, dt) # [nT/s]
    delta_Bz = np.diff(Bz)
    
    dBz_dt = np.insert(dBz_dt,-1,0)
    delta_Bz = np.insert(delta_Bz,-1,0)
    
    df_ls = []
    
    df_ls_filtered = []
    
    latest_time = -1
    
    for itime,time in enumerate(t[:-1]):
        
        if time<latest_time:
            continue
        
        df_match = False
        df_match_prev = False
        df_end = False
        
        # Only look for events longer than Delta_t_min
        for jtime in range(int(Delta_t_min//dt),len(t[itime:-1])):
            
            # Save memory of whether the previous jtime was a DF
            if df_match:
                df_match_prev = True
            
            # Assume this step does not fulfill criteria
            df_match = False
            
            # Require Bz rises by enough
            if Bz[itime+jtime]-Bz[itime] > Delta_Bz_min:
                
                # Require dBz_dt is large enough
                #if min(dBz_dt[itime:itime+jtime])>dBz_dt_min:
                if (Bz[itime+jtime]-Bz[itime])/(jtime*dt) > dBz_dt_min:
                    
                    # Require Bz never dips below the Bz at itime
                    if np.min(Bz[itime+1:itime+jtime]) > Bz[itime]:
                        
                        # Require coherent gradients:
                        if np.mean(delta_Bz[itime:itime+jtime]) > np.std(delta_Bz[itime:itime+jtime]):
                            
                            # Require average planetward flow
                            #if np.mean(ux[itime:itime+jtime])>0:
                    
                            df_match = True
            
            # Check if we've come to the end of the DF
            if df_match_prev and (not df_match):
                df_end = True
        
                # Update latest time to avoid double searching
                if t[itime+jtime] > latest_time:
                    latest_time = t[itime+jtime]
                    
                break
            
        if df_end:
            print("#####DF CANDIDATE FOUND######")
            print("At t =",time,r", $\Delta B_z>$",Delta_Bz_min,"after",jtime*dt,"s")
            print(r"$\Delta B_z =",Bz[itime+jtime]-Bz[itime])
            print(r"$\delta B_z / \delta t =",(Bz[itime+jtime]-Bz[itime])/(jtime*dt))
            print(" ")
            
            df_ls.append([time,t[itime+jtime]])
    
    # Check whether to apply the more complex filtering steps            
    if advFilter:
        for df in df_ls:
            
            # Declare parameter values
            alpha = 1.75
            gamma = 1.5
            eta = 1.75
            epsilon = 1
            xi = 2
            nu = 0.3
            
            # Extract variables, named as they appear in Dewey et al. 2020
            t1 = df[0]
            t2 = df[1]
            Delta_t_DF = t2-t1
            t1i = np.where(t == t1)[0][0]
            t2i = np.where(t == t2)[0][0]
            
            # Test1 checks if the Bz rise is significantly larger than other oscillations before and after the DF
            test1 = (np.mean(Bz[np.where((t>t2) & (t<t2 + gamma*Delta_t_DF))[0]]) - np.mean(Bz[np.where((t>t1-alpha*Delta_t_DF) & (t<t1))[0]])) / np.sqrt(np.std(Bz[np.where((t>t2) & (t<t2 + gamma*Delta_t_DF))[0]])**2 + np.std(Bz[np.where((t>t1-alpha*Delta_t_DF) & (t<t1))[0]])**2)
            
            # Test2 checks if the magnetic field enhancement is long enough
            indices = np.where((t > t2) & (Bz <= np.median(Bz[t1i:t2i])))[0]
            if indices.size > 0:
                tau_2 = t[indices[-1]] - t2
            else:
                tau_2 = 0
            test2 = tau_2
            
            # Apply tests
            if test1 > eta and test2 > epsilon*Delta_t_DF:
                df_ls_filtered.append(df)
            
    else:
        df_ls_filtered = df_ls
        
    if showPlot:
        
        fig,axs = plt.subplots(figsize=(20,8),nrows=2,height_ratios=[3,1])
        
        axs[0].plot(t,Bx,color='red',label = 'Bx')
        axs[0].plot(t,By,color='green',label = 'By')
        axs[0].plot(t,Bz,color='blue',label = 'Bz')
        
        axs[1].plot(t,ux)
        
        for df in df_ls:
            
            for axi in axs:
                
                axi.axvline(x=df[0],color='red')
                axi.axvspan(df[0],df[1],color='red',alpha=0.1)
            
                axi.set_xlim(t[0],t[-1])
                
        if advFilter:
            
            for df in df_ls_filtered:
                
                for axi in axs:
                    
                    axi.axvline(x=df[0],color='blue')
                    axi.axvspan(df[0],df[1],color='blue',alpha=0.1)
        
        axs[0].legend()
        axs[0].grid()
        axs[1].grid()
        
        axs[0].set_xlabel('t [s]')
        axs[0].set_ylabel('B [nT]')
        
        axs[1].set_xlabel('t [s]')
        axs[1].set_ylabel(r'$u_{e,x}$ [km/s]')
        
        # Save figure
        fig.savefig(str(output+"dewey_DFs.png"),bbox_inches='tight',dpi=300)
                
        plt.show()
        plt.close()
    
    return df_ls_filtered

def get_point_ls(xlim,ylim,z,dx,dy):
    '''Get a pointgrid for timeseries sampling
    e.g. p_ls = get_point_ls([-3.5,-1.1],[-1.1,1.1],0.18,0.1,0.1)
    '''
    
    x_ax = np.arange(xlim[0],xlim[1],dx)
    y_ax = np.arange(ylim[0],ylim[1],dy)
    
    xx,yy = np.meshgrid(x_ax,y_ax)
    
    x_ls = np.ravel(xx)
    y_ls = np.ravel(yy)
    z_ls = np.zeros_like(x_ls)+z
    
    return np.array([x_ls,y_ls,z_ls]).T

def generate_df_heatmap(timeseries,p_ls,advFilter = True):
    '''For an array of timeseries (generated from e.g. timeseries = timeseries_at_loc([120,180],p_ls)),
    Find all of the DFs in each time series and return a heatmap of counts'''
    
    x_ax = np.unique(p_ls[:,0])
    y_ax = np.unique(p_ls[:,1])
    
    xx,yy = np.meshgrid(x_ax,y_ax)
    
    df_count = np.zeros((len(p_ls),1))
    
    for pi,p in enumerate(p_ls):
        
        t = timeseries["t"]
        Bx = timeseries["Bx"][:,pi]
        By = timeseries["By"][:,pi]
        Bz = timeseries["Bz"][:,pi]
        ux = timeseries["uxS0"][:,pi]
        
        df_ls = dewey_DF_identification(t,Bx,By,Bz,ux,advFilter = advFilter)
        
        df_count[pi] = len(df_ls)
    
    df_heatmap = np.reshape(df_count, xx.shape)
    
    return df_heatmap

def get_tracer(X,Y,Z,Ax,Ay,Az,nsteps = 10000,step_size = 1e-3, cell_size = 0.01562501):
    ''' Sets up field line tracing based off 3 3d arrays for the vectorfield '''
    ny,nx,nz = Ax.shape
    field = np.zeros((nx,ny,nz,3))
    field[:,:,:,0] = np.transpose(Ax,axes=[1,0,2])
    field[:,:,:,1] = np.transpose(Ay,axes=[1,0,2])
    field[:,:,:,2] = np.transpose(Az,axes=[1,0,2])
    grid_spacing = [cell_size,cell_size,cell_size]
    grid = VectorGrid(field, grid_spacing, origin_coord = [X.min(),Y.min(),Z.min()])
    return StreamTracer(nsteps, step_size), grid

def nearest_2d_mask(X2d, Y2d, points):
    '''For a list of points, create a mask of the single nearest cells in X2d,Y2d'''
    # Initialize a mask with False values
    mask = np.zeros(X2d.shape, dtype=bool)

    for point in points:
        # Find the index in the grid closest to this point
        idx = (
            np.abs(Y2d[:, 0] - point[1]).argmin(),
            np.abs(X2d[0, :] - point[0]).argmin()
        )
        #print("error:",np.sqrt((point[0]-X[idx])**2+(point[1]-Y[idx])**2+(point[2]-Z[idx])**2))

        # Set the mask at the nearest cell to True
        mask[idx] = True
    
    return mask

''' Plotting support functions '''

def add_colorbar(fig,plot,ax,label,shrink = 0.5, fontsize = 12):
    ''' Add colorbar to fig, ax with stated label '''
    
    clb = fig.colorbar(plot,ax=ax,shrink=shrink)
    clb.ax.set_title(label,fontsize=fontsize)
    clb.ax.tick_params(labelsize=fontsize)
    
    return clb

def add_mercury(fig,ax):
    ''' Add Mercury to the given ax'''
    core = plt.Circle((0, 0), 0.8, color='black')
    surf = plt.Circle((0, 0), 1, color='grey')
    
    ax.add_patch(surf)
    ax.add_patch(core)
    
    return fig, ax

def plot_oblique_xz_slice(p1, p2, time, var, vmin=-50, vmax=50, nlevels = 51, logScale = False, fieldlines = True, Uquivers = False, Uiquivers = False, Equivers = False, show_recon = False, cmap = 'bwr', units = 'nT', figsize = (15,9), fig = None, ax = None, savefig = True, showfig = True):
    '''Generate plot of variable 'var' at time 'time' along the plane normal to XZ and that intersects p1, p2 (where p1_x<p2_x)
    Other parameters control plot output
    Will write onto axes "ax" if it is not none'''
    
    # Read in data
    data3d = read_3d_data(time)
    
    # Unpack coordinate arrays
    X3d = data3d["X"]
    Y3d = data3d["Y"]
    Z3d = data3d["Z"]
    
    # Work out the number of cells for the length and height of the plane    
    nL = int(compute_distance(p1, p2)//dx)
    nH = int((np.max(Z3d)-np.min(Z3d))//dx)
    
    # Create x,y,z axes for plane
    x_ax = np.linspace(p1[0], p2[0], nL)
    y_ax = np.linspace(p1[1], p2[1], nL) 
    z_ax = np.linspace(np.min(Z3d), np.max(Z3d), nH)
    
    # Compute plane coordinate, s, representing distance from p_2
    s_ax = np.sqrt((x_ax-p2[0])**2+(y_ax-p2[1])**2)
    
    # Create 2d meshgrids
    xx,_ = np.meshgrid(x_ax,z_ax)
    yy,_ = np.meshgrid(y_ax,z_ax)
    ss,zz = np.meshgrid(s_ax,z_ax)
    
    # Compute axial distance meshgrid
    rhorho = np.sqrt(xx**2+yy**2)
    
    # Create an interpolating function for the variable of interest
    if var in data3d.keys():
        A3d = data3d[var]
    else:
        A3d = compute_var(data3d,var)
    interp_A3d = get_interpolator(X3d, Y3d, Z3d, A3d)
    
    # Interpolate onto the meshgrid
    A2d = interpolate_onto_2d_mesh(xx, yy, zz, interp_A3d)
    
    # Create figure
    if ax == None:
        fig,ax = plt.subplots(figsize = figsize)
    
    # Show sliced data
    if logScale:
        A2d_plot = np.log10(A2d)
    else:
        A2d_plot = A2d
    levels = np.linspace(vmin,vmax,nlevels)
    contourplot = ax.contourf(xx, zz, A2d_plot, levels = levels, cmap = cmap, extend = 'both')
    
    # Import magnetic field data if fieldlines are required
    if fieldlines:
        Bx3d = data3d["Bx"]
        By3d = data3d["By"]
        Bz3d = data3d["Bz"]
        
        # Get interpolated B components
        Bx2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, Bx3d))
        By2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, By3d))
        Bz2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, Bz3d))
        
        # Define unit vector pointing from p1 to p2
        phat = np.array([p2[0]-p1[0],p2[1]-p1[1]])/compute_distance(p1, p2)
        
        # Project Bx,By into plane
        Brho = Bx2d*phat[0] + By2d*phat[1]
        
        # Add to figure
        #ax.streamplot(ss[:,::-1], zz, -Brho[:,::-1], Bz2d[:,::-1], broken_streamlines = False, color = 'black', linewidth = 0.4, density = 2, maxlength = 2.3,arrowsize=0.5)
        ax.streamplot(xx, zz, Brho, Bz2d, broken_streamlines = False, color = 'black', linewidth = 0.2, density = 1.5, maxlength = 2.3,arrowsize=0.5)

    if Uquivers:
        ux3d = data3d["uxS0"]
        uy3d = data3d["uyS0"]
        uz3d = data3d["uzS0"]
        
        # Get interpolated B components
        ux2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, ux3d))
        uy2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, uy3d))
        uz2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, uz3d))
        
        # Define unit vector pointing from p1 to p2
        phat = np.array([p2[0]-p1[0],p2[1]-p1[1]])/compute_distance(p1, p2)
        
        # Project Bx,By into plane
        urho = ux2d*phat[0] + uy2d*phat[1]
        
        # Add to figure
        qs = 2
        ax.quiver(xx[::qs,::qs], zz[::qs,::qs], -urho[::qs,::qs], uz2d[::qs,::qs], color  ='white')
        
    if Uiquivers:
        ux3d = data3d["uxS1"]
        uy3d = data3d["uyS1"]
        uz3d = data3d["uzS1"]
        
        # Get interpolated B components
        ux2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, ux3d))
        uy2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, uy3d))
        uz2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, uz3d))
        
        # Define unit vector pointing from p1 to p2
        phat = np.array([p2[0]-p1[0],p2[1]-p1[1]])/compute_distance(p1, p2)
        
        # Project Bx,By into plane
        urho = ux2d*phat[0] + uy2d*phat[1]
        
        # Add to figure
        qs = 2
        ax.quiver(xx[::qs,::qs], zz[::qs,::qs], -urho[::qs,::qs], uz2d[::qs,::qs], color  ='white')
        
    if Equivers:
        Ex3d = data3d["Ex"]
        Ey3d = data3d["Ey"]
        Ez3d = data3d["Ez"]
        
        # Get interpolated B components
        Ex2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, Ex3d))
        Ey2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, Ey3d))
        Ez2d = interpolate_onto_2d_mesh(xx, yy, zz, get_interpolator(X3d, Y3d, Z3d, Ez3d))
        
        # Define unit vector pointing from p1 to p2
        phat = np.array([p2[0]-p1[0],p2[1]-p1[1]])/compute_distance(p1, p2)
        
        # Project Bx,By into plane
        Erho = Ex2d*phat[0] + Ey2d*phat[1]
        
        # Add to figure
        qs = 2
        ax.quiver(xx[::qs,::qs], zz[::qs,::qs], -Erho[::qs,::qs], Ez2d[::qs,::qs], color ='white')
        
    if show_recon:
        c1 = compute_var(data3d, 'c1')
        c2 = compute_var(data3d, 'c2')
        
        #Interpolate c1, c2, and then apply cutoffs (should not interpolate recon_score, since it is a binary value)
        interp_c1 = get_interpolator(X3d, Y3d, Z3d, c1)
        interp_c2 = get_interpolator(X3d, Y3d, Z3d, c2)
        
        # Interpolate onto the meshgrid
        c1_2d = interpolate_onto_2d_mesh(xx, yy, zz, interp_c1)
        c2_2d = interpolate_onto_2d_mesh(xx, yy, zz, interp_c2)
        
        # Define recon_score
        recon2d = np.zeros_like(c1_2d)+1
        recon2d[c1_2d<c1_min] = 0
        recon2d[c2_2d<c2_min] = 0
        
        ax.contour(xx,zz,recon2d, [0.9], colors = 'magenta')


    # Set axes limits
    ax.set_xlim(-1.1,p1[0])
    ax.set_ylim(*z_region)
    
    # Fix aspect ratio
    ax.set_aspect(1)
    
    # Add axes labels
    ax.set_xlabel(r"X [$R_M$]")
    ax.set_ylabel(r"Z [$R_M$]")
    
    # Add colorbar
    if showfig or savefig:
        clb = add_colorbar(fig, contourplot, ax, units)
    
    # Add grid
    ax.grid()
    
    if savefig:
        # Title
        ax.set_title(str(r"$"+var+r"$ at t = "+round_to_dt(time)))
        
        # Save figure
        fig.savefig(str(output+var+"_"+str(p1)+"_"+str(p2)+"_"+round_to_dt(time)+".png"),
                    bbox_inches='tight',dpi=300)
        
    if showfig:
        plt.show()
    
    if showfig or savefig:
        plt.close(fig)
    
    return rhorho,zz,A2d,fig,ax,contourplot
    
def plot_cs(time, var, data2d=None, vmin=-100, vmax=100, nlevels = 51, cmap = 'bwr', units = 'nT', xlim = [-4,0], ylim = [-1.2,1.2], logScale = False, savefig = False, showfig = True, figsize = (11,7), fig = None, ax = None):
    '''Plot the current sheet in projected x,y coords, shaded by variable 'var'.
    If the 2d cs file is already loaded, it can be passed to the function as "csdata" to avoid reloading data unnecessarily.
    Fig and ax objects can be passed to the function, to add the plot to. If they remain "None", a new plot will be made 
    Returns the coordinat meshgrids, variable array, and fig and ax objects
    e.g.
    plot_cs(150, 'Bz')
    '''
    
    if data2d == None:
        # Read in data
        data2d = read_2d_data(time)
        
    if fig == None or ax == None:
        fig,ax = plt.subplots(figsize = figsize)
        
    X2d = data2d["X"]
    Y2d = data2d["Y"]
        
    if var in data2d.keys():
        A2d = data2d[var]
    else:
        A2d = compute_var(data2d, var)
        
    # Add plot
    levels = np.linspace(vmin,vmax,nlevels)
    contourplot = ax.contourf(X2d,Y2d, A2d, levels = levels, cmap = cmap, extend = 'both')
    
    # Add colorbar
    clb = add_colorbar(fig,contourplot,ax,units,shrink = 0.5, fontsize = 12)
    
    # Add grid
    ax.grid()
    
    # Set aspect
    ax.set_aspect(1)
    
    # Show Mercury
    add_mercury(fig,ax)
    
    # Set axes limits
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    
    # Add axes labels
    ax.set_xlabel(r"$X_{cs}$ [$R_M$]",fontsize = 12)
    ax.set_ylabel(r"$Y_{cs}$ [$R_M$]",fontsize = 12)
    ax.tick_params(axis='both', labelsize=12)
    
    if showfig:
        plt.show()
    
        plt.close()
    
    return X2d,Y2d,A2d,fig,ax,contourplot

''' Analysis wrapper functions '''

def DF_example(t_ls,p1,p2=[0,0,0],xlim=[-1.1,-2.75]):
    '''Analysis script written to show slices of example DFs.
    Example input: for DDF, 
    DF_example([149,149.5,150],[-3,1,0])
    DF_example([149,149.5,150],[-2.8,1.2,0],p2=[-0.8,-0.1,0])
    For FRDF:
    DF_example([150.25,151.25,152.25],[-3,0.6,0])
    DF_example([150.25,151.25,152.25],[-3.3,0.6,0],p2=[-0.8,-0.1,0])'''
    
    var_ls = ["By","Bz","uxS1","uxS0","rhoS1"]
    
    # Setup figure
    fig,axs = plt.subplots(figsize = (26,20),ncols=len(t_ls),nrows=len(var_ls), sharex=True, sharey=True, gridspec_kw={'wspace': 0.01, 'hspace': 0.02})
    
    for itime,time in enumerate(t_ls):
        for ivar,var in enumerate(var_ls):
            
            if var == "By":
                vmin = -30
                vmax = 30
                cmap = 'bwr'
                label = r"$B_y$ [nT]"
                logScale = False
            elif var == "Bz":
                vmin = -60
                vmax = 60
                label = r"$B_z$ [nT]"
            elif var == "uxS1":
                vmin = -1500
                vmax = 1500
                cmap = "PuOr"
                label = r"$u_{i,x}$ [km/s]"
            elif var == "uxS0":
                vmin = -3000
                vmax = 3000
                label = r"$u_{e,x}$ [km/s]"
            elif var == "rhoS1":
                vmin = -1
                vmax = 1
                label = r"log$(n)$ [amu/cc]"
                cmap = 'rainbow'
                logScale = True
            
            print("Plotting t =",time,", and var =",var)
            _,_,_,fig,axs[ivar,itime],contourplot = plot_oblique_xz_slice(p1, p2, time, var, fig = fig, ax = axs[ivar,itime], logScale=logScale, 
                                                              vmin = vmin, vmax = vmax, cmap = cmap, savefig = False, showfig = False)
            
            axs[ivar,itime].set_xlim(*xlim)
            
            # Overwrite axes labels
            axs[ivar,itime].set_xlabel("")
            axs[ivar,itime].set_ylabel("")
            
            if ivar == 0:
                axs[ivar,itime].set_title(str("t ="+round_to_dt(time)),fontsize=18)
            if ivar == len(var_ls)-1:
                axs[ivar,itime].set_xlabel(r"X [$R_M$]",fontsize=18)
                axs[ivar,itime].tick_params(axis='x', labelsize=18)
            if itime == 0:
                axs[ivar,itime].set_ylabel(r"Z [$R_M$]",fontsize=18)
                axs[ivar,itime].tick_params(axis='y', labelsize=18)
            if itime == len(t_ls)-1:
                clb = add_colorbar(fig, contourplot, axs[ivar,:], label, shrink = 0.85, fontsize=15)
           
    # Save figure
    fig.savefig(str(output+str(t_ls)+"_"+str(p1)+"_"+str(p2)+".png"),
                bbox_inches='tight',dpi=300)
    
    plt.show()
    plt.close()
    
def DF_temperature(t_ls,p1,p2=[0,0,0],xlim=[-1.1,-2.75]):
    '''Analysis script written to show temperature of example DFs.
    Example input: for DDF, 
    DF_temperature([149,149.5,150],[-2.5,0.6667,0])
    For FRDF:
    DF_temperature([150.25,151.25,152.25],[-2.5,0.4,0])
    DF_temperature([151.50,151.75,152.00],[-3.3,0.6,0],p2=[-0.8,-0.1,0])'''
    
    var_ls = ["T_i","T_e"]
    
    # Setup figure
    fig,axs = plt.subplots(figsize = (20,12),ncols=3,nrows=2, sharex=True, gridspec_kw={'hspace': -0.3})
    
    for itime,time in enumerate(t_ls):
        for ivar,var in enumerate(var_ls):
            
            if var == "T_i":
                vmin = 0
                vmax = 6000
                cmap = 'rainbow'
                label = r"$T_i$ [eV]"
            elif var == "T_e":
                vmin = 0
                vmax = 3000
                label = r"$T_e$ [eV]"

            
            print("Plotting t =",time,", and var =",var)
            _,_,_,fig,axs[ivar,itime],contourplot = plot_oblique_xz_slice(p1, p2, time, var, fig = fig, ax = axs[ivar,itime], 
                                                              vmin = vmin, vmax = vmax, cmap = cmap, savefig = False, showfig = False)
            
            axs[ivar,itime].set_xlim(*xlim)
            
            # Overwrite axes labels
            axs[ivar,itime].set_xlabel("")
            axs[ivar,itime].set_ylabel("")
            
            if ivar == 0:
                axs[ivar,itime].set_title(str("t ="+round_to_dt(time)))
            if ivar == 1:
                axs[ivar,itime].set_xlabel(r"Axial distance $\rho$ [$R_M$]")
            if itime == 0:
                axs[ivar,itime].set_ylabel(r"Z [$R_M$]")
            if itime == 2:
                clb = add_colorbar(fig, contourplot, axs[ivar,:], label)
           
    # Save figure
    fig.savefig(str(output+str(t_ls)+"_"+str(p1)+"_"+str(p2)+".png"),
                bbox_inches='tight',dpi=300)
    
    plt.show()
    plt.close()
    
def DF_rereconnection(t_ls,p1,p2=[0,0,0],xlim=[-1.1,-2.2]):
    '''Analysis script written to show T_e, u_e, B_y, and recon_score.
    Example input: for FRDF:=
    DF_rereconnection([151.20,151.50,151.80,152.10,152.30],[-3.3,0.6,0],p2=[-0.8,-0.1,0])
    '''
    
    var_ls = ["T_i","T_e"]
    
    # Setup figure
    fig,axs = plt.subplots(figsize = (20,10),ncols=len(t_ls),nrows=len(var_ls), sharex=True, sharey = True, gridspec_kw={'wspace':0.05,'hspace': -0.3})
    
    for itime,time in enumerate(t_ls):
        for ivar,var in enumerate(var_ls):
            
            if var == "E_mag":
                vmin=0
                vmax=50
                cmap = 'viridis'
                label = r'$E$ [mV/m]'
                Equivers = True
                Uquivers = False
                show_recon = False
            if var == "c1":
                vmin = 0
                vmax = 2*c1_min
                cmap = 'rainbow'
                label = 'c1'
                Equivers = True
                Uquivers = False
                show_recon = True
            if var == "c2":
                vmin = 0
                vmax = 2*c2_min
                cmap = 'rainbow'
                label = 'c2'
                Equivers = True
                Uquivers = False
                show_recon = True
            if var == "recon_score":
                 vmin = 0
                 vmax = 1
                 cmap = 'viridis'
                 label = 'recon_score'
                 Equivers = True
                 Uquivers = False
                 show_recon = True
            if var == "J_perp":
                vmin = -250
                vmax = 250
                label = r"$J_\perp$ [nA/m$^2$]"
                cmap = 'bwr'
                Equivers = True
                Uquivers = False
                show_recon = False
            if var == "T_i":
                vmin = 0
                vmax = 5000
                label = r"$T_i$ [eV]"
                cmap = 'nipy_spectral'
                Equivers = False
                Uiquivers = True
                Uquivers = False
                show_recon = False
            if var == "T_e":
                vmin = 0
                vmax = 5000
                label = r"$T_e$ [eV]"
                cmap = 'nipy_spectral'
                Equivers = False
                Uiquivers = False
                Uquivers = True
                show_recon = False
            
            print("Plotting t =",time,", and var =",var)
            _,_,_,fig,axs[ivar,itime],contourplot = plot_oblique_xz_slice(p1, p2, time, var, fig = fig, ax = axs[ivar,itime], Equivers = Equivers, 
                                                                          Uquivers = Uquivers, Uiquivers = Uiquivers, show_recon = show_recon,
                                                                          vmin = vmin, vmax = vmax, cmap = cmap, savefig = False, showfig = False)
            
            axs[ivar,itime].set_xlim(*xlim)
            
            # Overwrite axes labels
            axs[ivar,itime].set_xlabel("")
            axs[ivar,itime].set_ylabel("")
            
            if ivar == 0:
                axs[ivar,itime].set_title(str("t ="+round_to_dt(time)))
            if ivar == len(var_ls)-1:
                axs[ivar,itime].set_xlabel(r"X [$R_M$]")
            if itime == 0:
                axs[ivar,itime].set_ylabel(r"Z [$R_M$]")
            if itime == len(t_ls)-1:
                clb = add_colorbar(fig, contourplot, axs[ivar,:], label)
           
    # Save figure
    fig.savefig(str(output+str(t_ls)+"_"+str(p1)+"_"+str(p2)+".png"),
                bbox_inches='tight',dpi=300)
    
    plt.show()
    plt.close()
            
def plot_df_heatmap(df_heatmap,p_ls,cmap='viridis'):
    '''Make plots of all the DFs found in the timeseries data.
    "df_heatmap" should be made using the generate_df_heatmap() function, 
    which needs timeseries generated using timeseries_at_loc
    '''
    x_ax = np.unique(p_ls[:,0])
    y_ax = np.unique(p_ls[:,1])
    
    # Show heatmap
    
    fig,ax = plt.subplots(figsize = (12,10))
    
    plot = ax.imshow(df_heatmap, origin = 'lower', extent = [x_ax[0],x_ax[-1],y_ax[0],y_ax[-1]],cmap = cmap)
    
    ax.grid()
    
    ax.set_xlabel(r"$X_{cs}$ [$R_M$]")
    ax.set_ylabel(r"$Y_{cs}$ [$R_M$]")
    
    ax.set_xlim(-4,-0.5)
    ax.set_ylim(-1.2,1.2)
    
    ax.set_aspect(1)
    
    clb = add_colorbar(fig, plot, ax, "DFs/min")
    
    fig,ax = add_mercury(fig, ax)
    
    return fig, ax

def mean_max_timeseries(t_bound,xlim,ylim,zlim):
    ''' Plot a timeseries of the mean and max values of variables in var list within the specified volume
    eg: mean_max_timeseries([140,160],[-1.75,-1],[-0.5,0.5],[-0.2,0.6])
    '''
    
    t_ls = np.arange(t_bound[0],t_bound[1],0.05)
    
    var_ls = ["rhoS1","T_i","T_e","J_para"]
    
    # Set up dictionary for saving values
    var_dict = {"t":t_ls}
    for var in var_ls:
        var_dict[var] = np.zeros((len(t_ls),3)) # Col 1 is max, col 2 is mean, col 3 is std
        
    for itime,time in enumerate(t_ls):
        
        print("Loading time",round(time,3))
        data3d = read_3d_data(time)
        
        for ivar,var in enumerate(var_ls):
            max_index, max_val = argmax_in_region(var, data3d, xlim, ylim, zlim)
            mean, std = mean_in_region(var, data3d, xlim, ylim, zlim)
            var_dict[var][itime,:] = np.array([max_val,mean,std])
        
    return var_dict

def plot_timeseries(var_dict, overlay_dfs = False):
    '''Plot the timeseries generated by mean_max_timeseries'''
    
    var_ls = list(var_dict.keys())
    t_ls = var_dict["t"]
    
    fig,axs = plt.subplots(figsize = (20,12), nrows = len(var_ls)-1, ncols = 1, sharex=True, gridspec_kw={'wspace': 0.01, 'hspace': 0.02})
    
    for i,var in enumerate(var_ls[1:]): 
        
        max_vals = var_dict[var][:,0]
        mean_vals = var_dict[var][:,1]
        axs[i].plot(t_ls,max_vals,color='red')
        axs[i].plot(t_ls,mean_vals,color='black')
        #axs[i].fill_between(t_ls,var_dict[var][:,1]-var_dict[var][:,2],var_dict[var][:,1]+var_dict[var][:,2],color='black',alpha=0.1)
        
        axs[i].set_xlim(t_ls[0],t_ls[-1])
        #axs[i].set_ylim(0.1 * np.min(mean_vals),1.1*np.max(max_vals))
        
        #axs[i].set_yscale('log')
        
        if var == 'rhoS1':
            axs[i].set_ylabel(r"$n$ [amu/cc]")
        elif var == 'T_i':
            axs[i].set_ylabel(r"$T_i$ [eV]")
        elif var == 'T_e':
            axs[i].set_ylabel(r"$T_e$ [eV]")
        elif var == 'J_para':
            axs[i].set_ylabel(r"$J_\parallel$ [nA/m$^2$]")
   
    axs[i].set_xlabel("t [s]")
    
    if overlay_dfs:
        
        min_df_length = 20
        
        # Load dfs
        df_data = load_df_data()
        
        # Go through each
        for key in df_data.keys():
            
            df = df_data[key]
            
            # Filter the dfs
            if len(df)>min_df_length and np.mean(df["area"])*64**2>40 and np.max(np.abs(np.diff(df['time'])))<1 and np.median(2*(df["Bz"]*df['area'])**0.6)>0 and np.mean(df["X"].diff())>0:
                
                print("Plotting DF #",key,"from t =",df['time'][0])
                
                # Loop through each plotted variable and see if its saved in the df data
                for i,var in enumerate(var_ls[1:]):  
                    
                    df_var = None
                    scale = 1
                    
                    if var == 'rhoS1':
                        df_var = 'n'
                        scale = 1e6
                    elif var == 'T_i':
                        df_var = 'Ti'
                        scale = 1e-3
                    elif var == 'T_e':
                        df_var = 'Te'
                        scale = 1e-3
                    if df_var != None:
                        
                        axs[i].plot(df['time'],df[df_var]/scale,color='blue')
            
    plt.show()
    plt.close()

def make_movie(t_start,t_stop,functionName,funcInput1=None, funcInput2=None, funcInput3=None, funcInput4=None, funcInput5=None, funcInput6=None, funcInput7=None):
    '''Loop over the time range, calling "functionName" to make many plots
    eg make_movie(150,153,"plot_oblique_xz_slice",funcInput1=[-2.4,0.35,0],funcInput12=[-0.8,-0.1,0],funcInput3='T_i',funcInput4=0,funcInput5=6000,funcInput6='rainbow',funcInput7 = r'$T_i$ [eV]')'''

    for itime, time in enumerate(np.arange(t_start,t_stop,dt)):
        
        if functionName == "plot_oblique_xz_slice":
            # funcInput1 is p1, funcInput2 is p2, funcInput3 is var, 4 is vmin, 5 is vmax, 6 is cmap, 7 is unit label
            _,_,_,_,_,_ = plot_oblique_xz_slice(funcInput1, funcInput2, time, funcInput3, vmin = funcInput4, vmax = funcInput5, 
                                                cmap = funcInput6, units = funcInput7, showfig = False)

def plot_tracked_DFs(df_heatmap,p_ls,debug=False):
        
    min_df_length = 30# time steps
    
    # Load tracked DFs
    df_data = load_df_data()

    # Create figure
    fig,ax = plot_df_heatmap(df_heatmap,p_ls,cmap='viridis')
    my_dict = {}

    FR_DF_count = 0
    FR_DF_R = 0
    D_DF_count = 0
    D_DF_R = 0
    
    for key, value in df_data.items():
        if len(value)>min_df_length and np.mean(value["area"])*64**2>40 and np.max(np.abs(np.diff(value['time'])))<1 and np.median(2*(value["Bz"]*value['area'])**0.6)>0 and np.mean(value["X"].diff())>0:
        #if key==3698:
            # Create temporary df
            temp = value.copy()
            # Choose the order of the spline
            k = 3 # cubic spline
            # Choose interior knots, you can use np.linspace for uniformly distributed knots
            knots = np.linspace(value['time'].min(), value['time'].max(), 4)[1:-1]
            # Fit splines to each dimension
            spline_x = LSQUnivariateSpline(value['time'], value['X'], knots, k=k)
            spline_y = LSQUnivariateSpline(value['time'], value['Y'], knots, k=k)
            spline_z = LSQUnivariateSpline(value['time'], value['Z'], knots, k=k)
            # Evaluate the splines over a fine grid
            temp['X_smooth'] = spline_x(value['time'])
            temp['Y_smooth'] = spline_y(value['time'])
            temp['Z_smooth'] = spline_z(value['time'])
            z1 = value["frac_closed"]-value["frac_external"]
            z2 = value["frac_ropey"]
            a = (value["Bz"]*value['area'])**0.6*2
            #a = df_data[key]["wake_Jy"]
            # Create a set of line segments so that each segment is colored separately
            points = np.array([temp['X_smooth'], temp['Y_smooth']]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Create a LineCollection for the topology segments
            #cmap = Colormap('PRGn_10_r') .to_mpl()
            cmap = plt.get_cmap('coolwarm_r')#('winter') 
            norm = plt.Normalize(0, 1)
            lc2 = LineCollection(segments, cmap=cmap, norm=norm)
           # Set the values used for colormapping
            lc2.set_array(z2)
    
            # Set the widths
            lc2.set_linewidths(2)
    
            # Show line segments
            ax.add_collection(lc2)
    
            # Save data
            my_dict[key] = temp

            # Update counters for planetward moving guys
            if np.mean(value["X"].diff())>0: 
                if np.max(value['frac_ropey'])>0.5:
                    FR_DF_count += 1
                    FR_DF_R += np.sqrt(value["X"][0]**2+value["Y"][0]**2+value["Z"][0]**2)
                    if debug:
                        ax.scatter(temp['X_smooth'][0],temp['Y_smooth'][0],marker="*",color='blue')
                        #ax.text(temp['X_smooth'][0], temp['Y_smooth'][0],str(key), fontsize = 10, color='blue')
                        print("Found FRDF#",key,"starting at t =",value["time"][0])
                else:
                    D_DF_count += 1
                    D_DF_R += np.sqrt(value["X"][0]**2+value["Y"][0]**2+value["Z"][0]**2)
                    if debug:
                        ax.scatter(temp['X_smooth'][0],temp['Y_smooth'][0],marker="*",color='red')
                        #ax.text(temp['X_smooth'][0], temp['Y_smooth'][0], str(key), fontsize = 10, color='red')
                        print("Found DDF#",key,"starting at t =",value["time"][0])
    
    # Reformat the df_data object into a single dataframe
    my_df = pd.concat(my_dict.values(), axis=0)
    my_df.reset_index(drop=True, inplace=True)

    #plot1 = gaussian_filter(avgs["Bz"],sigma=5)
    #contour=ax.contour(avgs["X"],avgs["Y"],plot1,[0],zorder=0.1,colors=['black'],linestyles='dashed')

    # Show quivers
    #qskip = 6
    #ue_quiver = ax.quiver(avgs['X'][::qskip,::qskip],avgs['Y'][::qskip,::qskip],avgs['uxS0'][::qskip,::qskip],avgs['uyS0'][::qskip,::qskip],color='black')

    ax.set_title(str("Field line topology\n DDFs ($R_0\simeq$"+str(round(D_DF_R/D_DF_count,2))+"): "+str(D_DF_count)+"  FRDFs ($R_0\simeq$"+str(round(FR_DF_R/FR_DF_count,2))+"): "+str(FR_DF_count)),fontsize=15)
    
    # Custom colorbar
    cbar1 = plt.colorbar(lc2, ax=ax,shrink=0.5)   
    cbar1.ax.tick_params(labelsize=15)
    cbar1.set_ticks([0,1])
    cbar1.set_ticklabels(['Unlooped','Looped'],rotation='vertical')
    cbar1.ax.get_yticklabels()[0].set_va('bottom')
    cbar1.ax.get_yticklabels()[1].set_va('top')
    
    ax.set_xlim(-2.5,-0.6)
    ax.set_ulim(-1.1,1.1)
    ax.set_aspect(1)
    
    fig.savefig(str(output+"df_heatmap_with_tracked_dfs.png"),bbox_inches='tight',dpi=300)

def parcel_adiabaticity(t0, delta_t, p, r = 1/32, xlim = [-2.5,-0.5], debug = False):
    ''' Starting at point p = (x,y) (in current sheet) at time  = t0, track all cells within r backwards in time according to the electron velocity.
    Follow this plasma parcel backwards for delta_t seconds, computing quantities of interest.
    e.g.
    data_log, var_ls = parcel_adiabaticity(153,3,[-1.3,0])
    '''
    
    var_ls = ["Bx","By","Bz","rhoS1","T_i","T_e","uxS1","uyS1","uzS1","uxS0","uyS0","uzS0","S"]
    
    for itime, time in enumerate(np.arange(t0,t0-delta_t,-dt)):
        
        time = float(round_to_dt(time))
        print("#############################\nt=",time)
        
        # Load data
        data2d = read_2d_data(time)
        
        fig,ax = plt.subplots(figsize=(8,8))
        
        # Plot the cs
        X2d,Y2d,A2d,fig,ax,contourplot = plot_cs(time, "Bz", data2d=data2d, fig=fig, ax=ax, showfig=False, xlim = xlim, vmin=-150,vmax=150)
        
        # Get the current sheet z coords
        Z2d = data2d["Z"]
        z_interp = get_2d_interpolator(X2d, Y2d, Z2d)
        
        
        
        # Initialize the tracked cells for itime=0
        if itime==0:
            cells = np.sqrt((X2d - p[0])**2 + (Y2d - p[1])**2) <= r
            points = np.array([X2d[cells],Y2d[cells],Z2d[cells]]).T
            
            # Data format: a dictionary, with each entry being a point's time history
            data_log = {}
            for point in range(len(points)):
                data_log[point] = np.zeros((len(np.arange(t0,t0-delta_t,-dt)), len(var_ls)+4))
                data_log[point][0,0] = time
                data_log[point][0,1] = X2d[cells][point]
                data_log[point][0,2] = Y2d[cells][point]
                data_log[point][0,3] = data2d["Z"][cells][point]
                
            if debug:
                print("Initial points (x,y,z):")
                for point in data_log.keys():
                    print("point",point,":",data_log[point][0,1:4])
                    
        # At all other times, get the points and cells from the data_log
        else:
            for ipoint, point in enumerate(points):
                
                # Update the current sheet position for all the points at this time
                data_log[ipoint][itime,0] = time
                data_log[ipoint][itime,3] = z_interp( data_log[ipoint][itime,1:3])
                points[ipoint,:] = data_log[ipoint][itime,1:4]
                
                # Update the cell masks
                cells = nearest_2d_mask(X2d, Y2d, points)
                
            if debug:
                print("points (x,y,z) at t +",time,":")
                print(points)
                
        if debug:
            print("Point list center:",np.mean(points[:,0]),np.mean(points[:,1]))
            print("Cell mask center (should be similar):",np.mean(X2d[cells]),np.mean(Y2d[cells]))
                
            
                
        # At it=itime, we interpolate the values of each item in var_ls to each point using the data at t=time
        for ivar,var in enumerate(var_ls):
            
            # Special case for fieldline-integrated quantities
            if var == "S":
                
                # Load the 3d data
                data3d = read_3d_data(time)
                X3d = data3d["X"]
                Y3d = data3d["Y"]
                Z3d = data3d["Z"]
                Bx3d = data3d["Bx"]
                By3d = data3d["By"]
                Bz3d = data3d["Bz"]
                S_integrand = compute_var(data3d, "S_integrand")
                
                # Set up field line tracing
                tracer,grid = get_tracer(X3d,Y3d,Z3d,Bx3d,By3d,Bz3d,nsteps = 4000,step_size = 5e-4)
                
                # Get interpolating function
                S_interpolator = get_interpolator(X3d, Y3d, Z3d, S_integrand)
                
                # Trace field lines from points
                tracer.trace(points, grid)
                
                # Interate through each trace
                for ipoint, point in enumerate(points):
                    
                    # Get x,y,z
                    trace = tracer.xs[ipoint]
                    trace_x = trace[:,0]
                    trace_y = trace[:,1]
                    trace_z = trace[:,2]
            
                    # Filter to closed field lines
                    if (np.sum(trace[0,:]**2)<2**2) and (np.sum(trace[-1,:]**2)<2**2):
                        
                        # INTERPOLATE: interpolate the volume value at each point along the field line
                        # Compute distances between points
                        diffs = np.diff(trace, axis=0)
                        l = np.sqrt(np.sum(diffs**2, axis=1))
                        S_interp = np.nansum(S_interpolator(trace))
                        S_integrated = np.sum(l*S_interp) #[(nPa)^5/3 R_M / nT]
                        
                        # Save this result to the correct point, at time itime, to the column corresponding to this variable (plus 4, since the first 4 places are t,x,y,z)
                        data_log[ipoint][itime,ivar+4] = S_integrated
            
            # All other variables can just be pulled from the 2d array data.
            else:
                if var in data2d.keys():
                    var_array = data2d[var]
                    
                else:
                    var_array = compute_var(data2d, var)
                    
                # Get interpolating function
                var_interp = get_2d_interpolator(X2d, Y2d, var_array)
                
                # Save interpolated value for each point
                for ipoint, point in enumerate(points):
                    
                    #Save this result to the correct point, at time itime, to the column corresponding to this variable (plus 4, since the first 4 places are t,x,y,z)
                    data_log[ipoint][itime,ivar+4] = var_interp(point[0:2])
            
                    
        
        # Propogate points for the next time step using the electron bulk velocity, for all times except the last one
        if itime<len(data_log[0])-1:
            
            uex = data2d['uxS0']
            uey = data2d['uyS0']
            
            # Get u_e interpolating functions to get electron velocity at each point
            uex_interp = get_2d_interpolator(X2d, Y2d, uex)
            uey_interp = get_2d_interpolator(X2d, Y2d, uey)
            
            #if debug:
            #    print("Interpolation test: error between cell and interpolated point values of uex:",uex[cells] - uex_interp(points[:,:2]))

            # propogate each point using interpolated velocity
            for ipoint, point in enumerate(points):
                
                # Determine the velocity for this point, in R_M/s
                iuex = uex_interp(point[0:2])*1e3/R_M
                iuey = uey_interp(point[0:2])*1e3/R_M
                
                # Compute the change in positiion
                idx = -iuex*dt
                idy = -iuey*dt
                
                # Compute the new positions
                new_x = point[0] + idx
                new_y = point[1] + idy
                
                # Update x and y coords for this point at the next time
                data_log[ipoint][itime+1,1] =  new_x
                data_log[ipoint][itime+1,2] =  new_y
                
                # FWe will get the z coord at the next time step.
        
        
        
        # Update the diagnostic plot
        cells_plot = np.zeros_like(X2d)
        cells_plot[cells] = 1
        ax.contour(X2d,Y2d,cells_plot,[0.99], colors = ['black'], linewidths = [0.5], zorder = 100)
        ax.scatter(points[:,0],points[:,1],marker='x',color='lime',zorder=200)
        
        # Show the bulk reverse velocity vector
        ax.quiver(np.mean(X2d[cells]),np.mean(Y2d[cells]), -np.mean(uex[cells]), -np.mean(uey[cells]))
        
        ax.set_title(str("T = "+str(time)+"    Parcel originating at: "+str(p)))
        
        plt.show()
        plt.close()
        
        
    return [data_log, var_ls]

def plot_parcel_adiabaticity(data_log, var_ls):
    ''' Having run parcel_adiabaticity to get data_log and var_ls, now show the actual timeseries data for the parcel in question'''
    
    fig,axs = plt.subplots(figsize=(20,16),nrows = 6, sharex=True, gridspec_kw={'height_ratios': [1,3,2,2,2,2], 'wspace': 0.01, 'hspace': 0.02})
    
    axs2_2 = axs[2].twinx()
    
    for ipoint,data in data_log.items():
        
        t = np.arange(153,150,-dt)#data[:,0]
        
        print(t)
        
        axs[0].set_ylabel(r"$X,Y,Z$ [$R_M$]", fontsize=15)
        axs[0].plot(t, data[:,1], color = 'red')
        axs[0].plot(t, data[:,2], color = 'green')
        axs[0].plot(t, data[:,3], color = 'blue')
        
        axs[1].set_ylabel(r"B [nT]", fontsize=15)
        axs[1].plot(t, data[:,var_ls.index("Bx")+4], color = 'red')
        axs[1].plot(t, data[:,var_ls.index("By")+4], color = 'green')
        axs[1].plot(t, data[:,var_ls.index("Bz")+4], color = 'blue')
        
        axs[2].set_ylabel("T [eV]", fontsize=15)
        axs[2].plot(t, data[:,var_ls.index("T_i")+4],color='tab:blue')
        axs[2].plot(t, data[:,var_ls.index("T_e")+4],color='tab:orange')
        
        axs2_2.set_ylabel("n [amu/cc]", fontsize=15)
        axs2_2.plot(t, data[:,var_ls.index("rhoS1")+4], color = 'black')
        
        axs[3].set_ylabel(r"$u_i$ [km/s]", fontsize=15)
        axs[3].plot(t, data[:,var_ls.index("uxS1")+4], color = 'red')
        axs[3].plot(t, data[:,var_ls.index("uyS1")+4], color = 'green')
        axs[3].plot(t, data[:,var_ls.index("uzS1")+4], color = 'blue')
        
        axs[4].set_ylabel(r"$u_e$ [km/s]", fontsize=15)
        axs[4].plot(t, data[:,var_ls.index("uxS0")+4], color = 'red')
        axs[4].plot(t, data[:,var_ls.index("uyS0")+4], color = 'green')
        axs[4].plot(t, data[:,var_ls.index("uzS0")+4], color = 'blue')
        
        axs[5].set_ylabel(r"Entropy ([(nPa)$^{5/3}/R_M$]", fontsize=15)
        axs[5].set_yscale("log")
        axs[5].plot(t, data[:,var_ls.index("S")+4],color='black')
        
    for ax in axs:
        ax.grid()
        ax.set_xlim(t[-1],t[0])
        ax.tick_params(axis='both', labelsize=14)
        
    axs[-1].set_xlabel("Time [s]")
    
    axs[0].set_title(str("Lagrangian parcel ending at: "+str(data_log[0][0,1:4])), fontsize = 15)
    
    fig.savefig(str(output+"parcel_adiabaticity.png"),bbox_inches='tight',dpi=300)
    
    plt.show()
    plt.close()

# Load all the relevant file names
files3D = get_files(folder3D,key="3d\_fluid.*numpy\_t\_...\...",read_time = True)
files2D = get_files(folder2D,key="3d\_fluid.*csdata\_t\_...\...",read_time = True)

