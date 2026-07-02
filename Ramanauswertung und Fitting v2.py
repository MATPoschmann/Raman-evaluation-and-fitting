# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 2026
requires:
Python 3.12    
pandas: 2.2.2
matplotlib: 3.9.2
scipy: 1.14.0
or compatible
@author: poschmann

Spectral Data Processing and Visualization Script

This script provides a set of utility functions to read, process, and visualize Raman spectral data of Carbon materials. 

It includes:
- Interactive user input functions with type validation and exit options.
- A plotting function that visualizes the mean spectrum with shaded standard deviation,
  normalizing intensities and saving the plot as a high-resolution PNG file.
- function to find peaks in the mean spectrum using the `scipy.signal.find_peaks` method, with customizable parameters for peak detection.
- Fitting functions for Lorentzian and Breit-Wigner-Fano profiles, as well as combined fitting functions for multiple peaks.

The script is designed to support data analysis workflows in materials science or related fields,
enabling robust handling of spectral files and intuitive visualization of spectral characteristics.

Dependencies:
- pandas
- matplotlib
- sys (for program exit)

Usage:
- Customize default calibration parameters as needed.
- Use the input functions to safely gather user parameters.
- Call `header_info_spectra` to read spectral metadata.
- Use `plot_and_save` to generate and save visual summaries of spectral data.


"""
#fontsize in grafics
axislabelsize=20
axisticksize=20
legendsize=16
# needed python packages
import os
import pandas as pd
import glob
#from scipy.signal import argrelextrema
from scipy.signal import find_peaks
#import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
import json



plt.rcParams.update({'font.size': 10})

#initial Peak parameters
peak_width_init = 20
Lorenz_height_init = 50
BWF_height_init = 0.5
BWF_asymmetry_init = 0.02
base_BWF_init = 0
D_prime_center_init = 1610

#peak parameter limits
peak_width_limits = [1, 250]
Lorenz_height_limits = [30, 10000]
BWF_height_limits = [0, 5]
BWF_asymmetry_limits = [0, 10]
D_doubleprime_center_limits = [1050, 1250]
D_prime_center_limits = [1605, 1630]
D_center_limits = [1300, 1400]
G_center_limits = [1550, 1610]


#fitting functions
def fourpeak_3L_BWF(x, L_width1, L_width2, L_width3, widthBWF, asymBWF, L_center1, L_center2, L_center3, L_height1, L_height2, L_height3, centerBWF, heightBWF, baseBWF = 0):
    return (L_height1*  (1/ ((x - L_center1)**2 + (L_width1 / 2)**2)))+ (L_height2*  (1/ ((x - L_center2)**2 + (L_width2 / 2)**2)))+ (L_height3*  (1/ ((x - L_center3)**2 + (L_width3 / 2)**2))) + baseBWF + (heightBWF * ((1 + (centerBWF-x)/ (asymBWF * widthBWF))**2 / (1 + ((centerBWF-x)/widthBWF)**2)))
def threepeak_2L_BWF(x, L_width1, L_width2, widthBWF, asymBWF, L_center1, L_center2, L_height1, L_height2, centerBWF, heightBWF, baseBWF = 0):
    return (L_height1*  (1/ ((x - L_center1)**2 + (L_width1 / 2)**2)))+ (L_height2*  (1/ ((x - L_center2)**2 + (L_width2 / 2)**2))) + baseBWF + (heightBWF * ((1 + (centerBWF-x)/ (asymBWF * widthBWF))**2 / (1 + ((centerBWF-x)/widthBWF)**2)))
def lorenzandBWF(x, width, widthBWF, asymBWF, center, height, centerBWF, heightBWF, baseBWF = 0):
    return (height*  (1/ ((x - center)**2 + (width / 2)**2))) + baseBWF + (heightBWF * ((1 + (centerBWF-x)/ (asymBWF * widthBWF))**2 / (1 + ((centerBWF-x)/widthBWF)**2)))
def lorenz(x, width, center, height):
    return height * (1/ ((x - center)**2 + (width / 2)**2))
def BWF(x, widthBWF, asymBWF, centerBWF , heightBWF, baseBWF):
   return baseBWF + heightBWF * ((1 + (centerBWF-x)/ (asymBWF * widthBWF))**2 / (1 + ((centerBWF-x)/widthBWF)**2))

def remember_path():
    """
    Reads the last used path from a text file and returns it.
    
    Returns:
        str: The last used path as a string. If the file does not exist, returns an empty string.
    """
    try:
        path = os.path.realpath(__file__).rsplit("\\",1)[0]+'\\'
        with open(f'{path}last_used_paths.json', 'r') as f:
            json_dump = json.load(f)
            return json_dump.get(__file__, 'C:\\')  # Return the last used path or default if not found
    except (FileNotFoundError, json.JSONDecodeError):
        return ''  # Return empty string if file does not exist or is invalid

def save_path(pathdata_folder):
    """
    Saves the given path to a text file for future reference.
    
    Args:
        pathdata_folder (str): The path to save.
    """
    path = os.path.realpath(__file__).rsplit("\\",1)[0]+'\\'
    try:
        with open(f'{path}last_used_paths.json', 'r') as f:
            program_list = json.load(f)
            program_list[__file__] = pathdata_folder
        with open(f'{path}last_used_paths.json', 'w') as f:
            json.dump(program_list, f)
    except (FileNotFoundError, json.JSONDecodeError):
        with open(f'{path}last_used_paths.json', 'w') as f:
            json.dump({__file__: pathdata_folder}, f)


def select_folder(Titletext, initial_path):
    """Opens a folder selection dialog and allows the user to choose a folder.

    This function creates a hidden Tkinter root window and displays a native folder dialog.
    The selected folder path can then be used for further processing.

    Features:
    - Native folder dialog (operating system-specific).
    - Returns the folder path as a string or `None` if canceled.

    Args:
        None

    Returns:
        str:
            The absolute path of the selected folder as a string.
            Returns `None` if the user cancels the dialog.
            Returns seperate info if a folder is selected.

    Raises:
        None

    Examples:
        >>> selected_folder = select_folder()
        >>> if selected_folder:
        ...     print(f"Selected folder: {selected_folder}")
        Selected folder: /home/user/data.csv

    Notes:
        - The dialog does not support custom descriptions (e.g., as a label).
          Use the `title` parameter in the dialog for user instructions.
        - For a more detailed UI with descriptions, create a custom Tkinter window.
    """
    initial_path = remember_path()
    # Tkinter-Hauptfenster erstellen (wird für filedialog benötigt)
    root = tk.Tk()
    root.withdraw()  # Fenster verstecken, da wir nur den Dialog brauchen

    # Dateiauswahldialog anzeigen
    file_path = filedialog.askdirectory(
        title=Titletext,
        initialdir=initial_path,  # Standardmäßig das Wurzelverzeichnis
        mustexist=True
    )

    # Überprüfen, ob eine Datei ausgewählt wurde
    if file_path:
        print(f"Selected Folder: {file_path}")  # Optional für die Konsole
        file_true = True
        save_path(file_path)
    else:
        messagebox.showinfo("Break", "No folder selected.")
        file_true = False
    return file_path, file_true

def four_peak_fit_and_save(fitting_region):
    #least squares fitting algorithm 
    start_parameters =  [peak_width_init, 
                         peak_width_init, 
                         peak_width_init, 
                         peak_width_init, 
                         BWF_asymmetry_init, 
                         Peak_D['Wavenumber'].values[0]-220, 
                         Peak_D['Wavenumber'].values[0], 
                         D_prime_center_init,
                         Lorenz_height_init, 
                         Lorenz_height_init, 
                         Lorenz_height_init,
                         Peak_G['Wavenumber'].values[0], 
                         BWF_height_init, 
                         base_BWF_init] 
    min_bounds =        [peak_width_limits[0], 
                         peak_width_limits[0],
                         peak_width_limits[0],
                         peak_width_limits[0], 
                         BWF_asymmetry_limits[0],                     
                         D_doubleprime_center_limits[0],              
                         D_center_limits[0], 
                         D_prime_center_limits[0],
                         Lorenz_height_limits[0],    
                         Lorenz_height_limits[0],
                         Lorenz_height_limits[0],
                         G_center_limits[0],   
                         BWF_height_limits[0], 
                         base_BWF_init] 
    max_bounds =        [peak_width_limits[1], 
                         peak_width_limits[1],
                         peak_width_limits[1],
                         peak_width_limits[1], 
                         BWF_asymmetry_limits[1],                     
                         D_doubleprime_center_limits[1],              
                         D_center_limits[1], 
                         D_prime_center_limits[1],
                         Lorenz_height_limits[1],    
                         Lorenz_height_limits[1],
                         Lorenz_height_limits[1],
                         G_center_limits[1],   
                         BWF_height_limits[1],  
                         base_BWF_init + 0.0001]
    fit_params, pcov = scipy.optimize.curve_fit(fourpeak_3L_BWF, 
                                                fitting_region['Wavenumber'], 
                                                fitting_region['Mean_Intensity'], 
                                                p0 = start_parameters, 
                                                bounds = (min_bounds, max_bounds), 
                                                method = 'trf') 
    #generating fitting graph
    plt.close()
    plt.figure(figsize=(10,8))
    #mean Raman spectrum
    plt.plot(fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], label = 'Mean Intensity')
    #D-sideband curve
    plt.plot(fitting_region['Wavenumber'], lorenz(fitting_region['Wavenumber'], fit_params[0], fit_params[5], fit_params[8]), label = 'D" fit')
    #D-band curve
    plt.plot(fitting_region['Wavenumber'], lorenz(fitting_region['Wavenumber'], fit_params[1], fit_params[6], fit_params[9]), label = 'D fit')
    #D-band curve
    plt.plot(fitting_region['Wavenumber'], lorenz(fitting_region['Wavenumber'], fit_params[2], fit_params[7], fit_params[10]), label = "D' fit")
    #G-band curve
    plt.plot(fitting_region['Wavenumber'], BWF(fitting_region['Wavenumber'], fit_params[3], fit_params[4], fit_params[11], fit_params[12], fit_params[13]), label = 'G fit')
    #fit convolution
    plt.plot(fitting_region['Wavenumber'], fourpeak_3L_BWF(fitting_region['Wavenumber'], fit_params[0], fit_params[1], fit_params[2], fit_params[3], fit_params[4], fit_params[5], fit_params[6], fit_params[7], fit_params[8], fit_params[9], fit_params[10], fit_params[11], fit_params[12], fit_params[13]), label = 'Convolution Curve')
    #axis labels
    plt.xlabel(r'Wavenumber /$cm^{-1}$')
    plt.ylabel('Normalized Intensity /a. u.')
    #legend
    plt.legend(loc = 'upper left')
    #plot range
    plt.xlim([900, 1900])
    plt.ylim([-0.05, 1.2])
    #peak labeling
    plt.text(fit_params[5], 1.0, 'D"-Peak\n(Lorenz)', horizontalalignment='center')
    plt.text(fit_params[6], 1.1, 'D-Peak\n(Lorenz)', horizontalalignment='center')
    plt.text(fit_params[7], 1.0, "D'-Peak\n(Lorenz)", horizontalalignment='center')
    plt.text(fit_params[11], 1.1, 'G-Peak\n(BWF)', horizontalalignment='center')
    #plotsaving
    plt.savefig(samplename + 'Mean_Raman_data_fit.png', dpi = 300)
    #plot showing
    plt.show()
    #save fitting data to file
    with open(samplename + 'fit_data.txt', 'w') as fit_data_file:
        #D"-band data
        fit_data_file.writelines('D"-Band(Lorenzian): y = height / ((x - center)^2 + (width / 2)^2) = '+ str(round(fit_params[8],2)) + ' / ((x- ' + str(round(fit_params[5],2)) + ')^2 + (' + str(round(fit_params[0],2)) + '/2)^2)')
        #D-band data
        fit_data_file.writelines('\nD-Band (Lorenzian)): y = height / ((x - center)^2 + (width / 2)^2) = '+ str(round(fit_params[9],2)) + ' / ((x- ' + str(round(fit_params[6],2)) + ')^2 + (' + str(round(fit_params[1],2)) + '/2)^2)')
        #D'-band data
        fit_data_file.writelines("\nD'-Band (Lorenzian)): y = height / ((x - center)^2 + (width / 2)^2) = "+ str(round(fit_params[10],2)) + ' / ((x- ' + str(round(fit_params[7],2)) + ')^2 + (' + str(round(fit_params[2],2)) + '/2)^2)')
        #G-band data
        fit_data_file.writelines('\nG-Band (Breit-Wigner-Fano-Function): y = ' + str(round(fit_params[13], 2)) + ' + ' + str(round(fit_params[12],2)) + ' * ((1 + (' + str(round(fit_params[11],2)) + '-x)/ (' + str(round(fit_params[4],2)) + ' * ' + str(round(fit_params[3],2)) + '))^2 / (1 + ((' + str(round(fit_params[11],2)) + ' -x)/ ' + str(round(fit_params[3],2)) + ')^2)))\n\n')
        #generating dataframe for output
        fit_data_table = pd.DataFrame(fitting_region['Wavenumber'])
        fit_data_table['Mean_Spectrum'] = fitting_region['Mean_Intensity']
        fit_data_table['D"-Band'] = lorenz(fitting_region['Wavenumber'], fit_params[0], fit_params[5], fit_params[8])
        fit_data_table['D-Band'] = lorenz(fitting_region['Wavenumber'], fit_params[1], fit_params[6], fit_params[9])
        fit_data_table["D'-Band"] = lorenz(fitting_region['Wavenumber'], fit_params[2], fit_params[7], fit_params[10])
        fit_data_table['G-Band'] = BWF(fitting_region['Wavenumber'], fit_params[3], fit_params[4], fit_params[11], fit_params[12], fit_params[13])
        fit_data_table['Convolution'] = fourpeak_3L_BWF(fitting_region['Wavenumber'], fit_params[0], fit_params[1], fit_params[2], fit_params[3], fit_params[4], fit_params[5], fit_params[6], fit_params[7], fit_params[8], fit_params[9], fit_params[10], fit_params[11], fit_params[12], fit_params[13])
        #delete unneeded columns
        fit_data_table.set_index('Wavenumber', inplace = True )
        #show only significant digits
        fit_data_table = fit_data_table.round(4).reset_index(drop = False)
        #save to file
        fit_data_file.write(fit_data_table.to_string(index = False))
    return fit_params, pcov

def three_peak_fit_and_save(fitting_region):
    #least squares fitting algorithm 
    #fit_params, pcov = scipy.optimize.curve_fit(lorenz, fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], p0 = [100, maxima_numbers.iat[0,0], 2000], bounds = ([0, 1200, 1000], [200, 1400, 40000])) 
    #fit_params, pcov = scipy.optimize.curve_fit(BWF, fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], p0 = [1, 0.1, maxima_numbers.iat[1,0], 1, 0.001]) 
    start_parameters =  [peak_width_init, 
                         peak_width_init, 
                         peak_width_init, 
                         BWF_asymmetry_init, 
                         Peak_D['Wavenumber'].values[0]-150, 
                         Peak_D['Wavenumber'].values[0], 
                         Lorenz_height_init, 
                         Lorenz_height_init, 
                         Peak_G['Wavenumber'].values[0], 
                         BWF_height_init, 
                         base_BWF_init] 
    min_bounds =        [peak_width_limits[0], 
                         peak_width_limits[0],
                         peak_width_limits[0], 
                         BWF_asymmetry_limits[0],                                                           
                         D_doubleprime_center_limits[0],
                         D_center_limits[0], 
                         Lorenz_height_limits[0],    
                         Lorenz_height_limits[0],
                         G_center_limits[0],   
                         BWF_height_limits[0], 
                         base_BWF_init] 
    max_bounds =        [peak_width_limits[1], 
                         peak_width_limits[1],
                         peak_width_limits[1],
                         BWF_asymmetry_limits[1],                                
                         D_doubleprime_center_limits[1],                         
                         D_center_limits[1], 
                         Lorenz_height_limits[1],    
                         Lorenz_height_limits[1],
                         G_center_limits[1],   
                         BWF_height_limits[1],  
                         base_BWF_init + 0.0001]
    
    fit_params, pcov = scipy.optimize.curve_fit(threepeak_2L_BWF, 
                                                fitting_region['Wavenumber'], 
                                                fitting_region['Mean_Intensity'], 
                                                p0 = start_parameters, 
                                                bounds = (min_bounds, max_bounds), 
                                                method = 'trf') 
    #generating fitting graph
    plt.close()
    plt.figure(figsize=(10,8))
    #mean Raman spectrum
    plt.plot(fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], label = 'Mean Intensity')
    #D-sideband curve
    plt.plot(fitting_region['Wavenumber'], lorenz(fitting_region['Wavenumber'], fit_params[0], fit_params[4], fit_params[6]), label = 'D" fit')
    #D-band curve
    plt.plot(fitting_region['Wavenumber'], lorenz(fitting_region['Wavenumber'], fit_params[1], fit_params[5], fit_params[7]), label = 'D fit')
    #G-band curve
    plt.plot(fitting_region['Wavenumber'], BWF(fitting_region['Wavenumber'], fit_params[2], fit_params[3], fit_params[8], fit_params[9], fit_params[10]), label = 'G fit')
    #fit convolution
    plt.plot(fitting_region['Wavenumber'], threepeak_2L_BWF(fitting_region['Wavenumber'], fit_params[0], fit_params[1], fit_params[2], fit_params[3], fit_params[4], fit_params[5], fit_params[6], fit_params[7], fit_params[8], fit_params[9], fit_params[10]), label = 'Convolution Curve')
    #axis labels
    plt.xlabel(r'Wavenumber /$cm^{-1}$')
    plt.ylabel('Normalized Intensity /a. u.')
    #legend
    plt.legend(loc = 'upper left')
    #plot range
    plt.xlim([900, 1900])
    plt.ylim([-0.05, 1.2])
    #peak labeling
    plt.text(fit_params[4], 1.0, 'D"-Peak\n(Lorenz)', horizontalalignment='center')
    plt.text(fit_params[5], 1.1, 'D-Peak\n(Lorenz)', horizontalalignment='center')
    plt.text(fit_params[8], 1.1, 'G-Peak\n(BWF)', horizontalalignment='center')
    #plotsaving
    plt.savefig(samplename + 'Mean_Raman_data_fit.png', dpi = 300)
    #plot showing
    plt.show()
    #save fitting data to file
    with open(samplename + 'fit_data.txt', 'w') as fit_data_file:
        #D-band data
        fit_data_file.writelines('D"-Band(Lorenzian): y = height / ((x - center)^2 + (width / 2)^2) = '+ str(round(fit_params[6],2)) + ' / ((x- ' + str(round(fit_params[4],2)) + ')^2 + (' + str(round(fit_params[0],2)) + '/2)^2)')
        #D-band data
        fit_data_file.writelines('\nD-Band (Lorenzian)): y = height / ((x - center)^2 + (width / 2)^2) = '+ str(round(fit_params[7],2)) + ' / ((x- ' + str(round(fit_params[5],2)) + ')^2 + (' + str(round(fit_params[1],2)) + '/2)^2)')
        #G-band data
        fit_data_file.writelines('\nG-Band (Breit-Wigner-Fano-Function): y = ' + str(round(fit_params[10], 2)) + ' + ' + str(round(fit_params[9],2)) + ' * ((1 + (' + str(round(fit_params[8],2)) + '-x)/ (' + str(round(fit_params[3],2)) + ' * ' + str(round(fit_params[2],2)) + '))^2 / (1 + ((' + str(round(fit_params[8],2)) + ' -x)/ ' + str(round(fit_params[2],2)) + ')^2)))\n\n')
        #generating dataframe for output
        fit_data_table = pd.DataFrame(fitting_region['Wavenumber'])
        fit_data_table['Mean_Spectrum'] = fitting_region['Mean_Intensity']
        fit_data_table['D"-Band'] = lorenz(fitting_region['Wavenumber'], fit_params[0], fit_params[4], fit_params[6])
        fit_data_table['D-Band'] = lorenz(fitting_region['Wavenumber'], fit_params[1], fit_params[5], fit_params[7])
        fit_data_table['G-Band'] = BWF(fitting_region['Wavenumber'], fit_params[2], fit_params[3], fit_params[8], fit_params[9], fit_params[10])
        fit_data_table['Convolution'] = threepeak_2L_BWF(fitting_region['Wavenumber'], fit_params[0], fit_params[1], fit_params[2], fit_params[3], fit_params[4], fit_params[5], fit_params[6], fit_params[7], fit_params[8], fit_params[9], fit_params[10])
        #delete unneeded columns
        fit_data_table.set_index('Wavenumber', inplace = True )
        #show only significant digits
        fit_data_table = fit_data_table.round(4).reset_index(drop = False)
        #save to file
        fit_data_file.write(fit_data_table.to_string(index = False))
    return fit_params, pcov

def two_peak_fit_and_save(fitting_region):
    #least squares fitting algorithm 
    #fit_params, pcov = scipy.optimize.curve_fit(lorenz, fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], p0 = [100, maxima_numbers.iat[0,0], 2000], bounds = ([0, 1200, 1000], [200, 1400, 40000])) 
    #fit_params, pcov = scipy.optimize.curve_fit(BWF, fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], p0 = [1, 0.1, maxima_numbers.iat[1,0], 1, 0.001]) 
    start_parameters = [peak_width_init, 
                        peak_width_init, 
                        BWF_asymmetry_init, 
                        Peak_D['Wavenumber'].values[0],  
                        Lorenz_height_init, 
                        Peak_G['Wavenumber'].values[0], 
                        BWF_height_init, 
                        base_BWF_init]
    min_bounds = [peak_width_limits[0], 
                  peak_width_limits[0], 
                  BWF_asymmetry_limits[0],                                   
                  D_center_limits[0],   
                  Lorenz_height_limits[0],
                  G_center_limits[0],   
                  BWF_height_limits[0], 
                  base_BWF_init]
    max_bounds = [peak_width_limits[1], 
                  peak_width_limits[1],
                  BWF_asymmetry_limits[1],                                
                  D_center_limits[1], 
                  Lorenz_height_limits[1],
                  G_center_limits[1],   
                  BWF_height_limits[1],  
                  base_BWF_init + 0.0001]
    fit_params, pcov = scipy.optimize.curve_fit(lorenzandBWF, 
                                                fitting_region['Wavenumber'], 
                                                fitting_region['Mean_Intensity'], 
                                                p0 = start_parameters, 
                                                bounds = (min_bounds, max_bounds), 
                                                method = 'trf') 
    #generating fitting graph
    plt.close()
    plt.figure(figsize=(10,8))
    #mean Raman spectrum
    plt.plot(fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], label = 'Mean Intensity')
    #D-band curve
    plt.plot(fitting_region['Wavenumber'], lorenz(fitting_region['Wavenumber'], fit_params[0], fit_params[3], fit_params[4]), label = 'D fit')
    #G-band curve
    plt.plot(fitting_region['Wavenumber'], BWF(fitting_region['Wavenumber'], fit_params[1], fit_params[2], fit_params[5], fit_params[6], fit_params[7]), label = 'G fit')
    #fit convolution
    plt.plot(fitting_region['Wavenumber'], lorenzandBWF(fitting_region['Wavenumber'], fit_params[0], fit_params[1], fit_params[2], fit_params[3], fit_params[4], fit_params[5], fit_params[6], fit_params[7]), label = 'Convolution Curve')
    #axis labels
    plt.xlabel(r'Wavenumber /$cm^{-1}$')
    plt.ylabel('Normalized Intensity /a. u.')
    #legend
    plt.legend(loc = 'upper left')
    #plot range
    plt.xlim([900, 1900])
    plt.ylim([-0.05, 1.2])
    #peak labeling
    plt.text(fit_params[3], 1.1, 'Lorenz', horizontalalignment='center')
    plt.text(fit_params[5], 1.1, 'BWF', horizontalalignment='center')
    #plotsaving
    plt.savefig(samplename + 'Mean_Raman_data_fit.png', dpi = 300)
    #plot showing
    plt.show()
    #save fitting data to file
    with open(samplename + 'fit_data.txt', 'w') as fit_data_file:
        #D-band data
        fit_data_file.writelines('Lorenzian Function: y = height / ((x - center)^2 + (width / 2)^2) = '+ str(round(fit_params[4],2)) + ' / ((x- ' + str(round(fit_params[3],2)) + ')^2 + (' + str(round(fit_params[0],2)) + '/2)^2)')
        #G-band data
        fit_data_file.writelines('\nBreit-Wigner-Fano-Function : y = ' + str(round(fit_params[7], 2)) + ' + ' + str(round(fit_params[6],2)) + ' * ((1 + (' + str(round(fit_params[5],2)) + '-x)/ (' + str(round(fit_params[2],2)) + ' * ' + str(round(fit_params[1],2)) + '))^2 / (1 + ((' + str(round(fit_params[5],2)) + ' -x)/ ' + str(round(fit_params[1],2)) + ')^2)))\n\n')
        #generating dataframe for output
        fit_data_table = pd.DataFrame(fitting_region['Wavenumber'])
        fit_data_table['Mean_Spectrum'] = fitting_region['Mean_Intensity']
        fit_data_table['D-Band'] = lorenz(fitting_region['Wavenumber'], fit_params[0], fit_params[3], fit_params[4])
        fit_data_table['G-Band'] = BWF(fitting_region['Wavenumber'], fit_params[1], fit_params[2], fit_params[5], fit_params[6], fit_params[7])
        fit_data_table['Convolution'] = lorenzandBWF(fitting_region['Wavenumber'], fit_params[0], fit_params[1], fit_params[2], fit_params[3], fit_params[4], fit_params[5], fit_params[6], fit_params[7])
        #delete unneeded columns
        fit_data_table.set_index('Wavenumber', inplace = True )
        #show only significant digits
        fit_data_table = fit_data_table.round(4).reset_index(drop = False)
        #save to file
        fit_data_file.write(fit_data_table.to_string(index = False))
    return fit_params, pcov


def getvalue(valueinputtext, valuetype):
    valueinput = ''
    while (type(valueinput) != int) and (type(valueinput) != float):
        print("\n'e' for exit")
        valueinput = input(valueinputtext)
        try: 
            valueinput = valuetype(valueinput)
        except:
            if valueinput == 'e':
                sys.exit()
            else:
                print('wrong input')
    return valuetype(valueinput)

def getyninput(inputquestion):
    print("\n'e' for exit")
    ok = input(inputquestion)
    while (ok != 'n') and (ok != 'y'):
        if ok == 'e':
            sys.exit()
        else:
            print('wrong input')
            ok = input(inputquestion)
    return ok


#Enter full path of data folder
path, selected = select_folder('Select the folder containing the Raman spectra as .csv files', '/')
if selected == False:
    print('No folder selected. Exiting program.')
    exit()
# Get the path to the data folder from the user
os.chdir(path)
#creates a list of file names
all_files = glob.glob(path + '/*.csv')
samplename = path.split('/')[-1]
seperator = ','
first_file = pd.read_csv(all_files[0], header = None, sep=seperator)
if len(first_file.columns) < 2:
    seperator = input('Which separator is used in the data files? ("," for comma, ";" for semicolon, " " for space): ')
# creates a big table of all data 
df = pd.concat((pd.read_csv(f, header = None, names = ['Wavenumber', 'Intensity'], sep=seperator) for f in all_files), axis = 1)
# cuts out the intensity values as table
intensities = df.pop('Intensity')
intensities.columns = range(1, int(intensities.columns.size) + 1)
#generate empty list to put in difference values later
diff_list = []
#generate plot of single spectra
for col in intensities.columns:
    plt.plot(range(len(intensities)), intensities[col]/intensities[col].max())  
    #calculate difference of each spectrum to mean spectrum and put it in to list
    diff = pd.DataFrame(intensities[col] - intensities.mean(axis=1))   
    #sum up difference spectrum as measure for possible outlier
    diff_list.extend(diff.sum())
#put out graph
plt.xticks([])
plt.yticks([])
plt.legend(range(1, len(intensities.columns) + 1))
plt.show(block=False)
#find spectrum that has highest difference to mean spectrum
max_diff_meas = diff_list.index(max(diff_list)) + 1
#give info about spectrum with highest difference to mean spectrum
print('\nThe spectrum ' + str(max_diff_meas) + ' has highest difference to mean spectrum and might be an outlier.')
graph_legend = intensities.columns.to_list()
#delete single measurements form data evaluation
delete_list = []
delete_spectra_q = getyninput('\nShould a spectrum be removed from evaluation (y/n): ') #input('Should a spectrum be removed from evaluation (y/n): ')
if delete_spectra_q == 'y' : 
    delete_mult_spectra = delete_spectra_q
    intensities.columns = graph_legend
    while delete_mult_spectra == 'y':
        delete_spectra_a = getvalue('\nWhich spectrum should be removed ' + str(graph_legend) + ': ', int) #int(input('Which spectrum should be removed from 1 - ' + str(intensities.columns.size) + ': '))
        delete_list.append(delete_spectra_a)
        if delete_spectra_a in graph_legend:
            intensities = intensities.drop(delete_spectra_a, axis=1)
            graph_legend = [element for element in graph_legend if element != delete_spectra_a]
            #empty list to put in difference values later
            diff_list = []
            spectra_list = []
            plt.close()
            #generate plot of single spectra    
            for col in intensities.columns:
                 plt.plot(range(len(intensities)), intensities[col]/intensities[col].max())  
                 #calculate difference of each spectrum to mean spectrum and put it in to list
                 diff = pd.DataFrame(intensities[col]/intensities[col].max() - intensities.mean(axis=1))   
                 #sum up difference spectrum as measure for possible outlier
                 diff_list.extend(abs(diff.sum()))
                 spectra_list.append(col)
            #put out graph
            plt.legend(graph_legend)
            #plt.xticks([])
            #plt.yticks([])
            plt.show(block=False)  
            #find spectrum that has highest difference to mean spectrum
            max_diff_index = diff_list.index(max(diff_list))
            max_diff_meas = spectra_list[max_diff_index]
            #give info about spectrum with highest difference to mean spectrum
            print('\nThe spectrum ' + str(max_diff_meas) + ' now has highest difference to mean spectrum.')
        else :
            print('\nwrong input')
        delete_mult_spectra = getyninput('Do you want to delete further spectra? (y/n): ') #input('Do you want to delete further spectra? (y/n): ')
spectra = intensities.copy(deep=True)
spectra.insert(0, 'Wavenumber', df.iloc[:,1])
spectra.to_csv(samplename + '_spectra.txt', sep = ' ')        

#suppress data frame copy operation to prevent warning
#pd.set_option('mode.chained_assignment',  None)
#normalization of Intensity of each dataset
def absolute_maximum_scale(series):
    return series / series.abs().max()
for col in intensities.columns:
    intensities[col] = absolute_maximum_scale(intensities[col])
# calculates the average intensity 
mean_intensity = intensities.mean(axis=1)
#calculates the std deviation of mean intensity
intensity_std_dev = intensities.std(axis=1)
#creates new dataframe with wavenumber as entry
mean_intensity_table = pd.DataFrame(df.iloc[:,0])
#adds mean intensity to dataframe
mean_intensity_table.insert(1, 'Mean_Intensity', mean_intensity)
#adds std dev as column to dataframe
mean_intensity_table.insert(2, 'StdDev_of_Intensity', intensity_std_dev)
plt.plot(mean_intensity_table['Wavenumber'], mean_intensity_table['Mean_Intensity'])
plt.show(block=False)

maxima_ok = 'n'    
noise = 0.05
distance = 50
iteration = 0
while maxima_ok == 'n':
    # number of data points to higher and lower X to search for local extrema
    #getvalue('Noise level (integer value) of mean spectrum (<10 for graphitized carbons, > 10 for amorphous carbons): ', int) *20 #int(input('Noise level (integer value) of mean spectrum (<10 for graphitized carbons, > 10 for amorphous carbons): ')) * 20 
    if iteration > 0:
        noise = getvalue('Minimum peak height in %:', float) /100
        distance = getvalue('Minimal distance of peaks in cm^-1:', float)
    #width = getvalue('peak width:', float)
    #slice in x in which search for maxima is done
    maxima = mean_intensity_table.iloc[20:3800]
    #search for local maxima
    maxima_pos,  maxima_int = find_peaks(maxima.Mean_Intensity.values, height = noise, distance = distance, width = (20, 500)) #
    maxima_wavenumber = []
    maxima_stddev   = []
    for entry in maxima_pos:
        maxima_wavenumber.append(maxima.iat[entry,0])
        maxima_stddev.append(maxima.iat[entry,2])
    maxima_table = pd.DataFrame(maxima_wavenumber, columns = ['Wavenumber'])
    maxima_table['Mean_Intensity'] = maxima_int['peak_heights']
    maxima_table['StdDev_of_Intensity'] = maxima_stddev
    #maxima_numbers = maxima.iloc[(argrelextrema(maxima.Mean_Intensity.values, np.greater_equal, order=n, mode='clip'))]
    #extract G and D Peak values
    Peak_2D = maxima_table.loc[(maxima_table['Wavenumber']<2750) & (maxima_table['Wavenumber']>2550)]
    Peak_G = maxima_table.loc[(maxima_table['Wavenumber']<1630) & (maxima_table['Wavenumber']>1500)]
    Peak_D = maxima_table.loc[(maxima_table['Wavenumber']<1400) & (maxima_table['Wavenumber']>1200)]
    if len(Peak_2D) == 0:
        print('\nNo 2D peak detected!\n')
    try:
        Ratio_2D_over_G = Peak_2D['Mean_Intensity'].values[0]/Peak_G['Mean_Intensity'].values[0]
        Error_2D_over_G = Peak_2D['StdDev_of_Intensity'].values[0] * 1/Peak_G['Mean_Intensity'].values[0] + Peak_G['StdDev_of_Intensity'].values[0] * (Peak_2D['Mean_Intensity'].values[0] / Peak_G['Mean_Intensity'].values[0]**2)
        print('2D Peak at ' + str(Peak_2D['Wavenumber'].values[0])+ ' cm^-1')
    except:
        Ratio_2D_over_G = None
        Error_2D_over_G = None
    try:
        Ratio_D_over_G = Peak_D['Mean_Intensity'].values[0]/Peak_G['Mean_Intensity'].values[0]
        Error_D_over_G = round(Peak_D['StdDev_of_Intensity'].values[0] * 1/Peak_G['Mean_Intensity'].values[0] + Peak_G['StdDev_of_Intensity'].values[0] * (Peak_D['Mean_Intensity'].values[0] / Peak_G['Mean_Intensity'].values[0]**2) ,3)
        print('G Peak at ' + str(Peak_G['Wavenumber'].values[0])+ ' cm^-1')
        print('D Peak at ' + str(Peak_D['Wavenumber'].values[0])+ ' cm^-1')
    except:
        Ratio_D_over_G = None
        Error_D_over_G = None

    #calculation of y-error range     
    y_lower = maxima['Mean_Intensity'] - maxima['StdDev_of_Intensity']
    y_upper = maxima['Mean_Intensity'] + maxima['StdDev_of_Intensity']  
    
    #plot creation
    plt.close()
    plt.figure(figsize=(10,8))
    #shaded standard deviation
    plt.fill_between(maxima['Wavenumber'], y_lower, y_upper, alpha=0.3, label = 'Standard Deviation')
    #mean raman spectrum
    plt.plot(maxima['Wavenumber'], maxima['Mean_Intensity'], label = 'Mean Intensity')
    #detected local maxima as scatter plot
    plt.scatter(maxima_table['Wavenumber'], maxima_table['Mean_Intensity'], label = 'Local Maxima', c='r')
    #axis labels
    plt.xlabel('Wavenumber [$cm^{-1}$]', fontsize=axislabelsize)
    plt.ylabel('Normalized Intensity [a. u.]', fontsize=axislabelsize)
    plt.xticks(fontsize=axisticksize)
    plt.yticks(fontsize=axisticksize)
    #legend
    plt.legend(fontsize=legendsize)
    #plot range
    plt.xlim([50, 3400])
    plt.ylim([-0.02, 1.05])
    #plotsaving
    plt.savefig(samplename + 'Mean_Raman_data.png', dpi = 300)
    #plot showing
    plt.show(block=False)
    iteration+=1
    maxima_ok = getyninput('Are the maxima found correct?(y/n): ') #input('Are the maxima found correct?(y/n): ')

#fitting region
fitting_region = maxima.loc[(maxima['Wavenumber']> 800) &(maxima['Wavenumber']<1900)]
#plot fitting by choice of number of peaks in the region of D and G peak
peak_number = getvalue("What numbers of Peaks you want to have fitted in the D- and G-Peak region? (2, 3, 4):", int)
if peak_number <3:
    fit_params, pcov = two_peak_fit_and_save(fitting_region)
elif peak_number >3:
    fit_params, pcov = four_peak_fit_and_save(fitting_region)
else:
    fit_params, pcov = three_peak_fit_and_save(fitting_region)

#delete unneeded columns in mean_intensity_table
mean_intensity_table.set_index('Wavenumber', inplace = True)
#show only significant digits
mean_intensity_table = mean_intensity_table.round(3)
#saves dataframe with mean intensity to file
mean_intensity_table.to_csv(samplename + 'Mean_Raman_Intensity.txt', sep = ' ')

#delete unneeded columns in maxima_numbers
maxima_table.set_index('Wavenumber', inplace = True)
#show only significant digits
maxima_table = maxima_table.round(3).reset_index(drop = False)
#save data of maximas to file
with open(samplename + 'maxima_data.txt', 'w') as maxima_data_file:
    maxima_table_string = maxima_table.to_string(index = False)
    maxima_data_file.write(maxima_table_string)
    maxima_data_file.write('\nI_D/I_G ratio is: ' + str(round(Ratio_D_over_G,3)) + ' +- ' + str(round(Error_D_over_G,3)))
    if Ratio_2D_over_G != None:
        maxima_data_file.write('\nI_2D/I_G ratio is: ' + str(round(Ratio_2D_over_G,3))+ ' +- ' + str(round(Error_2D_over_G,3)))
