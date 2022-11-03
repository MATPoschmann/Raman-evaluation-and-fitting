# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# needed python packages
import os
import pandas as pd
import glob
from scipy.signal import argrelextrema
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import sys

def getvalue(valueinputtext, valuetype):
    valueinput = ''
    while (type(valueinput) != int) and (type(valueinput) != float):
        print("'e' for exit")
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
    print("'e' for exit")
    ok = input(inputquestion)
    while (ok != 'n') and (ok != 'y'):
        if ok == 'e':
            sys.exit()
        else:
            print('wrong input')
            ok = input(inputquestion)
    return ok


#entry of filepath
os.chdir(input('Enter full path of data folder: '))

#Data folder with .csv-data of one sample to do math
path = os.getcwd()
#creates a list of file names
all_files = glob.glob(path + '/*.csv')
# creates a big table of all data 
df = pd.concat((pd.read_csv(f, header = None, names = ['Wavenumber', 'Intensity']) for f in all_files), axis = 1)
# cuts out the intensity values as table
intensities = df.pop('Intensity')
intensities.columns = range(1, int(intensities.columns.size) + 1)
#generate empty list to put in difference values later
diff_list = []
#generate plot of single spectra
for col in intensities.columns:
    plt.plot(range(len(intensities)), intensities[col])  
    #calculate difference of each spectrum to mean spectrum and put it in to list
    diff = pd.DataFrame(intensities[col] - intensities.mean(axis=1))   
    #sum up difference spectrum as measure for possible outlier
    diff_list.extend(diff.sum())
#put out graph
plt.xticks([])
plt.yticks([])
plt.legend(range(1, len(intensities.columns) + 1))
plt.show()
#find spectrum that has highest difference to mean spectrum
max_diff_meas = diff_list.index(max(diff_list)) + 1
#give info about spectrum with highest difference to mean spectrum
print('The spectrum ' + str(max_diff_meas) + ' has highest difference to mean spectrum and might be an outlier.')
graph_legend = intensities.columns.to_list()
#delete single measurements form data evaluation
delete_list = []
delete_spectra_q = getyninput('Should a spectrum be removed from evaluation (y/n): ') #input('Should a spectrum be removed from evaluation (y/n): ')
if delete_spectra_q == 'y' : 
    delete_mult_spectra = delete_spectra_q
    intensities.columns = graph_legend
    while delete_mult_spectra == 'y':
        delete_spectra_a = getvalue('Which spectrum should be removed ' + str(graph_legend) + ': ', int) #int(input('Which spectrum should be removed from 1 - ' + str(intensities.columns.size) + ': '))
        delete_list.append(delete_spectra_a)
        if delete_spectra_a in graph_legend:
            intensities = intensities.drop(delete_spectra_a, axis=1)
            graph_legend = [element for element in graph_legend if element != delete_spectra_a]
            #empty list to put in difference values later
            diff_list = []
            #generate plot of single spectra    
            for col in intensities.columns:
                 plt.plot(range(len(intensities)), intensities[col])  
                 #calculate difference of each spectrum to mean spectrum and put it in to list
                 diff = pd.DataFrame(intensities[col] - intensities.mean(axis=1))   
                 #sum up difference spectrum as measure for possible outlier
                 diff_list.extend(abs(diff.sum()))
            #put out graph
            plt.legend(graph_legend)
            #plt.xticks([])
            #plt.yticks([])
            plt.show()  
            #find spectrum that has highest difference to mean spectrum
            max_diff_meas = diff_list.index(max(diff_list)) + 1
            #give info about spectrum with highest difference to mean spectrum
            print('The spectrum ' + str(max_diff_meas) + ' now has highest difference to mean spectrum.')
        else :
            print('wrong input')
        delete_mult_spectra = getyninput('Do you want to delete further spectra? (y/n): ') #input('Do you want to delete further spectra? (y/n): ')
        

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
maxima_ok = 'n'
while maxima_ok == 'n':
    # number of data points to higher and lower X to search for local extrema
    #getvalue('Noise level (integer value) of mean spectrum (<10 for graphitized carbons, > 10 for amorphous carbons): ', int) *20 #int(input('Noise level (integer value) of mean spectrum (<10 for graphitized carbons, > 10 for amorphous carbons): ')) * 20 
    noise = getvalue('Minimum peak height in %:', float) /100
    distance = getvalue('Minimal distance of peaks in cm^-1:', float)
    #width = getvalue('peak width:', float)
    #slice in x in which search for maxima is done
    maxima = mean_intensity_table.iloc[500:3800]
    #search for local maxima
    maxima_pos,  maxima_int = find_peaks(maxima.Mean_Intensity.values, height = noise, distance = distance, width = (20, 200))
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
    G_D_Peak = maxima_table.loc[maxima_table['Mean_Intensity'] > 0.5] 
    #calculation of G over D Intensity ratio
    IGID = round(G_D_Peak.Mean_Intensity.iloc[1] / G_D_Peak.Mean_Intensity.iloc[0], 3)
    #calculation of Intensity ratio error
    IGID_Error = round(G_D_Peak.StdDev_of_Intensity.iloc[1] * 1/G_D_Peak.Mean_Intensity.iloc[0] + G_D_Peak.StdDev_of_Intensity.iloc[0] * (G_D_Peak.Mean_Intensity.iloc[1] / (G_D_Peak.Mean_Intensity.iloc[0])**2) ,3)

    #calculation of y-error range     
    y_lower = maxima['Mean_Intensity'] - maxima['StdDev_of_Intensity']
    y_upper = maxima['Mean_Intensity'] + maxima['StdDev_of_Intensity']  

    #plot creation
    plt.figure(figsize=(10,8))
    #shaded standard deviation
    plt.fill_between(maxima['Wavenumber'], y_lower, y_upper, alpha=0.3, label = 'Standard Deviation')
    #mean raman spectrum
    plt.plot(maxima['Wavenumber'], maxima['Mean_Intensity'], label = 'Mean Intensity')
    #detected local maxima as scatter plot
    plt.scatter(maxima_table['Wavenumber'], maxima_table['Mean_Intensity'], label = 'Local Maxima', c='r')
    #axis labels
    plt.xlabel('Wavenumber /$cm^{-1}$')
    plt.ylabel('Normalized Intensity /a. u.')
    #legend
    plt.legend()
    #plot range
    plt.xlim([600, 3400])
    plt.ylim([-0.02, 1.05])
    #plotsaving
    plt.savefig('Mean_Raman_data.png', dpi = 300)
    #plot showing
    plt.show()
    maxima_ok = getyninput('Are the maxima found correct?(y/n): ') #input('Are the maxima found correct?(y/n): ')


#fitting region
fitting_region = maxima.loc[831:1868,:]
#plot fitting
#fitting functions
def lorenzandBWF(x, width, widthBWF, asymBWF, center, height, centerBWF, heightBWF, baseBWF = 0):
    return (height*  (1/ ((x - center)**2 + (width / 2)**2))) + baseBWF + (heightBWF * ((1 + (centerBWF-x)/ (asymBWF * widthBWF))**2 / (1 + ((centerBWF-x)/widthBWF)**2)))
def lorenz(x, width, center, height):
    return height * (1/ ((x - center)**2 + (width / 2)**2))
def BWF(x, widthBWF, asymBWF, centerBWF , heightBWF, baseBWF):
   return baseBWF + heightBWF * ((1 + (centerBWF-x)/ (asymBWF * widthBWF))**2 / (1 + ((centerBWF-x)/widthBWF)**2))

#least squares fitting algorithm 
#fit_params, pcov = scipy.optimize.curve_fit(lorenz, fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], p0 = [100, maxima_numbers.iat[0,0], 2000], bounds = ([0, 1200, 1000], [200, 1400, 40000])) 
#fit_params, pcov = scipy.optimize.curve_fit(BWF, fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], p0 = [1, 0.1, maxima_numbers.iat[1,0], 1, 0.001]) 
fit_params, pcov = scipy.optimize.curve_fit(lorenzandBWF, fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], p0 = [2, 3, 0.02, G_D_Peak.iat[0,0], 0.02, G_D_Peak.iat[1,0], 0.1, 0], bounds = ([1, 2, 0, 1200, 0, 1550, 0, -0.2], [10000, 100000, 10, 1400, 100000, 1610, 1000, 0.001])) 
#generating fitting graph
plt.figure(figsize=(10,8))
#mean Raman spectrum
plt.plot(fitting_region['Wavenumber'], fitting_region['Mean_Intensity'], label = 'Mean Intensity')
#D-band curve
plt.plot(fitting_region['Wavenumber'], lorenz(fitting_region['Wavenumber'], fit_params[0], fit_params[3], fit_params[4]), label = 'D-Band fit')
#G-band curve
plt.plot(fitting_region['Wavenumber'], BWF(fitting_region['Wavenumber'], fit_params[1], fit_params[2], fit_params[5], fit_params[6], fit_params[7]), label = 'G-Band fit')
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
plt.savefig('Mean_Raman_data_fit.png', dpi = 300)
#plot showing
plt.show()

#save fitting data to file
with open('fit_data.txt', 'w') as fit_data_file:
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
    fit_data_table = fit_data_table.round(4)
    #save to file
    fit_data_file.write(fit_data_table.to_string())

#delete unneeded columns in mean_intensity_table
mean_intensity_table.set_index('Wavenumber', inplace = True)
#show only significant digits
mean_intensity_table = mean_intensity_table.round(3)
#saves dataframe with mean intensity to file
mean_intensity_table.to_csv('Mean_Raman_Intensity.txt')

#delete unneeded columns in maxima_numbers
maxima_table.set_index('Wavenumber', inplace = True)
#show only significant digits
maxima_table = maxima_table.round(3)
#save data of maximas to file
with open('maxima_data.txt', 'w') as maxima_data_file:
    maxima_data_file.write('I_G/I_D ratio is :' + str(IGID) + '+-' + str(IGID_Error) + '\n\n')
    maxima_table_string = maxima_table.to_string()
    maxima_data_file.write(maxima_table_string)
