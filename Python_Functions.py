# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 12:57:14 2022

@author: mjki4cr2
"""
"""include all regularly used packages to import into python file"""
import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np
from scipy.signal import savgol_filter
import math
from pandas import DataFrame

"""function for plotting single data set on single axis""" 
def single_dataset_plot(real_x_axis, real_y_axis, colours, labels, lower_xlim = 240, higher_xlim = 310,
                        lower_ylim = 0, higher_ylim = 1.35, figsize_x = 10, figsize_y = 7, 
                        x_label = 'Wavelength (nm)', y_label = 'Absorbance (A.U.)',
                        Title = 'Real and Predicted Absorbance', Fontsize = 16, grid = 'on', legend = 'yes', leg_fontsize = 12): 
    plt.figure(figsize=(figsize_x, figsize_y))
    for sample in range(len(real_y_axis)): 
        plt.plot(real_x_axis[sample], real_y_axis[sample], color = colours[sample], linewidth=2, label=labels[sample])
    plt.xlabel(x_label, fontsize = Fontsize)
    plt.ylabel(y_label, fontsize = Fontsize)
    plt.title(Title, fontsize = Fontsize)
    plt.ylim(lower_ylim, higher_ylim)
    plt.xlim(lower_xlim, higher_xlim)
    if grid == 'on': 
        plt.grid()
    if legend == 'yes':
        plt.legend(fontsize = leg_fontsize, frameon = False)
    plt.show()
    
"""function for plotting real and predicted data sets on the same set of axis"""
def figure_real_vs_predicted(real_x_axis, real_y_axis, predicted_x_axis, predicted_y_axis, colours, labels, 
                             lower_xlim = 240, higher_xlim = 310, lower_ylim = 0, higher_ylim = 1.35, 
                             figsize_x = 10, figsize_y = 7, x_label = 'Wavelength (nm)', y_label = 'Absorbance (A.U.)',
                            Title = 'Real and Predicted Absorbance', Fontsize = 16, grid = 'on', legend = 'yes', leg_fontsize = 12): 
    plt.figure(figsize=(figsize_x, figsize_y))
    for sample in range(len(predicted_y_axis)): 
        plt.plot(real_x_axis[sample], real_y_axis[sample], color = colours[sample], linewidth=2, label=labels[sample])
        plt.plot(predicted_x_axis, predicted_y_axis[sample], color = colours[sample], ls='dotted', linewidth = 2, label = 'predicted ' + labels[sample])
    plt.xlabel(x_label, fontsize = Fontsize)
    plt.ylabel(y_label, fontsize = Fontsize)
    plt.title(Title, fontsize = Fontsize)
    plt.ylim(lower_ylim, higher_ylim)
    plt.xlim(lower_xlim, higher_xlim)
    if grid == 'on': 
        plt.grid()
    if legend == 'yes':
        plt.legend(fontsize = leg_fontsize, frameon = False)
    plt.show()
    
"""function for plotting real and predicted data sets on the same set of axis where each predicted spectra is shifted by a different amount"""
def figure_real_vs_predicted_shifted(real_x_axis, real_y_axis, predicted_x_axis, predicted_y_axis, colours, labels, 
                             lower_xlim = 240, higher_xlim = 310, lower_ylim = 0, higher_ylim = 1.35, 
                             figsize_x = 10, figsize_y = 7, x_label = 'Wavelength (nm)', y_label = 'Absorbance (A.U.)',
                            Title = 'Real and Predicted Absorbance', Fontsize = 16, grid = 'on', legend = 'yes', leg_fontsize = 12): 
    plt.figure(figsize=(figsize_x, figsize_y))
    for sample in range(len(predicted_y_axis)): 
        plt.plot(real_x_axis[sample], real_y_axis[sample], color = colours[sample], linewidth=2, label=labels[sample])
        plt.plot(predicted_x_axis[sample], predicted_y_axis[sample], color = colours[sample], ls='dotted', linewidth = 2, label = 'predicted ' + labels[sample])
    plt.xlabel(x_label, fontsize = Fontsize)
    plt.ylabel(y_label, fontsize = Fontsize)
    plt.title(Title, fontsize = Fontsize)
    plt.ylim(lower_ylim, higher_ylim)
    plt.xlim(lower_xlim, higher_xlim)
    if grid == 'on': 
        plt.grid()
    if legend == 'yes':
        plt.legend(fontsize = leg_fontsize, frameon = False)
    plt.show()
    
def A280_maxima_shift(Abs_measured, Abs_pred, WL_measured, WL_pred): 
    
    chosen_point_pred = []
    chosen_point_measured = []
    
    #Find the indexes that correlate with around ~250 and 255 nm for the predicted and measured absorbance spectrum.
    #We want to cut off the spectra so that we search for the secondary absorbance maxima
    for point in range(len(WL_pred)): 
        if int(WL_pred[point]) == 250: 
            chosen_point_pred.append(point)
    for point in range(len(WL_measured)): 
        if int(WL_measured[point]) == 255: 
            chosen_point_measured.append(point)
            
    #Restrict selected wavelength range so it only includes the secondary absorbance maxima
    Measured_range = Abs_measured[chosen_point_measured[0]:]
    Pred_range =Abs_pred[chosen_point_pred[0]:]
    
    
    #Search for the index of the maximum in each absorbance data set
    max_pred_index = np.where(Pred_range == np.nanmax(Pred_range))[0][0]
    max_meas_index = Measured_range.index(max(Measured_range))
            
    #Find wavelength difference between maxima of predicted and measured spectra
    WL_difference = WL_measured[max_meas_index + chosen_point_measured[0]] - WL_pred[max_pred_index + chosen_point_pred[0]]
    return WL_difference

#Find the index of the absorbance maxima at 280 nm in each absorbance spectrum
#Return the absorbance at 280 nm for each sample 
def A280(Wavelength_data, Absorbance_data): 
    A280_index = []
    A280_absorbance = []
    for spectrum in range(len(Wavelength_data)): 
        wavelength = Wavelength_data[spectrum]
        WL_difference = []
        for point in range(len(wavelength)): 
            difference = np.sqrt((280 - wavelength[point])**2)
            WL_difference.append(difference)
        A280_position = np.where(WL_difference == np.nanmin(WL_difference))[0][0]
        A280_index.append(A280_position)
        
    for spectrum in range(len(Wavelength_data)): 
        Abs = Absorbance_data[spectrum]
        idx = A280_index[spectrum]
        A280 = Abs[idx]
        A280_absorbance.append(A280)
        
    return (A280_absorbance)

#Reset the index of python file so that files have regular intervals of x-axis (Wavelength)
def index_reset(input_spectrum, x_index_label = 'x', idx = np.arange(190, 320, 0.2)):
    new_spectrum = input_spectrum.set_index(x_index_label)
    if new_spectrum.index.duplicated().any() == 'TRUE':
        new_spectrum = new_spectrum.loc[~new_spectrum.index.duplicated(),:] 
    else: 
        new_spectrum = new_spectrum
    if len(new_spectrum.columns) == 1: 
        NEW = np.interp(idx, new_spectrum.index.values , new_spectrum.values[:,0])
        New_idx_data = pd.DataFrame({'x':idx,'Abs':NEW})
    else: 
        New_idx_data = []
        for spectrum in range(len(new_spectrum.columns)):
            New_1 = np.interp(idx, new_spectrum.index.values , new_spectrum.values[:,spectrum])
            New_real_data_1 = pd.DataFrame({'x':idx, new_spectrum.columns[spectrum]:New_1})
            New_idx_data.append(New_real_data_1)
            
    return New_idx_data

#label columns of a dataframe
def label_columns_abs(input_dataframe, base_name = 'Abs'): 
    dataframe_length = len(input_dataframe.columns)-1
    if dataframe_length == 1: 
        labels_list = ['x', 'Abs']
        input_dataframe.columns = labels_list
    else: 
        labels_list = ['x']
        for k in range(len(input_dataframe.columns)-1):
            label = 'Abs_' + str(k + 1)
            labels_list.append(label)
        input_dataframe.columns = labels_list
    return input_dataframe
                

        
#Convert a list containing dataframes to a dictionary 
def convert_dflist_to_dict(input_dataframe):
    if type(input_dataframe) == list and type(input_dataframe[0]) == pd.core.frame.DataFrame:
        dataframe_dict = {input_dataframe[0].columns[0]:input_dataframe[0].iloc[:,0]}
        dataframe_dict.update({input_dataframe[i].columns[1]: input_dataframe[i].iloc[:,1] for i in range(0, len(input_dataframe))})
        return dataframe_dict 
    elif type(input_dataframe)== pd.core.frame.DataFrame:
        if len(input_dataframe.columns) == 2:
            dataframe_dict = {input_dataframe.columns[0]: input_dataframe.iloc[:,0], input_dataframe.columns[1]: input_dataframe.iloc[:,1] }
        else: 
            dataframe_dict = {input_dataframe.columns[0]: input_dataframe.iloc[:,0]}
            dataframe_dict.update({input_dataframe.columns[i]: input_dataframe.iloc[:,i] for i in range(1, len(input_dataframe.columns)-1)})
        return dataframe_dict
    else:
        print('data_input must be a dataframe or a list of dataframes')
        
#Convert a dictionary to a pandas dataframe 
def convert_dict_to_df(input_data): 
    df = pd.concat([pd.Series(v, name=k) for k, v in input_data.items()], axis=1)
    return df

def prepare_for_horizontal_fit(input_predicted, input_measured):
    Add_to_measured = len(input_predicted) #append to start and end of measured data
    Add_to_predicted = 2*Add_to_measured #append to end of predicted data
    list_of_zeros = [0]*Add_to_measured
    zeros = pd.DataFrame({'Real_spectrum': list_of_zeros})
    
    list_of_df = []
    for k in range(1, len(input_measured.columns)): 
        x =pd.DataFrame({'Real_spectrum':input_measured.iloc[:,k]})
        list_of_df.append(x)
        
    New_real_spectra = []
    for k in range(len(list_of_df)):
        New_spectrum = pd.concat([zeros, list_of_df[k], zeros], ignore_index = True)
        New_real_spectra.append(New_spectrum)
    
    list_of_zeros_pred = [0]*Add_to_predicted
    zeros_predicted = pd.DataFrame({'Predicted_spectrum': list_of_zeros_pred})
    predicted = pd.DataFrame({'Predicted_spectrum': input_predicted.iloc[:,1]})
    New_predicted = pd.concat([predicted, zeros_predicted], ignore_index = True)
    
    return [New_real_spectra, New_predicted]

def plot_separate_plots(Input_list_of_spectra, labels):
    colour_list = []
    for k in range(len(Input_list_of_spectra)): 
        colour = 'C' + str(k)
        colour_list.append(colour)
    label_list = []
    
    for k in range(len(Input_list_of_spectra)): 
        label = labels[k] + ' mgml$^-$$^1$'
        label_list.append(label)
        
    rows = math.ceil(len(Input_list_of_spectra)/3)
    
    fig, axs = plt.subplots(rows, 3, figsize=(20,8))
    for i in np.arange(rows):
        row = i
        for j in range(0,3): 
            axs[row, j].plot(Input_list_of_spectra[j + row*3], colour_list[j + row*3], label = label_list[j + row*3])
    
    for ax in axs.flat:
         ax.set(xlabel='Number of data points', ylabel='Absorbance')
         leg = ax.legend(fontsize = 10, loc='upper right', frameon=True)
    
    
        
    
        
         

