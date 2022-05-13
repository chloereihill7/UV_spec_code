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
from ipywidgets import interactive
from IPython.display import clear_output
import ipywidgets as widgets
from IPython.display import display


"""function for loading in oceanview files"""
def read_oceanview(file_location):
    file_whole = pd.read_csv(file_location, header=None)
    if str(file_whole.iloc[12]) == '0    >>>>>Begin Spectral Data<<<<<\nName: 12, dtype: object':
        file = pd.read_csv(file_location, skiprows = 14, sep = '\t', header = None)
    else: 
        file = pd.read_csv(file_location, header=None)
    return file 


"""function for loading in oceanview files in times series"""
def read_oceanview_time_series(file_location, time_units = 'seconds'):
    file_whole = pd.read_csv(file_location, header=None)
    if str(file_whole.iloc[12]) == '0    >>>>>Begin Spectral Data<<<<<\nName: 12, dtype: object':
        file = pd.read_csv(file_location, skiprows = 14, sep = '\t', header = 0)
    else: 
        file = pd.read_csv(file_location, header=None)
        
    
    Timestamp = file.iloc[:,0]
        
    Time_units = []
    Time_set_to_zero = []
        
    for k in Timestamp:
        hours = float(k[11:13])
        mins = float(k[14:16])
        secs = float(k[17:])
        total_time = hours*60*60 + mins*60 +secs
        
        if time_units == 'seconds' or time_units == 'sec':
            Time_units.append(total_time)
        elif time_units == 'minutes' or time_units == 'min':
            Time_units.append(total_time/60)
        elif time_units == 'hours':
            Time_units.append(total_time/(60*60))
        else: 
            print('Please enter time units in form: minutes, min, hours, seconds, sec')
            
    for p in Time_units: 
        new_time = p - Time_units[0]
        Time_set_to_zero.append(new_time)
        
    Time = pd.DataFrame({'Timestamp':Timestamp, 'Time in ' + str(time_units): Time_units, 'Time in ' + str(time_units) +' from start': Time_set_to_zero})
    new_df = pd.concat([Time, file.iloc[:,2:]], axis=1)
    return(new_df)
           
    
"""function for loading in FLOWVPX files"""
def read_FLOWVPX(file_location, extinction_coefficient, time_units = 'seconds'):
        file = pd.read_excel(file_location)
        file['Concentration']=file.iloc[:,2]/(file.iloc[:,1]/extinction_coefficient)
        mean_conc = file.groupby('Cycle Count').mean().iloc[:,3]
        stdev_conc = file.groupby('Cycle Count').std().iloc[:,3]
        
        Time_stamps = []
        Date_stamps = []
        Time_units = []
        Time_set_to_zero = []
        
        for k in range(1, len(mean_conc) +1):
            time = file.iloc[((k)*5-1),5][12:]
            Time_stamps.append(time)
            date =  file.iloc[((k)*5-1),5][:12]
            Date_stamps.append(date)
       
        for i in range(len(Time_stamps)):
            hours_in_secs = (int(Time_stamps[i][:2]))*60*60
            mins_in_secs = (int(Time_stamps[i][3:5]))*60
            secs = int(Time_stamps[i][6:])
            total_time = hours_in_secs + mins_in_secs + secs
            
            if time_units == 'seconds' or time_units == 'sec':
                Time_units.append(total_time)
            elif time_units == 'minutes' or time_units == 'min':
                Time_units.append(total_time/60)
            elif time_units == 'hours':
                Time_units.append(total_time/(60*60))
            else: 
                print('Please enter time units in form: minutes, min, hours, seconds, sec')
        
        for p in range(len(Time_units)):
            time = (Time_units[p] - Time_units[0])
            Time_set_to_zero.append(time)
             
        new_dataframe = pd.DataFrame({'Date':Date_stamps, 'Time_stamps':Time_stamps, 'Time in ' + str(time_units): Time_units, 'Time in ' + str(time_units) +' from start': Time_set_to_zero, 'Concentration': mean_conc, 'Standard deviation': stdev_conc })
        return(new_dataframe)
    
 
"""function for loading in SYMPHONX SKID Files"""   
def read_SYMPHONX_skid(file_location, extinction_coefficient, METHOD_GUID_dict, time_units = 'seconds'): 
        file = pd.read_csv(file_location)
        METHOD_GUID = file.iloc[1,4]
        limit = int(np.argwhere(np.isnan(list(file.iloc[:,8])))[0])
        Time_stamp_raw = file.iloc[:limit,7]
        Cond_1 = file.iloc[:limit, 8]
        Cond_2 = file.iloc[:limit, 9]
        Flowrate_1 = file.iloc[:limit, 10]
        Flowrate_2 = file.iloc[:limit, 11]
        Pressure_1 = file.iloc[:limit, 12]
        Pressure_2 = file.iloc[:limit, 13]
        Temp_1 = file.iloc[:limit, 14]
        Temp_2 = file.iloc[:limit, 15]
        Abs_1 = file.iloc[:limit, 16]
        Abs_2 = file.iloc[:limit, 17]
        Abs_3 = file.iloc[:limit, 18]
        pH = file.iloc[:limit, 19]
        Vol = file.iloc[:limit, 21]
        Valve = file.iloc[:limit, 22]
        Time_units = []
        Time_set_to_zero = []
        
        for k in range(len(Time_stamp_raw)): 
            time = list(Time_stamp_raw)[k][11:]
            hours = int(time[:2])
            minutes = int(time[3:5])
            seconds = int(time[6:8])
            total_time = hours*60*60 + minutes*60 + seconds -14.5*60
             
            if time_units == 'seconds' or time_units == 'sec':
                 Time_units.append(total_time)
            elif time_units == 'minutes' or time_units == 'min':
                 Time_units.append(total_time/60)
            elif time_units == 'hours':
                 Time_units.append(total_time/(60*60))
            else: 
                 print('Please enter time units in form: minutes, min, hours, seconds, sec')   
        
        for p in range(len(Time_units)):
            time = (Time_units[p] - Time_units[0])
            Time_set_to_zero.append(time)
            
        new_dataframe = pd.DataFrame({'Time_stamp':Time_stamp_raw, 'Time in ' + str(time_units):Time_units,
                                      'Time in ' + str(time_units) +' from start': Time_set_to_zero, 'Conductivity 1':Cond_1,
                                      'Conductivity 2':Cond_2, 'Flowrate 1':Flowrate_1, 'Flowrate 2':Flowrate_2, 'Pressure 1':Pressure_1,
                                      'Pressure 2':Pressure_2, 'Temp 1':Temp_1, 'Temp 2':Temp_2, 'Abs 1': Abs_1, 'Abs 2': Abs_2, 'Abs 3':Abs_3, 
                                      'pH':pH, 'Volume':Vol,  'Valve':Valve})
        
        Type = file.iloc[limit+1:,3]
        Time_stamp = file.iloc[limit+1:,2]
        Volume = file.iloc[limit+1:,4]
        Stage_name = file.iloc[limit+1:,5]
        
        Time_set_zero_2 = []
        
        for k in range(len(Time_stamp)): 
            time = list(Time_stamp)[k][11:]
            hours = int(time[:2])
            minutes = int(time[3:5])
            seconds = int(time[6:8])
            total_time = hours*60*60 + minutes*60 + seconds -14.5*60
             
            if time_units == 'seconds' or time_units == 'sec':
                 Time_set_zero_2.append(total_time)
            elif time_units == 'minutes' or time_units == 'min':
                 Time_set_zero_2.append(total_time/60)
            elif time_units == 'hours':
                 Time_set_zero_2.append(total_time/(60*60))
            else: 
                 print('Please enter time units in form: minutes, min, hours, seconds, sec') 
        
        stage_dataframe = pd.DataFrame({'Stage number': Stage_name, 'Time stamp':Time_stamp,  'Time in ' + str(time_units):Time_set_zero_2, 'Volume':Volume, 'Description':Type})
        return(new_dataframe, stage_dataframe, print('RUN: ' + str(METHOD_GUID_dict[METHOD_GUID])), METHOD_GUID_dict[METHOD_GUID], METHOD_GUID)
    
        


"""function for plotting single data set on single axis""" 
def single_dataset_plot(real_x_axis, real_y_axis, colours, labels, lower_xlim = 240, higher_xlim = 310,
                        lower_ylim = 0, higher_ylim = 1.35, figsize_x = 10, figsize_y = 7, 
                        x_label = 'Wavelength (nm)', y_label = 'Absorbance (A.U.)',
                        Title = 'Real and Predicted Absorbance', Fontsize = 16, grid = 'on', legend = 'yes', leg_fontsize = 12): 
    if len(real_x_axis) < len(real_y_axis) and type(real_x_axis) == 'list': 
        print('Please input list of wavelength data for each absorbance dataset')
    elif type(real_x_axis)!= 'list': 
        print('Please input a list of wavelength data and absorbance data')
    else: 
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
    
    """Find the indexes that correlate with around ~250 and 255 nm for the predicted and measured absorbance spectrum.
    We want to cut off the spectra so that we search for the secondary absorbance maxima"""
    for point in range(len(WL_pred)): 
        if int(WL_pred[point]) == 250: 
            chosen_point_pred.append(point)
    for point in range(len(WL_measured)): 
        if int(WL_measured[point]) == 255: 
            chosen_point_measured.append(point)
            
    """Restrict selected wavelength range so it only includes the secondary absorbance maxima"""
    Measured_range = Abs_measured[chosen_point_measured[0]:]
    Pred_range =Abs_pred[chosen_point_pred[0]:]
    
    
    """Search for the index of the maximum in each absorbance data set"""
    max_pred_index = np.where(Pred_range == np.nanmax(Pred_range))[0][0]
    max_meas_index = Measured_range.index(max(Measured_range))
            
    """Find wavelength difference between maxima of predicted and measured spectra"""
    WL_difference = WL_measured[max_meas_index + chosen_point_measured[0]] - WL_pred[max_pred_index + chosen_point_pred[0]]
    return WL_difference

"""Find the index of the absorbance maxima at 280 nm in each absorbance spectrum
Return the absorbance at 280 nm for each sample """
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

"""Reset the index of python file so that files have regular intervals of x-axis (Wavelength)"""
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

"label columns of a dataframe"
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
                

        
"Convert a list containing dataframes to a dictionary "
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
        
"Convert a dictionary to a pandas dataframe" 
def convert_dict_to_df(input_data): 
    df = pd.concat([pd.Series(v, name=k) for k, v in input_data.items()], axis=1)
    return df

def prepare_for_horizontal_fit(input_predicted, input_measured):
    Add_to_measured = len(input_predicted)#append to start and end of measured data
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
         ax.legend(fontsize = 10, loc='upper right', frameon=True)
         
         
def old_horizontal_curve_fit(predicted_spectrum, measured_spectra, wavelengths):
    predicted_spectrum = predicted_spectrum['Predicted_spectrum']/1000
    predicted_spectrum = pd.DataFrame(predicted_spectrum, columns = ['Predicted_spectrum'])
    Length = int(len(predicted_spectrum)/3)
    wavelengths = np.arange(198, 320, 0.2)
    upper_limit = Length + 15/round((wavelengths[1]-wavelengths[0]),2)
    Total = Length*2
    RMSE = []
    WL_shift = []

    
    for k in range(len(measured_spectra)): 
        Measured_spectrum = measured_spectra[k]
        rmse_tot = []
        
        for i in np.arange(Length, upper_limit): 
            f = int(i)
            zeros_pred_1 = pd.DataFrame([0]*f, columns = ['Predicted_spectrum'])
            zeros_pred_2 = pd.DataFrame([0]*(Total-f), columns = ['Predicted_spectrum'])
            #predicted_abs = pd.DataFrame({'Predicted Spectrum':predicted_spectrum})
            New_pred_spectra = pd.concat([zeros_pred_1, predicted_spectrum, zeros_pred_2], ignore_index = True)
            SUM = 0
            count = 0
            
            for j in np.arange(915, 1073): 
                diff = Measured_spectrum['Real_spectrum'][j] - New_pred_spectra['Predicted_spectrum'][j]
                diff_sqrd = diff**2
                count +=1
                SUM += diff_sqrd
                MSS = SUM/count
                rmse = math.sqrt(MSS)
            rmse_tot.append(rmse)
            #print(rmse_tot)
            #print(rmse_tot)
        minimum = rmse_tot.index(min(rmse_tot))
        wavelength_shift = round(minimum*0.2,2)
        RMSE.append(minimum)
        WL_shift.append(wavelength_shift)
        
    return [RMSE,  WL_shift]
    
def extinction_dataframe(Data):
    x = Data['x']
    F = Data['F']
    W = Data['W']
    Y = Data['Y']
    CY = Data['CY'] #Cystine
    #C = Data['C']
    ext_dict = {'Phenylalanine':F, 'Tryptophan':W, 'Tyrosine':Y, 'Cystine':CY}
    
    df_list = []
    for k in range(len(ext_dict)): 
        df1 = pd.DataFrame(list(ext_dict.values())[k])
        df1.columns = [list(ext_dict.keys())[k]]
        df_list.append(df1)
    

    ext_df = pd.concat(df_list, axis=1)
    ext_df
    ext_dict = {'F':F, 'W':W, 'Y':Y, 'CY':CY}
    return [x, ext_df, ext_dict]

def number_of_residues_dict(AA_count):
    AA_count_new =AA_count.set_index('Amino Acid')
    C_n = float(AA_count_new.loc['Cysteine'])
    CY_n = float(AA_count_new.loc['Cystine'])
    W_n = float(AA_count_new.loc['Tryptophan'])
    Y_n = float(AA_count_new.loc['Tyrosine'])
    F_n = float(AA_count_new.loc['Phenylalanine'])
    MW = float(AA_count_new.loc['Protein MW'])
    ext_number_dict = {'C_n':C_n, 'CY_n':CY_n, 'W_n':W_n, 'Y_n':Y_n, 'F_n':F_n, 'MW':MW}
    
   
    return ext_number_dict


def widget_plot(wavelength, dataframe_labeled_extinction): 
    #%matplotlib inline


    #jtplot.style()
    df = dataframe_labeled_extinction
    x = wavelength
    
    opts = df.columns.values

    selector = widgets.SelectMultiple(
    options=opts,
    value=[opts[1]],
    rows=len(opts),
    description='Variables',
    disabled=False)
    
    
    output = widgets.Output()
    
    display(selector)
    display(output)
    
    def multiplot(widg):
        choices = widg['new']
        data = df.loc[:, choices] if choices else df
        output.clear_output(wait=True)
        with output: 
            if len(choices) > 1:
                plt.figure(figsize = (12, 8))
                plt.plot(x, data, label = choices, linewidth = 2)
            else: 
                plt.figure(figsize = (12, 8))
                plt.plot(x, data, label = data.columns[0], linewidth = 2)
            plt.xlabel('Wavelength (nm)', fontsize = 14)
            plt.ylabel('Molar Extinction Coefficient (M$^-$$^1$cm$^-$$^1$)', fontsize = 14)
            plt.legend(fontsize = 12)
            plt.xlim(190, 360)
            plt.grid()
            plt.show()
            
    
    selector.observe(multiplot, names='value')
    
    
def calc_molar_ext(Data, AA_count): 
   
    F = Data.iloc[:, 1]
    W = Data.iloc[:, 2]
    Y = Data.iloc[:, 3]
    CY = Data.iloc[:, 4] #Cystine
    
    
    AA_count_new =AA_count.set_index('Amino Acid')
    CY_n = float(AA_count_new.loc['Cystine'])
    W_n = float(AA_count_new.loc['Tryptophan'])
    Y_n = float(AA_count_new.loc['Tyrosine'])
    F_n = float(AA_count_new.loc['Phenylalanine'])
    
    NET_Mol_Ext = W*W_n + Y*Y_n + CY*CY_n + F*F_n #net molar extinction    
    return NET_Mol_Ext  

def calc_mass_ext(Data, AA_count): 
     F = Data.iloc[:, 1]
     W = Data.iloc[:, 2]
     Y = Data.iloc[:, 3]
     CY = Data.iloc[:, 4] #Cystine
     
     
     AA_count_new =AA_count.set_index('Amino Acid')
     CY_n = float(AA_count_new.loc['Cystine'])
     W_n = float(AA_count_new.loc['Tryptophan'])
     Y_n = float(AA_count_new.loc['Tyrosine'])
     F_n = float(AA_count_new.loc['Phenylalanine'])
     MW = float(AA_count_new.loc['Protein MW'])
     
     NET_Mol_Ext = W*W_n + Y*Y_n + CY*CY_n + F*F_n #net molar extinction
     NET_Mass_Ext = NET_Mol_Ext/MW
     return NET_Mass_Ext    
            
def read_FLOWVPX(file_location, extinction_coefficient, time_units = 'seconds'):
        file = pd.read_excel(file_location)
        file['Concentration']=file.iloc[:,2]/(file.iloc[:,1]/extinction_coefficient)
        mean_conc = file.groupby('Cycle Count').mean().iloc[:,3]
        stdev_conc = file.groupby('Cycle Count').std().iloc[:,3]
        
        Time_stamps = []
        Date_stamps = []
        Time_units = []
        Time_set_to_zero = []
        
        for k in range(1, len(mean_conc) +1):
            time = file.iloc[((k)*5-1),5][12:]
            Time_stamps.append(time)
            date =  file.iloc[((k)*5-1),5][:12]
            Date_stamps.append(date)
       
        for i in range(len(Time_stamps)):
            hours_in_secs = (int(Time_stamps[i][:2]))*60*60
            mins_in_secs = (int(Time_stamps[i][3:5]))*60
            secs = int(Time_stamps[i][6:])
            total_time = hours_in_secs + mins_in_secs + secs
            
            if time_units == 'seconds' or time_units == 'sec':
                Time_units.append(total_time)
            elif time_units == 'minutes' or time_units == 'min':
                Time_units.append(total_time/60)
            elif time_units == 'hours':
                Time_units.append(total_time/(60*60))
            else: 
                print('Please enter time units in form: minutes, min, hours, seconds, sec')
        
        for p in range(len(Time_units)):
            time = (Time_units[p] - Time_units[0])
            Time_set_to_zero.append(time)
             
        new_dataframe = pd.DataFrame({'Date':Date_stamps, 'Time_stamps':Time_stamps, 'Time in ' + str(time_units): Time_units, 'Time in ' + str(time_units) +' from start': Time_set_to_zero, 'Concentration': mean_conc, 'Standard deviation': stdev_conc })
        return(new_dataframe)
    
 
    
def read_SYMPHONX_skid(file_location, extinction_coefficient, METHOD_GUID_dict, time_units = 'seconds'): 
        file = pd.read_csv(file_location)
        METHOD_GUID = file.iloc[1,4]
        limit = int(np.argwhere(np.isnan(list(file.iloc[:,8])))[0])
        Time_stamp_raw = file.iloc[:limit,7]
        Cond_1 = file.iloc[:limit, 8]
        Cond_2 = file.iloc[:limit, 9]
        Flowrate_1 = file.iloc[:limit, 10]
        Flowrate_2 = file.iloc[:limit, 11]
        Pressure_1 = file.iloc[:limit, 12]
        Pressure_2 = file.iloc[:limit, 13]
        Temp_1 = file.iloc[:limit, 14]
        Temp_2 = file.iloc[:limit, 15]
        Abs_1 = file.iloc[:limit, 16]
        Abs_2 = file.iloc[:limit, 17]
        Abs_3 = file.iloc[:limit, 18]
        pH = file.iloc[:limit, 19]
        Vol = file.iloc[:limit, 21]
        Valve = file.iloc[:limit, 22]
        Time_units = []
        Time_set_to_zero = []
        
        for k in range(len(Time_stamp_raw)): 
            time = list(Time_stamp_raw)[k][11:]
            hours = int(time[:2])
            minutes = int(time[3:5])
            seconds = int(time[6:8])
            total_time = hours*60*60 + minutes*60 + seconds -14.5*60
             
            if time_units == 'seconds' or time_units == 'sec':
                 Time_units.append(total_time)
            elif time_units == 'minutes' or time_units == 'min':
                 Time_units.append(total_time/60)
            elif time_units == 'hours':
                 Time_units.append(total_time/(60*60))
            else: 
                 print('Please enter time units in form: minutes, min, hours, seconds, sec')   
        
        for p in range(len(Time_units)):
            time = (Time_units[p] - Time_units[0])
            Time_set_to_zero.append(time)
            
        new_dataframe = pd.DataFrame({'Time_stamp':Time_stamp_raw, 'Time in ' + str(time_units):Time_units,
                                      'Time in ' + str(time_units) +' from start': Time_set_to_zero, 'Conductivity 1':Cond_1,
                                      'Conductivity 2':Cond_2, 'Flowrate 1':Flowrate_1, 'Flowrate 2':Flowrate_2, 'Pressure 1':Pressure_1,
                                      'Pressure 2':Pressure_2, 'Temp 1':Temp_1, 'Temp 2':Temp_2, 'Abs 1': Abs_1, 'Abs 2': Abs_2, 'Abs 3':Abs_3, 
                                      'pH':pH, 'Volume':Vol,  'Valve':Valve})
        
        Type = file.iloc[limit+1:,3]
        Time_stamp = file.iloc[limit+1:,2]
        Volume = file.iloc[limit+1:,4]
        Stage_name = file.iloc[limit+1:,5]
        
        Time_set_zero_2 = []
        
        for k in range(len(Time_stamp)): 
            time = list(Time_stamp)[k][11:]
            hours = int(time[:2])
            minutes = int(time[3:5])
            seconds = int(time[6:8])
            total_time = hours*60*60 + minutes*60 + seconds -14.5*60
             
            if time_units == 'seconds' or time_units == 'sec':
                 Time_set_zero_2.append(total_time)
            elif time_units == 'minutes' or time_units == 'min':
                 Time_set_zero_2.append(total_time/60)
            elif time_units == 'hours':
                 Time_set_zero_2.append(total_time/(60*60))
            else: 
                 print('Please enter time units in form: minutes, min, hours, seconds, sec') 
        
        stage_dataframe = pd.DataFrame({'Stage number': Stage_name, 'Time stamp':Time_stamp,  'Time in ' + str(time_units):Time_set_zero_2, 'Volume':Volume, 'Description':Type})
        return(new_dataframe, stage_dataframe, print('RUN: ' + str(METHOD_GUID_dict[METHOD_GUID])), METHOD_GUID_dict[METHOD_GUID], METHOD_GUID)
    
        
