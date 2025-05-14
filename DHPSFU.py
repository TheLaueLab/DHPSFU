# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 20:43:01 2024

@author: Ziwei Zhang
"""

import os
import numpy as np
import pandas as pd
import math
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from datetime import datetime

# Paths and extension names
path = r'D:\DHPSFU\Data\Loc10' # path containing all peak-fit localisation data.
ext = r'.xls' # Extension of data file (.xls for GDSC SMLM)
calib_name = r"D:\DHPSFU\Data\Calib.xls" # calibration file name


# General parameters
px_size = 1 # nm 
precision = 100 # precision cutoff in nm
calib_step = 60 # step length of calibration in nm
poly_degree = 20; # Polynomial fitting degree for calibration
fitting_mode = 'Frame' # fitting mode, can be 'Frame', 'Angle', or 'Z'. Default is 'Frame'
range_to_fit = (1, 82) # range for fitting. Units: 'Z' mode in nm; 'Angle' mode in degrees; 'Frame' mode in number. Default is (1, 97) in frames
initial_distance_filter = (4.5, 10) # minimum and maximum distance between a pair of dots in px
frame_number = 10000
general_paras = [px_size, precision, initial_distance_filter]

# Filtering parameters
enable_filters = True # true if enable all filters
enable_filter_calib_range = True # remove localisations out of the angular range of calibration; if False, polynomial fit is extrapolated beyond the range.
enable_filter_distance = True # remove dots with unexpected distances
distance_dev = 0.2 # relative deviation of the distance between dots, compared to calibration
enable_filter_intensity_ratio = True # filter based on the ratio of intensities between the dots
intensity_dev = 1.5 # relative deviation of the intensity difference between dots, compared to calibration
filter_paras = [enable_filters, enable_filter_calib_range, 
                enable_filter_distance, distance_dev, 
                enable_filter_intensity_ratio, intensity_dev]

#%% Functions
def get_data(path, ext):
    datalist = [] 
    for root, dirs, files in os.walk(path):  
        for filespath in files: 
            if filespath[-len(ext):]==ext:
                datalist.append(os.path.join(root,filespath))
    report = pd.DataFrame(columns = ['Filename', 'Pos', 'Slice', 'Num_of_locs'])
    report.Filename = datalist
    poslist = []
    slicelist = []
    new_fname_list = []
    histo_list = []
    for fname in datalist:
        pos_Pos = fname.find('Pos')+3
        pos_Pos2 = fname.find('\\', pos_Pos)
        pos = fname[pos_Pos:pos_Pos2]
        pos_slice = fname.find('slice')+5
        pos_slice2 = fname.find('.', pos_slice)
        slice_n = fname[pos_slice:pos_slice2]
        new_fname = fname.replace(ext, '_processed.csv')
        new_histo = fname.replace(ext, '_histogram.png')
        poslist.append(pos)
        slicelist.append(slice_n)
        new_fname_list.append(new_fname)
        histo_list.append(new_histo)
    report.Pos = poslist
    report.Slice = slicelist
    return report, new_fname_list, histo_list

# Calibration function
def DH_calibration(calib_name, calib_step, initial_distance_filter, poly_degree, fitting_mode = 'Frame', range_to_fit = (1, 97)):
    cols = ['#Frame', 'X (px)', 'Y (px)', 'Intensity (photon)'] # useful column names
    calib_raw_data = pd.read_csv(calib_name, skiprows=8, sep = '\t') # read xls file
    calib_data = calib_raw_data[cols] # take useful data
    frame_num = np.max(calib_data['#Frame']) # take frame number
    frame_list = list(calib_data['#Frame'])
    
    # remove bad frames (with non-2 localisations)
    bad_frames = []
    for i in range(1, frame_num+1):
        if frame_list.count(i) == 1 | frame_list.count(i) > 2 :
            bad_frames.append(i)
    if bad_frames != []:
        for j in range(0, len(bad_frames)):
            calib_data = calib_data.drop(calib_data[calib_data['#Frame']==bad_frames[j]].index)
    
    # take data from dataframe for analysis
    frames = calib_data['#Frame']
    true_frame_list = list(frames.drop_duplicates())
    xs = calib_data['X (px)']
    ys = calib_data['Y (px)']
    ins = calib_data['Intensity (photon)']
    
    true_frames = np.asarray(true_frame_list)
    avg_x_all = np.zeros(len(true_frame_list))
    avg_y_all = avg_x_all.copy()
    avg_dist_all = avg_x_all.copy()
    avg_inten_all = avg_x_all.copy()
    avg_ratio_all = avg_x_all.copy()
    angle_all = avg_x_all.copy()
    
    # Calculate distance, relative intensity, average x, y, average intensity, and angles
    for i in range(len(true_frame_list)):
        frame = true_frame_list[i]
        index_frame = list(frames[frames==frame].index)
        fm1 = int(index_frame[0])
        fm2 = int(index_frame[1])
        xs1 = xs[fm1]
        xs2 = xs[fm2]
        ys1 = ys[fm1]
        ys2 = ys[fm2]
        avg_x = (xs1+xs2)/2
        avg_y = (ys1+ys2)/2
        dist_xy = np.sqrt((xs2-xs1)**2+(ys2-ys1)**2)
        avg_intensity = (ins[fm1] + ins[fm2])/2
        ratio_intensity = np.max([ins[fm1], ins[fm2]])/np.min([ins[fm1], ins[fm2]])
        angle = math.atan2(ys2-ys1, xs2-xs1)
        if angle<0:
            angle = angle+np.pi
        avg_x_all[i] = avg_x
        avg_y_all[i] = avg_y
        avg_dist_all[i] = dist_xy
        avg_inten_all[i] = avg_intensity
        avg_ratio_all[i] = ratio_intensity
        angle_all[i] = angle
    
    # Combine all data
    raw = pd.DataFrame({'Frame': true_frames, 'avg_x': avg_x_all, 
                        'avg_y': avg_y_all, 'avg_distance': avg_dist_all,
                        'avg_intensity': avg_inten_all,
                        'ratio': avg_ratio_all,
                        'Angle': angle_all})
    
    # Filter data, remove pairs with too short/long distances
    final = raw.drop(raw[raw['avg_distance']>initial_distance_filter[1]].index)
    final = final.drop(final[final['avg_distance']<initial_distance_filter[0]].index)
    
    # Filter by mode
    if fitting_mode=='Frame':
        final = final.drop(final[final['Frame']<range_to_fit[0]].index)
        final = final.drop(final[final['Frame']>range_to_fit[1]].index)
    elif fitting_mode == 'Angle':
        final = final.drop(final[final['Angle']<range_to_fit[0]].index)
        final = final.drop(final[final['Angle']>range_to_fit[1]].index)
    elif fitting_mode == 'Z':
        final = final.drop(final[final['Frame']*calib_step<range_to_fit[0]].index)
        final = final.drop(final[final['Frame']*calib_step>range_to_fit[1]].index)
    
    # Fitting
    dx = np.polyfit(final['Frame']*calib_step, final['avg_x'], poly_degree)
    dy = np.polyfit(final['Frame']*calib_step, final['avg_y'], poly_degree)
    dz = np.polyfit(final['Angle'], final['Frame']*calib_step, poly_degree)
    dd = np.polyfit(final['Angle'], final['avg_distance'], poly_degree)
    dr = np.polyfit(final['Angle'], final['ratio'],poly_degree)
    angle_range = (np.min(final['Angle']), np.max(final['Angle']))
    
    fitting_paras = [dx, dy, dz, dd, dr, angle_range]
    return fitting_paras

# Analysis function
def DH_analysis(data_name, general_paras, fitting_paras, filter_paras):
    cols = ['#Frame', 'X (px)', 'Y (px)', 'Intensity (photon)', 'Precision (nm)'] # useful column names
    raw_data = pd.read_csv(data_name, skiprows=8, sep = '\t') # read xls file
    data = raw_data[cols] # take useful data
    data = data[data['Precision (nm)'] < general_paras[1]] # filter by precision
    frame_list = data['#Frame']
    frames = list(frame_list.drop_duplicates())
    
    # Seperate data into chunks based on frame number
    data_chunks = []
    for frame in frames:
        data_chunks.append(data[data['#Frame'] == frame])
    
    Frames = []
    Distances = []
    X = []
    Y = []
    Angle = []
    Ratio = []
    Intensity = []
    
    # for each chunk of data: (for each frame)
    for i in range(0, len(frames)):
        frame = frames[i]
        frame_data = data_chunks[i]   # Take data of frame #i
        xs = frame_data['X (px)']
        ys = frame_data['Y (px)']
        ints = frame_data['Intensity (photon)']
        
        xs_arr = np.reshape(xs.values, [len(xs), 1]) # Change DataFrame into numpy arrays
        ys_arr = np.reshape(ys.values, [len(ys), 1])
        ints_arr = np.reshape(ints.values, [len(ints), 1])
        
        D = squareform(pdist(np.concatenate([xs_arr, ys_arr], axis = 1), 'euclidean')) # Calculate distances between any two locs
        D[D>general_paras[2][1]] = 0 # filter out pairs with distances over/below expected values
        D[D<general_paras[2][0]] = 0
        D = np.triu(D)   # Take upper triangle of the matrix
        [I, J] = np.nonzero(D)   # Take indices of chosen pairs of localisations
        dists = D[I,J]
        avg_x_all = []
        avg_y_all = []
        ratio_all = []
        intensity_all = []
        angles_all = []
        frame_all = list(np.ones(len(I))*frame)
        for k in range(0, len(I)):
            ratio = np.max([ints_arr[I[k]], ints_arr[J[k]]]) / np.min([ints_arr[I[k]], ints_arr[J[k]]]) # Calculate intensity ratio between two dots in a pair
            intensity = np.mean([ints_arr[I[k]], ints_arr[J[k]]]) # Calculate mean intensity of the two dots
            angle = math.atan2(ys_arr[J[k]]-ys_arr[I[k]], xs_arr[J[k]]-xs_arr[I[k]]) # Calculate angle of the DH 
            avg_x = (xs_arr[I[k]] + xs_arr[J[k]])/2 # Determine x pos of the fluorophore
            avg_x = avg_x[0]
            avg_y = (ys_arr[I[k]] + ys_arr[J[k]])/2 # Determine y pos of the fluorophore
            avg_y = avg_y[0]
            if angle < 0:
                angle = angle+np.pi
            ratio_all.append(ratio)
            intensity_all.append(intensity)
            angles_all.append(angle)
            avg_x_all.append(avg_x)
            avg_y_all.append(avg_y)
        
        Frames.append(frame_all)
        Distances.append(list(dists))
        X.append(avg_x_all)
        Y.append(avg_y_all)
        Angle.append(angles_all)
        Ratio.append(ratio_all)
        Intensity.append(intensity_all)
    
    Frame_arr = np.concatenate(Frames)
    Distances_arr = np.concatenate(Distances)
    X_arr = np.concatenate(X)
    Y_arr = np.concatenate(Y)
    Angle_arr = np.concatenate(Angle)
    Ratio_arr = np.concatenate(Ratio)
    Intensity_arr = np.concatenate(Intensity)
    
    zN = np.polyval(fitting_paras[2], Angle_arr) # Calculate all Z values of fluorophores in data using their angles
    #xN = (X_arr - np.polyval(fitting_paras[0], zN)) * general_paras[0] # Calibrate X by Z, convert into nm    
    #yN = (Y_arr - np.polyval(fitting_paras[1], zN)) * general_paras[0] # Calibrate Y by Z, convert into nm
    
    zmin = np.polyval(fitting_paras[2], fitting_paras[5][0]) # Minimum angle from calibration
    zN = np.polyval(fitting_paras[2], Angle_arr) # Calculate all Z values of fluorophores in data using their angles
    xN = (X_arr - (np.polyval(fitting_paras[0], zN) - np.polyval(fitting_paras[0], zmin) )) * general_paras[0] # Calibrate X by Z, convert into nm
    yN = (Y_arr - (np.polyval(fitting_paras[1], zN) - np.polyval(fitting_paras[1], zmin) )) * general_paras[0] # Calibrate Y by Z, convert into nm
    
    
    marker = np.ones((len(Frame_arr), 1))
    if filter_paras[0]:
        if filter_paras[1]:
            low_angle = [i for i in range(len(Angle_arr)) if Angle_arr[i] < fitting_paras[5][0]]
            high_angle = [i for i in range(len(Angle_arr)) if Angle_arr[i] > fitting_paras[5][1]]
            if low_angle:
                marker[np.asarray(low_angle)] = -1
            if high_angle:
                marker[np.asarray(high_angle)] = -1
        if filter_paras[2]:
            diff_dist = np.abs(Distances_arr-np.polyval(fitting_paras[3], Angle_arr))/Distances_arr
            large_distance = [i for i in range(len(diff_dist)) if diff_dist[i] > filter_paras[3]]
            if large_distance:
                marker[np.asarray(large_distance)] = -1
        if filter_paras[4]:
            diff_ratio = np.abs(Ratio_arr - np.polyval(fitting_paras[4], Angle_arr))
            bad_ratio = [i for i in range(len(diff_ratio)) if diff_ratio[i] > filter_paras[5]]
            bad_ratio2 = [i for i in range(len(Ratio_arr)) if Ratio_arr[i] > 4]
            if bad_ratio:
                marker[np.asarray(bad_ratio)] = -1
            if bad_ratio2:
                marker[np.asarray(bad_ratio2)] = -1
    
    raw = pd.DataFrame({'Frame': Frame_arr,
                        'Distance': Distances_arr,
                        'X': xN,
                        'Y': yN,
                        'Z': zN,
                        'Angle': Angle_arr,
                        'Ratio': Ratio_arr,
                        'Intensity': Intensity_arr
                        })
    
    drop_row = [i for i in range(len(marker)) if marker[i] < 0]
    filtered = raw.drop(drop_row)
    return filtered

def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],y[i], ha = 'center', va = 'bottom')

def run_analysis(path, report, frame_number, new_fname_list, histo_list, general_paras, fitting_paras, filter_paras, poly_degree):
    num_files = report.shape[0]
    for i in range(num_files):
        data_name = report.Filename[i]
        filtered_data = DH_analysis(data_name, general_paras, fitting_paras, filter_paras)
        threeD_data = filtered_data.copy()
        threeD_data = threeD_data[['X', 'Y', 'Z', 'Intensity', 'Frame']]
        report.Num_of_locs[i] = filtered_data.shape[0]
        #filtered_data.to_csv(new_fname_list[i], sep = "\t", index = None)
        threeD_data.to_csv(new_fname_list[i].replace('csv','3d'), sep='\t', index = None, header = None)
        plt.figure()
        plt.hist(filtered_data.Frame, bins = np.min((filtered_data.shape[0], 50))) # default number of bins = 50
        plt.title('Slice '+ report.Slice[i])
        plt.xlabel('Frame number')
        plt.ylabel('Number of localisations')
        #plt.savefig(histo_list[i])
    '''
    report_sort = report.sort_values(['Pos', 'Slice'], ascending = [True, True])
    report_sort['Num_locs_per_frame'] = report_sort.Num_of_locs / frame_number
    #report_sort.to_csv(path + '\\Report_'+ path.split('\\')[-1]+'.csv', index = None)
    pos = np.unique(np.asarray(report_sort.Pos))
    for p in pos:
        pos_imgname = path + '\\Pos ' + str(p)  +'.png'
        pos_imgname2 = path + '\\Pos ' + str(p)  +'_per_frame.png'
        
        pos_slice = report_sort[report_sort.Pos == str(p)].Slice
        num_locs = report_sort[report_sort.Pos == str(p)].Num_of_locs
        
        plt.figure()
        plt.bar(pos_slice, num_locs)
        addlabels(pos_slice, np.asarray(num_locs))
        plt.title('Position '+ str(p))
        plt.xlabel('Slice number')
        plt.ylabel('Number of localisations')
        #plt.savefig(pos_imgname)
        plt.close()
        
        plt.figure()
        plt.bar(pos_slice, np.asarray(num_locs)/frame_number)
        addlabels(pos_slice, np.asarray(num_locs)/frame_number)
        plt.title('Position '+ str(p))
        plt.xlabel('Slice number')
        plt.ylabel('Number of localisations / Frame')
        #plt.savefig(pos_imgname2)
        plt.close()
    '''
    return 0


#%%
report, new_fname_list, histo_list = get_data(path, ext)
fitting_paras = DH_calibration(calib_name, calib_step, initial_distance_filter,poly_degree, fitting_mode, range_to_fit)
run_analysis(path, report, frame_number, new_fname_list, histo_list, general_paras, fitting_paras, filter_paras,poly_degree)
