# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 07:25:58 2024

@author: leaga
"""



# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 08:43:08 2021

@author: leaga
"""

import numpy as np
from tkinter import filedialog
import os
import pandas as pd
from read_roi import read_roi_zip
from shapely.geometry import Polygon, Point
from skimage import io
from shapely.geometry import Polygon, MultiPolygon      


####################################################################################################################################
#
# Please enter your resolution!
#
####################################################################################################################################

global resolution 
resolution = 0.1221896

####################################################################################################################################
#
# Please ensure the channel are correctly set
#
####################################################################################################################################

AT8_channel = 1
totalTau_channel = 0

####################################################################################################################################
#
# Functions
#
####################################################################################################################################

#get the header of your file names to identify connected data
def fileNameGetter(s):
    name = s.split('.tif')
    new_name = name[0]+ ".tif"
    new_name.rstrip()
    return new_name

# Read the coordintates for each ROI in your ROIs
def readRoi(roi_zip_path):
    # Create a list to store the result
    roi_list = []
    
    #open the rois
    rois = read_roi_zip(roi_zip_path)

    for key, value in rois.items():
        xy_pairs = [[float(x), float(y)] for x, y in zip(value['x'], value['y'])]
        roi_list.append(xy_pairs)
        
    
    return roi_list


# Helper function to extract pixel values inside a polygon
def get_intensity_values(polygon, intensity_image):
    # Create a mask for the polygon
    min_x, min_y, max_x, max_y = map(int, polygon.bounds)  # Get bounding box
    
    # Clip bounds to ensure they fit within the image dimensions
    min_x = max(min_x, 0)
    min_y = max(min_y, 0)
    max_x = min(max_x, intensity_image.shape[1] - 1)  # Limit to image width (columns)
    max_y = min(max_y, intensity_image.shape[0] - 1)  # Limit to image height (rows)

    mask = np.zeros_like(intensity_image, dtype=bool)

    for x in range(min_x, max_x + 1):
        for y in range(min_y, max_y + 1):
            if polygon.contains(Polygon([(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)])):
                mask[y, x] = True  # Flip to match numpy's row-col indexing (y first)

    return intensity_image[mask]

####################################################################################################################################
#
# MAIN
#
####################################################################################################################################


root_path = filedialog.askdirectory(title="Please select your folder with images to analyze")


for condition_folder in os.listdir(root_path):
    
    #create dir for your condition results
    results_path_condition = os.path.join(root_path,"Analysis")
    if not os.path.exists(results_path_condition):
        os.makedirs(results_path_condition)
        
    if condition_folder == "Analysis":
        continue
    
    condition_path = os.path.join(root_path,condition_folder)
    print("Condition folder : ", condition_folder)
    
    #initialize analysis factors and set to zero for each condition --> summarized data  
    neuron_data_AT8 = {
            'file': [],
            'hyperphosPuncta_nucleus' : [],
            'hyperphosPuncta_soma' : [],
            'hyperphosPuncta_dendrite' : [],
            'hyperphos_area_nucleus': [],
            'hyperphos_area_soma' : [],
            'hyperphos_area_dendrite': [],
            'area_nucleus': [],
            'area_soma' : [],
            'area_dendrite': [],
            'total_area' : [],
            'totalhyperphosPuncta_norm' : [],
            'hyperphosPuncta_nucleus_norm': [],
            'hyperphosPuncta_soma_norm' : [],
            'hyperphosPuncta_dendrite_norm': [],
            'totalhyperphosArea_norm' : [],
            'hyperphosArea_nucleus_norm': [],
            'hyperphosArea_soma_norm' : [],
            'hyperphosArea_dendrite_norm': [],
            'mean_Int_nucleus_AT8' : [],
            'mean_Int_soma_AT8' : [],
            'mean_Int_dendrite_AT8' : [],
            'IntDen_nucleus_AT8' : [],
            'IntDen_soma_AT8' : [],
            'IntDen_dendrite_AT8' : [],
            'hyperphos_area_whole' : [],
            'hyperphos_area_whole_norm' : [],
            'mean_Int_whole_AT8' : [],
            'IntDen_whole_AT8' : []
        }
    
    neuron_data_tau = {
            'file': [],
            'tauPuncta_nucleus' : [],
            'tauPuncta_soma' : [],
            'tauPuncta_dendrite' : [],
            'tau_area_nucleus': [],
            'tau_area_soma' : [],
            'tau_area_dendrite': [],
            'area_nucleus': [],
            'area_soma' : [],
            'area_dendrite': [],
            'total_area' : [],
            'totaltauPuncta_norm' : [],
            'tauPuncta_nucleus_norm': [],
            'tauPuncta_soma_norm' : [],
            'tauPuncta_dendrite_norm': [],
            'totaltauArea_norm' : [],
            'tauArea_nucleus_norm': [],
            'tauArea_soma_norm' : [],
            'tauArea_dendrite_norm': [],
            'mean_Int_nucleus_tau' : [],
            'mean_Int_soma_tau' : [],
            'mean_Int_dendrite_tau' : [],
            'IntDen_nucleus_tau' : [],
            'IntDen_soma_tau' : [],
            'IntDen_dendrite_tau' : [],
            'tau_area_whole' : [],
            'tau_area_whole_norm' : [],
            'mean_Int_whole_tau' : [],
            'IntDen_whole_tau' : []
        }
    
    tau_AT8_data = {
            'file': [],
            'ratio_Puncta_nucleus' : [],
            'ratio_Puncta_soma' : [],
            'ratio_Puncta_dendrite' : [],
            'ratio_area_nucleus' : [],
            'ratio_area_soma' : [],
            'ratio_area_dendrite' : [],
            'ratio_mean_int_nucleus': [],
            'ratio_mean_int_soma' : [],
            'ratio_mean_int_dendrite' : [],
            'ratio_IntDen_nucleus' : [],
            'ratio_IntDen_soma': [],
            'ratio_IntDen_dendrite' : [],
            'ratio_area_whole' : [],
            'ratio_mean_Int_whole' : [],
            'ratio_IntDen_whole' : []
            
            }
    
    

    for round_folder in os.listdir(condition_path):
        
        if round_folder == "Analysis":
            continue

        cells_path = os.path.join(condition_path,round_folder)
            
  
        print("Round folder: ", round_folder)

        
        for files in os.listdir(cells_path):
            if files == "Analysis":
                continue
             
            # ensure to have each cell only once
            if not files.endswith("-nucleus_ROI.zip"):
                continue
            
            
            #get your files header to identify connected data
            file_header = fileNameGetter(files)
            print("file_header: ", file_header)
            
            # Load the 3-channel image to calculate the intensities of the staining
            file_image = file_header.split(".tif")[0]
            image = io.imread(os.path.join(cells_path, str(file_image + "_MIP.tif")))
            AT8_image = image[:,:,AT8_channel]
            tau_image = image[:,:,totalTau_channel]
            
            whole_img_intensity_AT8 = np.mean(AT8_image)
            whole_img_intensity_tau = np.mean(tau_image)
            
            area_img_AT8 = AT8_image.size 
            area_img_tau = tau_image.size
            
            #load rois, structure: [[[x,y][x,y]],[[x,y][x,y]]]
            nucleus_rois = readRoi(os.path.join(cells_path, str(file_header + "-nucleus_ROI.zip")))
            print("Nucleus rois: ", len(nucleus_rois))
            soma_rois = readRoi(os.path.join(cells_path, str(file_header + "-soma_ROI.zip")))
            print("Soma rois: ", len(soma_rois))
            dendrite_rois = readRoi(os.path.join(cells_path, str(file_header + "-dendrites_ROI.zip")))
            print("Dendrites rois: ", len(dendrite_rois))
            AT8_rois = readRoi(os.path.join(cells_path, str(file_header + "AT8_ROI.zip")))
            print("AT8 rois: ", len(AT8_rois))
            tau_rois = readRoi(os.path.join(cells_path, str(file_header + "totalTau_ROI.zip")))
            print("Tau rois: ", len(tau_rois))
            
            #calculate the area of the nucleus
            nucleus_polygon = Polygon(nucleus_rois[0]) # we only have one nucleus
            nucleus_area = nucleus_polygon.area * (resolution*resolution)
            #print("nucleus_area: ", nucleus_area)
            
            #calculate the area of the soma
            soma_polygon = Polygon(soma_rois[0])
            soma_region = soma_polygon.difference(nucleus_polygon)
            # Compute area: Sum all parts if it's a MultiPolygon, otherwise just take the area
            if isinstance(soma_region, MultiPolygon):
                soma_area = sum(p.area for p in soma_region.geoms) * (resolution * resolution)
            else:
                soma_area = soma_region.area * (resolution * resolution)
            
            #calculate the area of the dendrites
            
            # Create individual polygons for each ROI
            dendrite_polygons = [Polygon(roi) for roi in dendrite_rois]   
            
            dendrite_area = 0
            for dendrite in dendrite_polygons:
                dendrite_area += dendrite.area * (resolution*resolution)
                
            #print("Dendrite area: ", dendrite_area)
            
            # Define the minimum and maximum area limits for AT8 & total Tau ROIs
            min_size = 2
            max_size = 200
            
            # Convert AT8 ROIs to polygons or points and filter by area
            filtered_AT8_polygons = []
            for roi in AT8_rois:
                # Create Polygon or Point based on the ROI structure
                at8_polygon = Polygon(roi) if len(roi) > 2 else Point(roi[0])
                
                # Calculate area if it's a polygon
                at8_area = at8_polygon.area if isinstance(at8_polygon, Polygon) else 0
                
                # Filter based on size limits
                if min_size <= at8_area <= max_size:
                    filtered_AT8_polygons.append(at8_polygon)
                    
            # Convert total tau ROIs to polygons or points and filter by area
            filtered_tau_polygons = []
            for roi in tau_rois:
                # Create Polygon or Point based on the ROI structure
                tau_polygon = Polygon(roi) if len(roi) > 2 else Point(roi[0])
                
                # Calculate area if it's a polygon
                tau_area = tau_polygon.area if isinstance(tau_polygon, Polygon) else 0
                
                # Filter based on size limits
                if min_size <= tau_area <= max_size:
                    filtered_tau_polygons.append(tau_polygon)
                    
            
            
            # Initialize counters and area accumulators for AT8
            AT8_in_nucleus_count = 0
            AT8_in_nucleus_area = 0
            
            AT8_in_soma_count = 0
            AT8_in_soma_area = 0
            
            AT8_in_dendrites_count = 0
            AT8_in_dendrites_area = 0
            
            # Calculate the total AT8 area
            AT8_total_area = sum(polygon.area for polygon in filtered_AT8_polygons)
            
            # Check which AT8 ROIs are inside the nucleus and accumulate the area
            for at8 in filtered_AT8_polygons:
                if nucleus_polygon.contains(at8):
                    AT8_in_nucleus_count += 1
                    AT8_in_nucleus_area += at8.area  
                    
            # Get intensity values inside the ROI
            roi_intensity_values = get_intensity_values(nucleus_polygon, AT8_image)
            nucleus_intensities_AT8 = np.mean(roi_intensity_values)
            
            
            # Check which AT8 ROIs are inside the soma and accumulate the area
            for at8 in filtered_AT8_polygons:
                if soma_region.contains(at8):
                    AT8_in_soma_count += 1
                    AT8_in_soma_area += at8.area  
                    
            # Get intensity values inside the ROI
            roi_intensity_values = get_intensity_values(soma_region, AT8_image)
            soma_intensities_AT8 = np.mean(roi_intensity_values)
            
            #initialize lists for dendrite intensities (more than one polygon cannot be calculated and once and combining non-touching polygons would crate false data)
            dendrite_intensities_AT8 = []
            
            # Check which AT8 ROIs are inside the dendrites and accumulate the area
            for dendrite_polygon in dendrite_polygons:
                roi_intensity_values = get_intensity_values(dendrite_polygon, AT8_image)
                dendrite_intensities_AT8.append(np.mean(roi_intensity_values))
                for at8 in filtered_AT8_polygons:
                    try:
                        if dendrite_polygon.contains(at8):
                            AT8_in_dendrites_count += 1
                            AT8_in_dendrites_area += at8.area
                            
                            # Get intensity values inside the ROI
                    except:
                        print(f'In one of the dendrites of {file_header} is no hyperphos punct.')
                        
            
            tau_in_nucleus_count = 0
            tau_in_nucleus_area = 0
            
            tau_in_soma_count = 0
            tau_in_soma_area = 0
            
            tau_in_dendrites_count = 0
            tau_in_dendrites_area = 0
            
            
            # Calculate the total tau area
            tau_total_area = sum(polygon.area for polygon in filtered_tau_polygons)
            
            # Check which tau ROIs are inside the nucleus and accumulate the area
            for tau in filtered_tau_polygons:
                if nucleus_polygon.contains(tau):
                    tau_in_nucleus_count += 1
                    tau_in_nucleus_area += tau.area  
                    
            # Get intensity values inside the ROI
            roi_intensity_values = get_intensity_values(nucleus_polygon, tau_image)
            nucleus_intensities_tau = np.mean(roi_intensity_values)
            
            
            # Check which tau ROIs are inside the soma and accumulate the area
            for tau in filtered_tau_polygons:
                if soma_region.contains(tau):
                    tau_in_soma_count += 1
                    tau_in_soma_area += tau.area  
                    
            # Get intensity values inside the ROI
            roi_intensity_values = get_intensity_values(soma_region, tau_image)
            soma_intensities_tau = np.mean(roi_intensity_values)
            
            #initlalize list for dendrite intensities of total tau
            dendrite_intensities_tau = []
            
            # Check which tau ROIs are inside the dendrites and accumulate the area
            for dendrite_polygon in dendrite_polygons:
                roi_intensity_values = get_intensity_values(dendrite_polygon, tau_image)
                dendrite_intensities_tau.append(np.mean(roi_intensity_values))
                for tau in filtered_tau_polygons:
                    try:
                        if dendrite_polygon.contains(tau):
                            tau_in_dendrites_count += 1
                            tau_in_dendrites_area += tau.area
                    
                    except:
                        print(f'In one of the dendrites of {file_header} is no hyperphos punct.')
            
            
            # Multiply areas by resolution if necessary (for scaling units)
            AT8_in_nucleus_area *= (resolution * resolution)
            AT8_in_soma_area *= (resolution * resolution)
            AT8_in_dendrites_area *= (resolution * resolution)
            
            # Multiply areas by resolution if necessary (for scaling units)
            tau_in_nucleus_area *= (resolution * resolution)
            tau_in_soma_area *= (resolution * resolution)
            tau_in_dendrites_area *= (resolution * resolution)
            
            area_img_AT8 *= (resolution * resolution)
            area_img_tau *= (resolution * resolution)
            
            
            # Results
            print(f"Number of AT8 ROIs inside the nucleus: {AT8_in_nucleus_count}")
            print(f"Total area of AT8 ROIs inside the nucleus: {AT8_in_nucleus_area}")
            
            print(f"Number of AT8 ROIs inside the dendrites: {AT8_in_dendrites_count}")
            print(f"Total area of AT8 ROIs inside the dendrites: {AT8_in_dendrites_area}")
                
            print("IntDen dendrite: ", np.mean(dendrite_intensities_AT8)*dendrite_area)
 
            
            # fill the dicitionary with AT8 data for saving
            neuron_data_AT8['file'].append(file_header)
            neuron_data_AT8['hyperphosPuncta_nucleus'].append(AT8_in_nucleus_count)
            neuron_data_AT8['hyperphosPuncta_soma'].append(AT8_in_soma_count)
            neuron_data_AT8['hyperphosPuncta_dendrite'].append(AT8_in_dendrites_count) 
            neuron_data_AT8['hyperphos_area_nucleus'].append(AT8_in_nucleus_area)
            neuron_data_AT8['hyperphos_area_soma'].append(AT8_in_soma_area)
            neuron_data_AT8['hyperphos_area_dendrite'].append(AT8_in_dendrites_area)
            neuron_data_AT8['area_nucleus'].append(nucleus_area)
            neuron_data_AT8['area_soma'].append(soma_area)
            neuron_data_AT8['area_dendrite'].append(dendrite_area)
            neuron_data_AT8['total_area'].append(nucleus_area + soma_area + dendrite_area)
            neuron_data_AT8['totalhyperphosPuncta_norm'].append((AT8_in_nucleus_count+AT8_in_dendrites_count+AT8_in_soma_count)/(nucleus_area+dendrite_area+soma_area))
            neuron_data_AT8['hyperphosPuncta_nucleus_norm'].append(AT8_in_nucleus_count/nucleus_area)
            neuron_data_AT8['hyperphosPuncta_soma_norm'].append(AT8_in_soma_count/soma_area)
            neuron_data_AT8['hyperphosPuncta_dendrite_norm'].append(AT8_in_dendrites_count/dendrite_area) 
            neuron_data_AT8['totalhyperphosArea_norm'].append((AT8_in_nucleus_area+AT8_in_dendrites_area+AT8_in_soma_area)/(nucleus_area+dendrite_area+soma_area))
            neuron_data_AT8['hyperphosArea_nucleus_norm'].append(AT8_in_nucleus_area/nucleus_area)
            neuron_data_AT8['hyperphosArea_soma_norm'].append(AT8_in_soma_area/soma_area)
            neuron_data_AT8['hyperphosArea_dendrite_norm'].append(AT8_in_dendrites_area/dendrite_area)
            neuron_data_AT8['mean_Int_nucleus_AT8'].append(np.mean(nucleus_intensities_AT8))
            neuron_data_AT8['mean_Int_soma_AT8'].append(np.mean(soma_intensities_AT8))
            neuron_data_AT8['mean_Int_dendrite_AT8'].append(np.mean(dendrite_intensities_AT8))
            neuron_data_AT8['IntDen_nucleus_AT8'].append(np.mean(nucleus_intensities_AT8)*nucleus_area)
            neuron_data_AT8['IntDen_dendrite_AT8'].append(np.mean(dendrite_intensities_AT8)*dendrite_area)
            neuron_data_AT8['IntDen_soma_AT8'].append(np.mean(soma_intensities_AT8)*soma_area)
            neuron_data_AT8['hyperphos_area_whole'].append(AT8_total_area)
            neuron_data_AT8['hyperphos_area_whole_norm'].append(AT8_total_area/area_img_AT8)
            neuron_data_AT8['mean_Int_whole_AT8'].append(whole_img_intensity_AT8)
            neuron_data_AT8['IntDen_whole_AT8'].append(whole_img_intensity_AT8 * area_img_AT8)

            
            # fill the dicitionary with tau data for saving
            neuron_data_tau['file'].append(file_header)
            neuron_data_tau['tauPuncta_nucleus'].append(tau_in_nucleus_count)
            neuron_data_tau['tauPuncta_soma'].append(tau_in_soma_count)
            neuron_data_tau['tauPuncta_dendrite'].append(tau_in_dendrites_count) 
            neuron_data_tau['tau_area_nucleus'].append(tau_in_nucleus_area)
            neuron_data_tau['tau_area_soma'].append(tau_in_soma_area)
            neuron_data_tau['tau_area_dendrite'].append(tau_in_dendrites_area)
            neuron_data_tau['area_nucleus'].append(nucleus_area)
            neuron_data_tau['area_soma'].append(soma_area)
            neuron_data_tau['area_dendrite'].append(dendrite_area)
            neuron_data_tau['total_area'].append(nucleus_area + soma_area + dendrite_area)
            neuron_data_tau['totaltauPuncta_norm'].append((tau_in_nucleus_count+tau_in_dendrites_count+tau_in_soma_count)/(nucleus_area+dendrite_area+soma_area))
            neuron_data_tau['tauPuncta_nucleus_norm'].append(tau_in_nucleus_count/nucleus_area)
            neuron_data_tau['tauPuncta_soma_norm'].append(tau_in_soma_count/soma_area)
            neuron_data_tau['tauPuncta_dendrite_norm'].append(tau_in_dendrites_count/dendrite_area) 
            neuron_data_tau['totaltauArea_norm'].append((tau_in_nucleus_area+tau_in_dendrites_area+tau_in_soma_area)/(nucleus_area+dendrite_area+soma_area))
            neuron_data_tau['tauArea_nucleus_norm'].append(tau_in_nucleus_area/nucleus_area)
            neuron_data_tau['tauArea_soma_norm'].append(tau_in_soma_area/soma_area)
            neuron_data_tau['tauArea_dendrite_norm'].append(tau_in_dendrites_area/dendrite_area)
            neuron_data_tau['mean_Int_nucleus_tau'].append(np.mean(nucleus_intensities_tau))
            neuron_data_tau['mean_Int_soma_tau'].append(np.mean(soma_intensities_tau))
            neuron_data_tau['mean_Int_dendrite_tau'].append(np.mean(dendrite_intensities_tau))
            neuron_data_tau['IntDen_nucleus_tau'].append(np.mean(nucleus_intensities_tau)*nucleus_area)
            neuron_data_tau['IntDen_soma_tau'].append(np.mean(soma_intensities_tau)*soma_area)
            neuron_data_tau['IntDen_dendrite_tau'].append(np.mean(dendrite_intensities_tau)*dendrite_area)
            neuron_data_tau['tau_area_whole'].append(tau_total_area)
            neuron_data_tau['tau_area_whole_norm'].append(tau_total_area/area_img_tau)
            neuron_data_tau['mean_Int_whole_tau'].append(whole_img_intensity_tau)
            neuron_data_tau['IntDen_whole_tau'].append(whole_img_intensity_tau * area_img_tau)
            
            
        
   
            
            #fill data of ratio AT8 and total tau
            tau_AT8_data['file'].append(file_header)
            
            #somethimes there is no signal inside, so to still do the analysis, if there is nothing, append NaN
            
            #nucleus data
            try:
                tau_AT8_data['ratio_Puncta_nucleus'].append(AT8_in_nucleus_count/tau_in_nucleus_count)
                tau_AT8_data['ratio_area_nucleus'].append(AT8_in_nucleus_area/tau_in_nucleus_area)
                        
            except:
                tau_AT8_data['ratio_Puncta_nucleus'].append(np.nan)
                tau_AT8_data['ratio_area_nucleus'].append(np.nan)

            
            #soma data
            try:
                tau_AT8_data['ratio_Puncta_soma'].append(AT8_in_soma_count/tau_in_soma_count)
                tau_AT8_data['ratio_area_soma'].append(AT8_in_soma_area/tau_in_soma_area)       
            
            except:
                tau_AT8_data['ratio_Puncta_soma'].append(np.nan)
                tau_AT8_data['ratio_area_soma'].append(np.nan)
                
            #dendrite data
            try:
                tau_AT8_data['ratio_Puncta_dendrite'].append(AT8_in_dendrites_count/tau_in_dendrites_count)
                tau_AT8_data['ratio_area_dendrite'].append(AT8_in_dendrites_area/tau_in_dendrites_area)
                
            except:
                tau_AT8_data['ratio_Puncta_dendrite'].append(np.nan)
                tau_AT8_data['ratio_area_dendrite'].append(np.nan)
                
            
            #intensity data
            tau_AT8_data['ratio_mean_int_nucleus'].append(np.mean(nucleus_intensities_AT8)/np.mean(nucleus_intensities_tau))
            tau_AT8_data['ratio_IntDen_nucleus'].append((np.mean(nucleus_intensities_AT8)*nucleus_area)/(np.mean(nucleus_intensities_tau)*nucleus_area))
            
            tau_AT8_data['ratio_mean_int_soma'].append(np.mean(soma_intensities_AT8)/np.mean(soma_intensities_tau))
            tau_AT8_data['ratio_IntDen_soma'].append((np.mean(soma_intensities_AT8)*soma_area)/(np.mean(soma_intensities_tau)*soma_area))
            
            tau_AT8_data['ratio_mean_int_dendrite'].append(np.mean(dendrite_intensities_AT8)/np.mean(dendrite_intensities_tau))
            tau_AT8_data['ratio_IntDen_dendrite'].append((np.mean(dendrite_intensities_AT8)*dendrite_area)/(np.mean(dendrite_intensities_tau)*dendrite_area))
            
            
            tau_AT8_data['ratio_area_whole'].append(AT8_total_area/tau_total_area)
            tau_AT8_data['ratio_mean_Int_whole'].append(whole_img_intensity_AT8/whole_img_intensity_tau)
            tau_AT8_data['ratio_IntDen_whole'].append((whole_img_intensity_AT8 * area_img_AT8)/(whole_img_intensity_tau * area_img_tau))
            
        
    #write condition data to .csv --> summary of all data
    
    df_data1 = pd.DataFrame(neuron_data_AT8)
    df_data1.to_csv(results_path_condition + fr'/neuron_data_AT8{condition_folder}.csv', index=False)
    
    df_data2 = pd.DataFrame(neuron_data_tau)
    df_data2.to_csv(results_path_condition + fr'/neuron_data_tau{condition_folder}.csv', index=False)
    
    df_data3 = pd.DataFrame(tau_AT8_data)
    df_data3.to_csv(results_path_condition + fr'/neuron_data_ratios{condition_folder}.csv', index=False)
    
    
    

    

   


  
            
            
            
   





    
    
    

 