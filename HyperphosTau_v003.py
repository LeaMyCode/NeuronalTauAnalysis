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
import math
import os
import csv
import pandas as pd
import statistics
from read_roi import read_roi_zip
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
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
    neuron_data = {
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
            'IntDen_dendrite_AT8' : [],
            'IntDen_soma_AT8' : [],
            'hyperphos_area_whole' : [],
            'hyperphos_area_whole_norm' : [],
            'mean_Int_whole_AT8' : [],
            'IntDen_whole_AT8' : []
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
            try:
                file_image = file_header.split(".tif")[0]
                image = io.imread(os.path.join(cells_path, str(file_image + "_MIP.tif")))
            
            except:
                file_image = file_header
                image = io.imread(os.path.join(cells_path, str(file_image)))
                
            AT8_image = image[:,:,AT8_channel]
            
            whole_img_intensity_AT8 = np.mean(AT8_image)
            
            area_img_AT8 = AT8_image.size 
            area_img_AT8 *= (resolution * resolution)

            
            
            #load rois, structure: [[[x,y][x,y]],[[x,y][x,y]]]
            nucleus_rois = readRoi(os.path.join(cells_path, str(file_header + "-nucleus_ROI.zip")))
            print("Nucleus rois: ", len(nucleus_rois))
            soma_rois = readRoi(os.path.join(cells_path, str(file_header + "-soma_ROI.zip")))
            print("Soma rois: ", len(soma_rois))
            dendrite_rois = readRoi(os.path.join(cells_path, str(file_header + "-dendrites_ROI.zip")))
            print("Dendrites rois: ", len(dendrite_rois))
            AT8_rois = readRoi(os.path.join(cells_path, str(file_header + "AT8_ROI.zip")))
            print("AT8 rois: ", len(AT8_rois))
            
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
            
            # Define the minimum and maximum area limits for AT8 ROIs
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
            
            # Initialize counters and area accumulators
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
                    
            # Get intensity values inside nucleus
            roi_intensity_values = get_intensity_values(nucleus_polygon, AT8_image)
            nucleus_intensities_AT8 = np.mean(roi_intensity_values)
            
                    
            
            # Check which AT8 ROIs are inside the nucleus and accumulate the area
            for at8 in filtered_AT8_polygons:
                if soma_region.contains(at8):
                    AT8_in_soma_count += 1
                    AT8_in_soma_area += at8.area  
                    
            # Get intensity values inside the ROI
            roi_intensity_values = get_intensity_values(soma_region, AT8_image)
            soma_intensities_AT8 = np.mean(roi_intensity_values)
                    
            
            dendrite_intensities_AT8 = []
            
            # Check which AT8 ROIs are inside the dendrites and accumulate the area
            for dendrite_polygon in dendrite_polygons:
                # Get intensity values inside the ROI
                roi_intensity_values = get_intensity_values(dendrite_polygon, AT8_image)
                dendrite_intensities_AT8.append(np.mean(roi_intensity_values))
                
                for at8 in filtered_AT8_polygons:
                    try:
                        if dendrite_polygon.contains(at8):
                            AT8_in_dendrites_count += 1
                            AT8_in_dendrites_area += at8.area
                            
                    except:
                        print(f'In one of the dendrites of {file_header} is no hyperphos punct.')
                        
            
            # Multiply areas by resolution if necessary (for scaling units)
            AT8_in_nucleus_area *= (resolution * resolution)
            AT8_in_soma_area *= (resolution * resolution)
            AT8_in_dendrites_area *= (resolution * resolution)
            
            # Results
            print(f"Number of AT8 ROIs inside the nucleus: {AT8_in_nucleus_count}")
            print(f"Total area of AT8 ROIs inside the nucleus: {AT8_in_nucleus_area}")
            
            print(f"Number of AT8 ROIs inside the dendrites: {AT8_in_dendrites_count}")
            print(f"Total area of AT8 ROIs inside the dendrites: {AT8_in_dendrites_area}")
                
            print("IntDen dendrite: ", np. mean(dendrite_intensities_AT8)/dendrite_area)
            
            # fill the dicitionary for saving
            neuron_data['file'].append(file_header)
            neuron_data['hyperphosPuncta_nucleus'].append(AT8_in_nucleus_count)
            neuron_data['hyperphosPuncta_soma'].append(AT8_in_soma_count)
            neuron_data['hyperphosPuncta_dendrite'].append(AT8_in_dendrites_count) 
            neuron_data['hyperphos_area_nucleus'].append(AT8_in_nucleus_area)
            neuron_data['hyperphos_area_soma'].append(AT8_in_soma_area)
            neuron_data['hyperphos_area_dendrite'].append(AT8_in_dendrites_area)
            neuron_data['area_nucleus'].append(nucleus_area)
            neuron_data['area_soma'].append(soma_area)
            neuron_data['area_dendrite'].append(dendrite_area)
            neuron_data['total_area'].append(nucleus_area+soma_area+dendrite_area)
            neuron_data['totalhyperphosPuncta_norm'].append((AT8_in_nucleus_count+AT8_in_dendrites_count+AT8_in_soma_count)/(nucleus_area+dendrite_area+soma_area))
            neuron_data['hyperphosPuncta_nucleus_norm'].append(AT8_in_nucleus_count/nucleus_area)
            neuron_data['hyperphosPuncta_soma_norm'].append(AT8_in_soma_count/soma_area)
            neuron_data['hyperphosPuncta_dendrite_norm'].append(AT8_in_dendrites_count/dendrite_area) 
            neuron_data['totalhyperphosArea_norm'].append((AT8_in_nucleus_area+AT8_in_dendrites_area+AT8_in_soma_area)/(nucleus_area+dendrite_area+soma_area))
            neuron_data['hyperphosArea_nucleus_norm'].append(AT8_in_nucleus_area/nucleus_area)
            neuron_data['hyperphosArea_soma_norm'].append(AT8_in_soma_area/soma_area)
            neuron_data['hyperphosArea_dendrite_norm'].append(AT8_in_dendrites_area/dendrite_area)
            neuron_data['mean_Int_nucleus_AT8'].append(np.mean(nucleus_intensities_AT8))
            neuron_data['mean_Int_soma_AT8'].append(np.mean(soma_intensities_AT8))
            neuron_data['mean_Int_dendrite_AT8'].append(np.mean(dendrite_intensities_AT8))
            neuron_data['IntDen_nucleus_AT8'].append(np.mean(nucleus_intensities_AT8)*nucleus_area)
            neuron_data['IntDen_soma_AT8'].append(np.mean(soma_intensities_AT8)*soma_area)
            neuron_data['IntDen_dendrite_AT8'].append(np.mean(dendrite_intensities_AT8)*dendrite_area)
            neuron_data['hyperphos_area_whole'].append(AT8_total_area)
            neuron_data['hyperphos_area_whole_norm'].append(AT8_total_area/area_img_AT8)
            neuron_data['mean_Int_whole_AT8'].append(whole_img_intensity_AT8)
            neuron_data['IntDen_whole_AT8'].append(whole_img_intensity_AT8 * area_img_AT8)
            
            
    

        
    #write condition data to .csv --> summary of all data
    
    df_data = pd.DataFrame(neuron_data)
    df_data.to_csv(results_path_condition + fr'/neuron_data_{condition_folder}.csv', index=False)

    

   


  
            
            
            
   





    
    
    

 