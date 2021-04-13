#!/usr/bin/python
# -*- coding: utf-8 -*-

###################################################################################################################################################

# Automatic pixel-based alignment tool for radio images for use in the creation of spectral ageing and spectral index maps.
# If you have made use of this script, please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803.

# Note this is the DS9 region file version. If you are using DS9 style region files there is a separate version available.

# See the readme for further usage details.

####################################################################################################################################################

# List of raised exceptions:
# ValueError
# FileNotFoundError
# IndexError
# NotADirectoryError
# ValueError

# Modes
# 0 = Predefined image
# 1 = Fixed reference
# 2 = Mean 

# Reference image
# Mode 0 = The index of the input list to use a reference. If the index points to more than one file e.g. ./*.fits the first file will be used, ordered alphabetically. Default 0.
# Mode 1 = Not applicable.
# Mode 2 = The index of the input list for which the mean should be calculated. To use all index use -1. Default 0.

# Initial setup:
# bia.Setup(input_files:list, region_files:list, output_dir="./", mode=2, reference_image=0, reference_location=[1024.0, 1024.0], residual_region=0, overwrite_files=False):

import traceback
import bratsalign as bia

# List of paths and file names of the images to align. Use '*' for all files. First element of list used for reference pixel determination.
input = ['./testimages/*.fits'] # Example: var_input = ['./path/to/some/fits/files/*.fits',  './path/to/some/more/fits/files/*.fits', './path/to/some/even/more/fits/files/*.fits']

# Path to directory for the output images
output = './output' 

# List of regions files. DS9 (.reg) format, one per input (can either be the same file or different if significantly shifted)
region = ['J1206_gaussfit.reg'] # Example: var_region = ['aregionfile.reg', './pathtoanother/regionfile.reg', 'aregionfile.reg']

# Parameters
reftype = 2 # Choose if we should use a reference point set by a predefined image (0), a fixed value (1), or the mean (2) 
refimage = 0 # Choose which image to align to. Only used if reftype = 0
ref = [1024.0, 1024.0] # List in [x, y] coodinates to align to. Only used if var_fixedref = 1

try:
    example = bia.Setup(input, region, output, overwrite_files=True, mode=2, reference_image=0)
    example.align()
except Exception as e:
    traceback.print_exc() # With traceback enabled
    #print('\n*** Exception: ' + str(e) + ' ***\n') # Or without traceback
