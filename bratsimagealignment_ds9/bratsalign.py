#!/usr/bin/env python3 
# -*- coding: utf-8 -*-

# Created by Jeremy. J. Harwood (2019)
# Contact: Jeremy.Harwood@physics.org or J.Harwood3@herts.ac.uk

# Please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803 if you have made use of this class.

# Automatically fits Gaussians and aligns images in pixel space. Particularly useful for spectral index and spectral age fitting on resolved scales.
# Change the var_ variables at the top of the script to suit your setup
# The script should be run in CASA using either exec() or the -C command line argument

####################################################################################################################################################

import os, glob, scipy.ndimage, sys, subprocess, pickle
import numpy as np
from astropy.io import fits

'''
TODO List:
Decide where log files should go. Same or different output folder?
If overwrite is false, decide how we check if output files will be overwritten early on, rather than at align time. Calculate the names and check at startup?
Write descriptive output for the log files
Check if aligned image name already exists (append int to end?)
'''


class _Constants:
    _VERSION = "v1.0.1"
    _MAX_MODES = 2
    _REFERENCE_ROUNDING_ACCURACY = 2
    _CASA_COMMANDS_FILE_NAME = "brats_casa_commands.py"
    _CASA_DATA_FILES_PREFIX = "brats_"
    _DEFAULT_OUTPUT_PATH = "./"
    _DEFAULT_CASA_PATH = "/soft/casapy/bin/casa"
    _CASA_HEADER = '''\
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This temporary script was created by the BRATS alignment tool and can be safely deleted once a run has completed.

# Created by Jeremy. J. Harwood (2019)
# Contact: Jeremy.Harwood@physics.org or J.Harwood3@herts.ac.uk

# Please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803 if you have made use of this script.

import os, sys, pickle
'''


class Setup(_Constants):

    def __init__(self, input_files:list, region_files:list, output_dir=_Constants._DEFAULT_OUTPUT_PATH, mode=2, reference_image=0, reference_location=[1024.0, 1024.0], residual_region=0, casa_path=_Constants._DEFAULT_CASA_PATH, overwrite_files=False):
        _Setters._set_input(self, input_files)
        _Setters._set_regions(self, region_files)
        _Setters._set_output(self, output_dir)
        _Setters._set_mode(self,mode)
        _Setters._set_reference_image(self,reference_image)
        _Setters._set_reference_location(self,reference_location)
        _Setters._set_casa_path(self,casa_path)
        self._residual_region=residual_region # THIS NEEDS A SETTER!
        self._overwrite = overwrite_files

    def set_mode(self, mode:int):
        _Setters._set_mode(self, mode)

    def set_reference_image(self, reference_image:int):
        _Setters._set_reference_image(self, reference_image)

    def set_reference_location(self, reference_location:int):
        _Setters._set_reference_location(self, reference_location)

    def set_casa_path(self, casa_path:str):
        self._casa_path = _set_casa_path(self, casa_path)

    def overwrite(self, overwrite_files:bool):
        self._overwrite = overwrite_files
    
    def align(self):
        _Core._align(self)
  

# Setter class for setting and updating parameters
class _Setters(Setup):

    def _set_input(self, __input_files):
        if len(__input_files) == 0:
            raise IndexError('Input files list cannot be empty!')
        else:
            self._input = []
            for path in __input_files:
                if not path.strip():
                   raise ValueError('The input files list contains a blank element. All elements must contain a valid value.')
                elif len(glob.glob(path)) == 0:
                    raise FileNotFoundError('Unable to locate any valid files matching %s.' % (file))
                else:
                    self._input.append(sorted(glob.glob(path)))

    def _set_output(self, __output_directory):
        if not __output_directory:
            raise EOFError('Output directory string cannot be empty!')
        elif not os.path.isdir(__output_directory):
            print('The output directory %s does not exist. Attempting to create it...' % (__output_directory))
            os.mkdir(__output_directory)
            #raise NotADirectoryError('The output directory %s does not exist. Please ensure the path is correct.' % (__output_directory))

        if __output_directory[len(__output_directory)-1] != '/':   # Standardise the path format to make maniulation safer and easier to manage
            __output_directory += '/'

        self._output = __output_directory

    def _set_regions(self, __region_files):
        if len(__region_files) == 1:
            if not __region_files[0].strip():
                   raise ValueError('The region files list contains a blank element. All elements must contain a valid value.')
            elif not os.path.isfile(__region_files[0]):
                   raise FileNotFoundError('The region file %s does not exist. Please ensure all paths are correct.' % (__region_files[0]))
            self._regions = __region_files * len(self._input)
        elif len(__region_files) == len(self._input):
            for file in __region_files:
               if not file.strip():
                   raise ValueError('The region files list contains a blank element. All elements must contain a valid value.')
               elif not os.path.isfile(file):
                   raise FileNotFoundError('The region file %s does not exist. Please ensure all paths are correct.' % (file))
            self._regions = __region_files
        elif len(__region_files) == 0:
            raise IndexError('Region files list cannot be empty!')
        else:
            raise IndexError('Region files list must be equal to either the number of images or to one [Regions: %i Images: %i]' % (len(__region_files), len(self._input)) )

    def _set_mode(self, __mode):
        if not 0 <= __mode <= self._MAX_MODES:
            raise ValueError('Modes must be in the range 0 and %i' % (self._MAX_MODES))
        self._mode = __mode

    def _set_reference_image(self,__reference_image):
        if self._mode == 2 and not -1 <= __reference_image < len(self._input):
            raise IndexError('Reference image must be an index between -1 and %d' % (len(self._input)-1))
        elif not 0 <= __reference_image < len(self._input):
            raise IndexError('Reference image must be an index between 0 and %d' % (len(self._input)-1))
        else:
            self._reference_image = __reference_image

    def _set_reference_location(self,__reference_location):
        if not len(__reference_location) == 2:
            raise IndexError('Reference location must be a list containing 2 element (currently %d)' % (len(__reference_location)))
        self._reference_location = __reference_location

    def  _set_casa_path(self, __casa_path): # This needs error checking!
        self._casa_path = __casa_path


class _Core(Setup):

    def _align (self):
        self.__commands = []
        self._reference_position = []
        self._aligned_files = []
        _Core._check_mode(self)
        # Do alignment
        _Core._define_position_parameters(self, self._input)
        _Core._get_reference_position(self)
        _Core._get_positions(self, "positions", True)
        _Core._casa_write(self)
        _Core.__casa_run(self)
        _Core._perform_shift(self)
        # Check alignment
        _Core._define_position_parameters(self, self._aligned_files)
        _Core._get_residual_positions(self, "aligned_positions", False)
        _Core._casa_write(self)
        _Core.__casa_run(self)
        _Core._get_residuals(self)


    def __casa_run(self):
        __casa_log_file = self._output + "/brats_align_casa.log"
        __command_line_arguements = [self._casa_path , "--nogui", "--nologger", "-c", self._CASA_COMMANDS_FILE_NAME, ">>", __casa_log_file]
        __casa_process = subprocess.Popen(__command_line_arguements, universal_newlines=True)
        stdout, stderr = __casa_process.communicate() # Use communicate rather than run or similar to prevent the risk of buffer overflows

    def _check_mode(self):
        if self._mode == 0 and self._reference_image == -1:
            raise IndexError('A reference image of -1 cannot be used for mode 1 (predefined image)')

    def _define_position_parameters(self, __inputfiles):
        self.__commands.append("__input = " + str(__inputfiles))
        self.__commands.append("__regions = " + str(self._regions))        

    def _get_reference_position(self):
        if self._mode == 0:
            print('Using a predefined image for the reference coordinates')
            __reference_fitting_file = self._input[self._reference_image][0]
            self.__commands.append("__gauss_fit = imfit(imagename=\'" + __reference_fitting_file + "\', region=\'" + str(self._regions[self._reference_image]) + "\', logfile=\'" + self._output + "reference_coords.txt\', dooff=True, append=False)")
            self.__commands.append("__reference_x, __reference_y = __gauss_fit['results']['component0']['pixelcoords'].round(" + str(self._REFERENCE_ROUNDING_ACCURACY) + ")")           
        elif self._mode == 1:
            print('Using a fixed reference of {:.2f}'.format(_reference_location[0]) + ', {:.2f}'.format(_reference_location[1]) + " for alignment")
            self.__commands.append("__reference_x, __reference_y = " + str(_reference_location))
        elif self._mode == 2:
            print('Using the mean value of all images for the reference coordinates')
            self.__commands.extend(["__sum_x = 0.0", "__sum_y = 0.0"])
            self.__commands.append("__sum_images = 0")

            if (self.set_reference_image == -1):
                self.__commands.append("for index in range(len(__input)):")
                self.__commands.append("\t__sum_images+=len(fits_location)")
                self.__commands.append("\tfor file in __input[index]:")
                self.__commands.append("\t\t__gauss_fit = imfit(imagename=file, region=__regions[index], logfile=\'" + self._output + "reference_coords.txt\', dooff=True)")
                self.__commands.append("\t\t__sum_x += __gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][0]")
                self.__commands.append("\t\t__sum_y += __gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][1]")
            else:
                self.__commands.append("__sum_images+=len(__input[" + str(self._reference_image) + "])")
                self.__commands.append("for file in __input[" + str(self._reference_image) + "]:")
                self.__commands.append("\t__gauss_fit = imfit(imagename=file, region=\'" + str(self._regions[self._reference_image]) + "\', logfile=\'" + self._output + "reference_coords.txt\', dooff=True)")
                self.__commands.append("\t__sum_x += __gauss_fit[\'results\'][\'component0\'][\'pixelcoords\'][0]")
                self.__commands.append("\t__sum_y += __gauss_fit[\'results\'][\'component0\'][\'pixelcoords\'][1]")
                
            self.__commands.extend(["__reference_x = (__sum_x/__sum_images)", "__reference_y = (__sum_y/__sum_images)"])
    
        self.__commands.append("print(\'Aligning to reference coordinates: {:.2f}\'.format(__reference_x) + ', {:.2f}\'.format(__reference_y))")

    def _get_positions(self, __return_file_name, __get_reference):
        self.__commands.append("print(\'Performing initial Gaussian fitting...\')")
        self.__commands.extend(["__x_positions = []", "__y_positions = []"])
        self.__commands.append("for index in range(len(__input)):")
        self.__commands.extend(["\t__x_positions.append([])", "\t__y_positions.append([])"])
        self.__commands.append("\tfor file in __input[index]:")
        self.__commands.append("\t\tif not file.upper().endswith(\'.FITS\'):")
        self.__commands.append("\t\t\tprint(\'The file %s does not contain a valid FITS extension. Skipping...\' % (file))")
        self.__commands.append("\t\t\tcontinue")
        __log_name = self._output + "original_coordinates"
        self.__commands.append("\t\t__gauss_fit = imfit(imagename=file, region=__regions[index], logfile=\'" + __log_name + ".log\', dooff=True)")
        self.__commands.append("\t\t__x_positions[index].append(__gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][0])")
        self.__commands.append("\t\t__y_positions[index].append(__gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][1])")
        _Core._return_data(self, __return_file_name)
        _Core._return_reference(self)

    def _get_residual_positions(self, __return_file_name, __get_reference):
        self.__commands.append("print(\'Performing residual Gaussian fitting...\')")
        self.__commands.extend(["__x_positions = []", "__y_positions = []"])
        self.__commands.append("for file in __input:")
        self.__commands.append("\tif not file.upper().endswith(\'.FITS\'):")
        self.__commands.append("\t\tprint(\'The file %s does not contain a valid FITS extension. Skipping...\' % (file))")
        self.__commands.append("\t\tcontinue")
        __log_name = self._output + "aligned_coordinates"
        __residual_region = self._regions[int(self._residual_region)] if isinstance(self._residual_region, int) else residual_region
        self.__commands.append("\t__gauss_fit = imfit(imagename=file, region=\'" + __residual_region + "\', logfile=\'" + __log_name + ".log\', dooff=True)")
        self.__commands.append("\t__x_positions.append(__gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][0])")
        self.__commands.append("\t__y_positions.append(__gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][1])")
        _Core._return_data(self, __return_file_name)

    def _return_data(self, __return_file_name):
        self.__commands.extend(["__positions = []", "__positions.append(__x_positions)", "__positions.append(__y_positions)"])
        self.__commands.append("with open(\'" + self._output + self._CASA_DATA_FILES_PREFIX + __return_file_name +".data\', \'ab\') as __positions_file:")
        self.__commands.append("\tpickle.dump(__positions, __positions_file)")

    def _return_reference(self):
        self.__commands.extend(["__reference_positions = []", "__reference_positions.append(__reference_x)", "__reference_positions.append(__reference_y)"])
        self.__commands.append("with open(\'"+ self._output + self._CASA_DATA_FILES_PREFIX + "references.data\', \'ab\') as __reference_file:")
        self.__commands.append("\tpickle.dump(__reference_positions, __reference_file)")

    ##### EXPORT INDIVIDUAL RESIDUALS TO FILE? PER IMAGE? PER SET? ######
    def _get_residuals(self):
        __min_residual_x = __min_residual_y = 1e9
        __max_residual_x = __max_residual_y =0.0
        __sum_residual_x = __sum_residual_y = 0.0

        with open(self._output + self._CASA_DATA_FILES_PREFIX + "aligned_positions.data", 'rb') as __aligned_positions_file:
            __aligned_positions = pickle.load(__aligned_positions_file, fix_imports=True, encoding="bytes")

        for index in range(len(self._aligned_files)):
            __residual_x = self._reference_position[0] - __aligned_positions[0][index] # _reference_position set during _perform_shift
            __residual_y = self._reference_position[1] - __aligned_positions[1][index]
            __sum_residual_x += abs(__residual_x)
            __sum_residual_y += abs(__residual_y)
            __min_residual_x = abs(__residual_x) if abs(__residual_x) < abs(__min_residual_x) else __min_residual_x
            __min_residual_y = abs(__residual_y) if abs(__residual_y) < abs(__min_residual_y) else __min_residual_y
            __max_residual_x = abs(__residual_x) if abs(__residual_x) > abs(__max_residual_x) else __max_residual_x
            __max_residual_y = abs(__residual_y) if abs(__residual_y) > abs(__max_residual_y) else __max_residual_y

        __mean_residual_x = __sum_residual_x/len(__aligned_positions[0])
        __mean_residual_y = __sum_residual_y/len(__aligned_positions[1])
        print('Overall minimum residuals ({:.3g}:'.format(__min_residual_x) + ', {:.3g})'.format(__min_residual_y))
        print('Overall maximum residuals ({:.3g}:'.format(__max_residual_x) + ', {:.3g})'.format(__max_residual_y))
        print('Overall mean residuals ({:.3g}:'.format(__mean_residual_x) + ', {:.3g})'.format(__mean_residual_y))

    #def _track_residuals(self):
    #    self.__mean_residual_x.append(self.__mean_offset_x)
    #    self.__mean_residual_y.append(self.__mean_offset_y)
    #    self.__mean_weighting.append(len(self.__input[index]))
    #    self.__min_residual_x = self.__min_offset_x if abs(self.__min_offset_x) < abs(self.__min_residual_x) else self.__min_residual_x
    #    self.__min_residual_y = self.__min_offset_x if abs(self.__min_offset_x) < abs(self.__min_residual_y) else self.__min_residual_y
    #    self.__max_residual_x = self.__max_offset_x if abs(self.__max_offset_x) > abs(self.__max_residual_x) else self.__max_residual_x
    #    self.__max_residual_y = self.__max_offset_y if abs(self.__max_offset_y) > abs(self.__max_residual_y) else self.__max_residual_y

    def _perform_shift(self):
        ####### TODO: CHECK THESE FILES EXIST; OUTPUT OFFSET RESULTS TO FILE #######
        print("Aligning images...")
        with open(self._output + self._CASA_DATA_FILES_PREFIX + "positions.data", 'rb') as __positions_file:
            __positions = pickle.load(__positions_file, fix_imports=True, encoding="bytes")

        with open(self._output + self._CASA_DATA_FILES_PREFIX + "references.data", 'rb') as __reference_file:
            self._reference_position = pickle.load(__reference_file, fix_imports=True, encoding="bytes")

        for index in range(len(self._input)):
            __sum_offset_x = __sum_offset_y = 0.0
            __min_offset_x = __min_offset_y = 1e9
            __max_offset_x = __max_offset_y = 0.0
            __mean_offset_x = __mean_offset_y = 0.0

            for file_index in range(len(self._input[index])):
                __offset_x = self._reference_position[0] - __positions[0][index][file_index]
                __offset_y = self._reference_position[1] - __positions[1][index][file_index]
                __sum_offset_x += abs(__offset_x)
                __sum_offset_y += abs(__offset_y)
                __min_offset_x = abs(__offset_x) if abs(__offset_x) < abs(__min_offset_x) else __min_offset_x
                __min_offset_y = abs(__offset_y) if abs(__offset_y) < abs(__min_offset_y) else __min_offset_y
                __max_offset_x = abs(__offset_x) if abs(__offset_x) > abs(__max_offset_x) else __max_offset_x
                __max_offset_y = abs(__offset_y) if abs(__offset_y) > abs(__max_offset_y) else __max_offset_y

                __image_data, __image_hdr = fits.getdata(self._input[index][file_index], 0, header=True)
                np.nan_to_num(__image_data, copy=False) # Replace any NaNs so we can interpolate

                __degenerate_axis = __image_data.ndim - 2 # See how many degenerate axis we have so we can add them back later
                __squeezed_image = __image_data.squeeze() # We need to squeeze any degenerate axis. This shouldn't be a problem in general but may in specific cases.
                __aligned_image = np.zeros(__squeezed_image.shape)
                scipy.ndimage.interpolation.shift(__squeezed_image, np.array([__offset_y, __offset_x]), __aligned_image) # Note the x,y inversion

                for i in range(__degenerate_axis): # Add back any degenerate axis
                    __aligned_image = np.expand_dims(__aligned_image, 0)

                __extension_index = self._input[index][file_index].rfind('.')
                __directory_index = self._input[index][file_index].rfind('/')
                __output_file = self._output + self._input[index][file_index][__directory_index+1:__extension_index] + '_aligned' + self._input[index][file_index][__extension_index:]
                self._aligned_files.append(__output_file)
                fits.writeto(filename=__output_file, data=__aligned_image, header=__image_hdr, overwrite=True)
                
            __mean_offset_x = __sum_offset_x/len(self._input[index])
            __mean_offset_y = __sum_offset_y/len(self._input[index])

            __directory_index = self._input[index][file_index].rfind('/') 
            __file_descriptor = (self._input[index][file_index][:__directory_index+1] + "*.fits") if len(self._input[index]) > 1 else self._input[index][0]
            print('Minimum absolute offset for ' + __file_descriptor + ' ({:.3g}'.format(__min_offset_x) + ', {:.3g})'.format(__min_offset_y))
            print('Maximum absolute offset for ' + __file_descriptor + ' ({:.3g}'.format(__max_offset_x) + ', {:.3g})'.format(__max_offset_y))
            print('Mean absolute offset for ' + __file_descriptor + ' ({:.3g}'.format(__mean_offset_x) + ', {:.3g})'.format(__mean_offset_y))


    def _casa_write(self):
        if os.path.isfile(self._CASA_COMMANDS_FILE_NAME) and self._overwrite != True:
            raise FileExistsError('The temporary CASA command file %s already exists (overwrite_files=False).' % (self._CASA_COMMANDS_FILE_NAME))
        self.__commands_file = open(self._CASA_COMMANDS_FILE_NAME, "w")
        self.__commands_file.write(self._CASA_HEADER + "\n")
        for command in self.__commands:
            self.__commands_file.write("\n" + command)
        self.__commands_file.write("\n")
        self.__commands_file.close()

