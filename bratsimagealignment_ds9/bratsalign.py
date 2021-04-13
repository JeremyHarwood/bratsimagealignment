#!/usr/bin/env python3 
# -*- coding: utf-8 -*-

# Created by Jeremy. J. Harwood (2019)
# Contact: Jeremy.Harwood@physics.org or J.Harwood3@herts.ac.uk

# Please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803 if you have made use of this class.

# Automatically fits Gaussians and aligns images in pixel space. Particularly useful for spectral index and spectral age fitting on resolved scales.

####################################################################################################################################################

import os, glob, scipy.ndimage, sys, pyregion
import numexpr as ne
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

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
    _DEFAULT_OUTPUT_PATH = "./"
    _GFACTOR=2.0*np.sqrt(2.0*np.log(2.0))


class Setup(_Constants):

    def __init__(self, input_files:list, region_files:list, output_dir=_Constants._DEFAULT_OUTPUT_PATH, mode=2, reference_image=0, reference_location=[1024.0, 1024.0], residual_region=0, overwrite_files=False):
        _Setters._set_input(self, input_files)
        _Setters._set_regions(self, region_files)
        _Setters._set_output(self, output_dir)
        _Setters._set_mode(self,mode)
        _Setters._set_reference_image(self,reference_image)
        _Setters._set_reference_location(self,reference_location)
        self._residual_region=residual_region # THIS NEEDS A SETTER!
        self._overwrite = overwrite_files

    def set_mode(self, mode:int):
        _Setters._set_mode(self, mode)

    def set_reference_image(self, reference_image:int):
        _Setters._set_reference_image(self, reference_image)

    def set_reference_location(self, reference_location:int):
        _Setters._set_reference_location(self, reference_location)

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
            elif not __region_files[0].upper().endswith('.REG'):
                raise ValueError('The region file must contain be in a valid DS9 style format (.reg).')
            self._regions = __region_files * len(self._input)
        elif len(__region_files) == len(self._input):
            print("Here")
            for file in __region_files:
               if not file.strip():
                   raise ValueError('The region files list contains a blank element. All elements must contain a valid value.')
               elif not os.path.isfile(file):
                   raise FileNotFoundError('The region file %s does not exist. Please ensure all paths are correct.' % (file))
               elif not file.upper().endswith('.REG'):
                   raise ValueError('The region files list must contain only DS9 style regions (.reg). Incorrect file format for: ' % (file))
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


class _Core(Setup):

    def _align (self):
        self._aligned_files = []
        self._x_residual_pos = []
        self._y_residual_pos = []
        _Core._check_mode(self)

        # Do alignment
        _Core._get_reference_position(self)
        _Core._get_positions(self)
        _Core._perform_shift(self)

        # Check alignment
        _Core._get_residual_positions(self)
        _Core._get_residuals(self)


    def _check_mode(self):
        if self._mode == 0 and self._reference_image == -1:
            raise IndexError('A reference image of -1 cannot be used for mode 1 (predefined image)') 

    def _get_reference_position(self):
        if self._mode == 0:
            print('Using a predefined image for the reference coordinates')
            __reference_norm, __reference_x, __reference_y, __reference_ra, __reference_dec = GFit(self._input[self._reference_image][0], self._regions[self._reference_image])._fit()
        elif self._mode == 1:
            print('Using a fixed reference of {:.2f}'.format(_reference_location[0]) + ', {:.2f}'.format(_reference_location[1]) + " for alignment")
            __reference_x, __reference_y = _reference_location
        elif self._mode == 2:
            print('Using the mean value of all images for the reference coordinates')
            __sum_x = 0.0
            __sum_y = 0.0
            __sum_images = 0

            if (self.set_reference_image == -1):
                for index in range(len(self._input)):
                    __sum_images+=len(fits_location)
                    for file in self._input[index]:
                        __reference_norm, __reference_x, __reference_y, __reference_ra, __reference_dec = GFit(file, self._regions[index])._fit()
                        __sum_x += __reference_x
                        __sum_y += __reference_y
            else:
                __sum_images+=len(self._input[self._reference_image])
                for file in self._input[self._reference_image]:
                    __reference_norm, __reference_x, __reference_y, __reference_ra, __reference_dec = GFit(file, self._regions[self._reference_image])._fit()
                    __sum_x += __reference_x
                    __sum_y += __reference_y

            self._reference_x = (__sum_x/__sum_images)
            self._reference_y = (__sum_y/__sum_images)
    
        print('Aligning to reference coordinates: {:.2f}'.format(self._reference_x) + ', {:.2f}'.format(self._reference_y))

    def _get_positions(self):
        print('Performing initial Gaussian fitting...')
        self._x_positions = []
        self._y_positions = []
        for index in range(len(self._input)):
            self._x_positions.append([])
            self._y_positions.append([])
            for file in self._input[index]:
                if not file.upper().endswith('.FITS'):
                    print('The file %s does not contain a valid FITS extension. Skipping...' % (file))
                    continue

                __norm, __x_pos, __y_pos, __ra, __dec = GFit(file, self._regions[index])._fit()
                self._x_positions[index].append(__x_pos)
                self._y_positions[index].append(__y_pos)

    def _get_residual_positions(self):
        print('Performing residual Gaussian fitting...')

        for file in self._aligned_files:
            if not file.upper().endswith('.FITS'):
                print('The file %s does not contain a valid FITS extension. Skipping...' % (file))
                continue

            __residual_region = self._regions[int(self._residual_region)] if isinstance(self._residual_region, int) else residual_region
            __norm, __x_pos, __y_pos, __ra, __dec = GFit(file, __residual_region)._fit()
            self._x_residual_pos.append(__x_pos)
            self._y_residual_pos.append(__y_pos)

    ##### EXPORT INDIVIDUAL RESIDUALS TO FILE? PER IMAGE? PER SET? ######
    def _get_residuals(self):
        __min_residual_x = __min_residual_y = 1e9
        __max_residual_x = __max_residual_y =0.0
        __sum_residual_x = __sum_residual_y = 0.0

        for index in range(len(self._aligned_files)):
            __residual_x = self._reference_x - self._x_residual_pos[index] # _reference_position set during _perform_shift
            __residual_y = self._reference_y - self._y_residual_pos[index]
            __sum_residual_x += abs(__residual_x)
            __sum_residual_y += abs(__residual_y)
            __min_residual_x = abs(__residual_x) if abs(__residual_x) < abs(__min_residual_x) else __min_residual_x
            __min_residual_y = abs(__residual_y) if abs(__residual_y) < abs(__min_residual_y) else __min_residual_y
            __max_residual_x = abs(__residual_x) if abs(__residual_x) > abs(__max_residual_x) else __max_residual_x
            __max_residual_y = abs(__residual_y) if abs(__residual_y) > abs(__max_residual_y) else __max_residual_y

        __mean_residual_x = __sum_residual_x/len(self._x_residual_pos)
        __mean_residual_y = __sum_residual_y/len(self._y_residual_pos)
        print('Overall minimum residuals ({:.3g}:'.format(__min_residual_x) + ', {:.3g})'.format(__min_residual_y))
        print('Overall maximum residuals ({:.3g}:'.format(__max_residual_x) + ', {:.3g})'.format(__max_residual_y))
        print('Overall mean residuals ({:.3g}:'.format(__mean_residual_x) + ', {:.3g})'.format(__mean_residual_y))

    def _perform_shift(self):
        ####### TODO: CHECK THESE FILES EXIST; OUTPUT OFFSET RESULTS TO FILE #######
        print('Aligning images...')

        for index in range(len(self._input)):
            __sum_offset_x = __sum_offset_y = 0.0
            __min_offset_x = __min_offset_y = 1e9
            __max_offset_x = __max_offset_y = 0.0
            __mean_offset_x = __mean_offset_y = 0.0

            for file_index in range(len(self._input[index])):
                __offset_x = self._reference_x- self._x_positions[index][file_index]
                __offset_y = self._reference_y - self._y_positions[index][file_index]
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

    def _gaussian(self, xsize, ysize, x0, y0, sx, sy, pa):
        X, Y = np.meshgrid(np.arange(0,xsize,1.0), np.arange(0,ysize,1.0))
        pa*=np.pi/180.0
        a=0.5*((np.cos(pa)/sx)**2.0+(np.sin(pa)/sy)**2.0)
        b=0.25*((-np.sin(2*pa)/sx**2.0)+(np.sin(2*pa)/sy**2.0))
        c=0.5*((np.sin(pa)/sx)**2.0+(np.cos(pa)/sy)**2.0)
    
        return ne.evaluate('exp(-(a*(X-x0)**2.0+2*b*(X-x0)*(Y-y0)+c*(Y-y0)**2.0))')


class GFit(_Core):

    def __init__(self, fitsfile, region_name):
        # Read in the files and store the subimage
        __image_data, __image_hdr = fits.getdata(fitsfile, 0, header=True)

        r=pyregion.open(region_name).as_imagecoord(__image_hdr)

        __image_data=__image_data.squeeze()

        mask=pyregion.get_mask(r, hdu=__image_data)

        # get bounding box
        xc=np.arange(__image_data.shape[1])
        yc=np.arange(__image_data.shape[0])

        xv, yv = np.meshgrid(xc, yc, sparse=False, indexing='xy')

        xmin=np.min(xv[mask])
        xmax=np.max(xv[mask])+1
        ymin=np.min(yv[mask])
        ymax=np.max(yv[mask])+1
        self.subim=__image_data[ymin:ymax,xmin:xmax]
        self.smask=mask[ymin:ymax,xmin:xmax]

        self.ysize,self.xsize=self.subim.shape
        pixscale=__image_hdr['CDELT2']
        self.bmaj=__image_hdr['BMAJ']/pixscale/self._GFACTOR
        self.bmin=__image_hdr['BMIN']/pixscale/self._GFACTOR
        self.bpa=__image_hdr['BPA']
        self.xmin=xmin
        self.ymin=ymin
        self.wcs=WCS(__image_hdr, naxis=2)
        
    def _feval(self,norm,xp,yp):
        return norm*self._gaussian(self.xsize, self.ysize, xp, yp, self.bmin, self.bmaj, self.bpa)

    def _resid(self,norm,xp,yp):
        return ((self.subim-self._feval(norm,xp,yp))**2.0)[self.smask]

    def _chi2(self,norm,xp,yp):
        return np.sum(self._resid(norm,xp,yp))

    def _sresid(self,p):
        # for passing to scipy
        return self._resid(p[0],p[1],p[2])

    def _fit(self):
        # initial guess is centre of region
        # return norm, xpos, ypos, ra, dec
        xp=self.xsize/2.0
        yp=self.ysize/2.0
        norm=np.max(self.subim)
        self.result=least_squares(self._sresid,(norm,xp,yp),method='lm')
        xpos=self.result.x[1]+self.xmin
        ypos=self.result.x[2]+self.ymin
        wpos=self.wcs.wcs_pix2world(xpos,ypos,0)
        ra=float(wpos[0])
        dec=float(wpos[1])
        return [self.result.x[0],xpos,ypos,ra,dec]

