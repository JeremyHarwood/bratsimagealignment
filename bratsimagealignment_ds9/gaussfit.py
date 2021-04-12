# take a fits image, define a sub-region with ds9 regions, fit a
# beam-shaped Gaussian, and return

from __future__ import print_function
from astropy.io import fits
from astropy.wcs import WCS
import pyregion
import numpy as np
import numexpr as ne
from scipy.optimize import least_squares

gfactor=2.0*np.sqrt(2.0*np.log(2.0))

def gaussian(xsize,ysize,x0,y0,sx,sy,pa):
    X, Y = np.meshgrid(np.arange(0,xsize,1.0), np.arange(0,ysize,1.0))
    pa*=np.pi/180.0
    a=0.5*((np.cos(pa)/sx)**2.0+(np.sin(pa)/sy)**2.0)
    b=0.25*((-np.sin(2*pa)/sx**2.0)+(np.sin(2*pa)/sy**2.0))
    c=0.5*((np.sin(pa)/sx)**2.0+(np.cos(pa)/sy)**2.0)
    
    return ne.evaluate('exp(-(a*(X-x0)**2.0+2*b*(X-x0)*(Y-y0)+c*(Y-y0)**2.0))')

class GFit(object):

    def __init__(self,fitsfile,region_name):
        # Read in the files and store the subimage
        hdu=fits.open(fitsfile)
        r=pyregion.open(region_name).as_imagecoord(hdu[0].header)

        mask=r.get_mask(hdu=hdu[0])

        # get bounding box
        xc=np.arange(hdu[0].data.shape[1])
        yc=np.arange(hdu[0].data.shape[0])

        xv, yv = np.meshgrid(xc, yc, sparse=False, indexing='xy')

        xmin=np.min(xv[mask])
        xmax=np.max(xv[mask])+1
        ymin=np.min(yv[mask])
        ymax=np.max(yv[mask])+1

        self.subim=hdu[0].data[ymin:ymax,xmin:xmax]
        self.smask=mask[ymin:ymax,xmin:xmax]

        self.ysize,self.xsize=self.subim.shape

        pixscale=hdu[0].header['CDELT2']
        self.bmaj=hdu[0].header['BMAJ']/pixscale/gfactor
        self.bmin=hdu[0].header['BMIN']/pixscale/gfactor
        self.bpa=hdu[0].header['BPA']
        self.xmin=xmin
        self.ymin=ymin
        self.hdu=hdu
        self.wcs=WCS(self.hdu[0])

    def feval(self,norm,xp,yp):
        return norm*gaussian(self.xsize,self.ysize,xp,yp,self.bmin,self.bmaj,self.bpa)

    def resid(self,norm,xp,yp):
        return ((self.subim-self.feval(norm,xp,yp))**2.0)[self.smask]

    def chi2(self,norm,xp,yp):
        return np.sum(self.resid(norm,xp,yp))

    def sresid(self,p):
        # for passing to scipy
        return self.resid(p[0],p[1],p[2])

    def fit(self):
        # initial guess is centre of region
        # return norm, xpos, ypos, ra, dec
        xp=self.xsize/2.0
        yp=self.ysize/2.0
        norm=np.max(self.subim)
        self.result=least_squares(f.sresid,(norm,xp,yp),method='lm')
        xpos=self.result.x[1]+self.xmin
        ypos=self.result.x[2]+self.ymin
        wpos=self.wcs.wcs_pix2world(xpos,ypos,0)
        ra=float(wpos[0])
        dec=float(wpos[1])
        return [self.result.x[0],xpos,ypos,ra,dec]

if __name__=='__main__':
    f=GFit('output.fits','ds9.reg')
    print(f.fit())

