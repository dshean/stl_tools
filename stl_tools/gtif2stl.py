#! /usr/bin/env python

#David Shean
#dshean@gmail.com
#2/19/14

#Utility to convert geotiff to scaled stl for printing

#Note: Inputs must have proper nodata definition
#Combined surface+bed using gdalwarp to form continuous piece for difference - makes scaling easier
#Masking in numpy2stl is inconsistent with the numpy ma and gdal nodata approaches

#For isolating pieces w/ two tif surfaces
#Could clip surface to ice extent, clip bed to ice extent, then cat ascii output - but need sides, add edge facets 

import sys
import os
import subprocess

import gdal
import numpy as np
from scipy.spatial import Delaunay
import stl_tools

import geolib
import malib

#Dimensions of printer (mm)
max_width = 152.4
max_length = 294.6
max_height = 154.9

#Vertical exaggeration of input geospatial data
vexagg = 100.0
base_thickness = 3.0 
stride = 8 
ascii=False

#Print resolution in mm
#Min thickness 4*res (4 layers)
res = 0.02

src_fn = sys.argv[1]

src_ds = gdal.Open(src_fn, gdal.GA_ReadOnly)
src_gt = src_ds.GetGeoTransform()
src_extent = geolib.get_extent(src_ds)
src_b = src_ds.GetRasterBand(1)
src_ndv = src_b.GetNoDataValue()

src_a = malib.ds_getma(src_ds)

#Smooth input
import dem_filter
#f = dem_filter.gauss_fltr_astropy(src_a, size=3, sigma=1)
f_g = dem_filter.gauss_fltr_astropy(src_a, size=17, sigma=3)
#f_m = dem_filter.median_fltr_skimage(src_a, radius=5, erode=0)
#from skimage.measure import block_reduce
#f = block_reduce(src_a, block_size=(3,3), func=np.ma.mean)

src_a = f_g

#ant_dim_m = np.array([6667000.0,6667000.0,11026.0])
ant_dim_m = np.array([6667000.0,6667000.0,11136.0])
x_dim = float(src_ds.RasterXSize*src_gt[1])
y_dim = float(src_ds.RasterYSize*-src_gt[5])
z_dim = float(src_a.ptp())

#This scales the output dimensions based on actual valid area
#Note: this will be centered, so won't necessarily align with bed
src_a_trim = malib.masktrim(src_a)
x_dim = float(src_a_trim.shape[1]*src_gt[1])
y_dim = float(src_a_trim.shape[0]*-src_gt[5])

min_dim_mm = min(max_width, max_length)
max_width = min_dim_mm * (x_dim/ant_dim_m[0])
max_length = min_dim_mm * (y_dim/ant_dim_m[1])
#This should scale appropriately with the above term
max_height = max(max_width, max_length) * vexagg * (z_dim/max(x_dim,y_dim))
#max_height = max(max_width, max_length) * (z_dim/max(x_dim,y_dim))

#Remove z offset and scale, then set units to output mm
#src_a -= src_a.min()
src_a = max_height*(src_a - src_a.min()) / src_a.ptp()

#Use the stl_tools lib for stl output
#values < mask_val are not included
if ascii:
    out_stl_fn = os.path.splitext(src_fn)[0]+'_ascii.stl'
else:
    out_stl_fn = os.path.splitext(src_fn)[0]+'.stl'

#Downsample
out_a = src_a[::stride,::stride]

#Mask properly
#mask_val = src_ndv  
#mask_val = -9999
mask_val = 0.0 
out_a = out_a.filled(mask_val)

#Scale has been lobotomized in numpy2stl
scale = 1.0 
min_thickness_percent = 0.0

stl_tools.numpy2stl(out_a, out_stl_fn, scale=scale, solid=True, mask_val=mask_val+0.00000001, max_width=max_width, max_depth=max_length, max_height=max_height, min_thickness_percent=min_thickness_percent, ascii=ascii)

sys.exit()

#The following was used to spit out xyz points from input grid, with the hope of creating mesh externally

def edgepts(a):
    a = malib.checkma(a)
    px = np.ma.notmasked_edges(a, axis=0)
    edges0 = np.ma.notmasked_edges(a, axis=0)
    edges1 = np.ma.notmasked_edges(a, axis=1)
    edges = np.array([np.concatenate([edges0[0][0], edges0[1][0], edges1[1][0], edges1[0][0]]), np.concatenate([edges0[0][1], edges0[1][1], edges1[1][1], edges1[0][1]])])
    #x = np.concatenate((px[0][1][::1], px[1][1][::-1], [px[0][1][0]]))
    #y = np.concatenate((px[0][0][::1], px[1][0][::-1], [px[0][0][0]]))
    x = edges[0]
    y = edges[1]
    return zip(x, y)

def tif2pts(src_fn, dst_fn, outsize=100):
    src_ds = gdal.Open(src_fn, gdal.GA_ReadOnly)
    src_b = src_ds.GetRasterBand(1)
    src_ndv = src_b.GetNoDataValue()
    cmd=['gdal_translate', '-a_nodata', str(src_ndv), '-of', 'XYZ', '-outsize', '%s%%' % str(outsize), '%s%%' % str(outsize), src_fn, dst_fn] 
    print ' '.join(cmd)
    subprocess.call(cmd)
    #This is much less efficient than grep
    xyz = np.genfromtxt(dst_fn)
    return xyz[xyz[:,2] != src_ndv] 

def scalepts(xyz):
    xlim = (xyz[:,0].min(), xyz[:,0].max()) 
    ylim = (xyz[:,1].min(), xyz[:,1].max()) 
    zlim = (xyz[:,2].min(), xyz[:,2].max()) 

dst_fn = os.path.splitext(src_fn)[0]+'_10p.xyz'
out_fn = os.path.splitext(src_fn)[0]+'_10p_scaled.xyz'

#Spit out points
src_xyz = tif2pts(src_fn, dst_fn, outsize=10)

#Upper left corner is origin
#Should be positive down and positive right
x_offset = src_gt[0]
y_offset = src_gt[3]
#Lowest point
zlim = malib.calcperc(src_xyz[:,2], (0.1, 99.9))
zlim = (src_xyz[:,2].min(), src_xyz[:,2].max())
#clamp to zlim
src_xyz[src_xyz[:,2] < zlim[0]] = zlim[0]
src_xyz[src_xyz[:,2] > zlim[1]] = zlim[1]
z_offset = zlim[0]

#Just the top surface
#out_xyz = np.array([(src_xyz[:,0] - x_offset)*x_scale, (src_xyz[:,1] - y_offset)*y_scale, (src_xyz[:,2] - z_offset)*z_scale]).T
#np.savetxt(out_fn, out_xyz, fmt='%0.3f %0.3F %0.3f')

#Both top and bottom
out_xyz_top = np.array([(src_xyz[:,0] - x_offset)*x_scale, (src_xyz[:,1] - y_offset)*y_scale, (src_xyz[:,2] - z_offset)*z_scale+base_thickness]).T
out_xyz_bottom = np.array([(src_xyz[:,0] - x_offset)*x_scale, (src_xyz[:,1] - y_offset)*y_scale, src_xyz[:,2]*0]).T
out_xyz = np.concatenate([out_xyz_top, out_xyz_bottom])
np.savetxt(out_fn, out_xyz, fmt='%0.3f %0.3F %0.3f')

#tri = Delaunay(out_xyz)
