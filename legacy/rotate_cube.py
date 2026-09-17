from astropy.io import fits
import numpy as np
from scipy import ndimage
from mpdaf.obj import Cube
import os

## TO FIX
# hard-coded paths when inputing and outputing files
# repeated patterns: the same code for cube and mask
# one single big for loop

fields = np.array(['RO0959','RO1001','SXDF-LAB1','CL1449','RO0958'])
add = 50
for f in range(len(fields)):
    print(fields[f])
    if f < 3:
        #continue
        path = '../'+fields[f]
        dict = '../'+fields[f]+'/esoreflex/reflex_tmp_products/muse/muse_scipost_1'
    if f >= 3:
        path = '/Volumes/LaCie 1/muse/'+fields[f]
        dict = '/Volumes/LaCie 1/muse/'+fields[f]+'/esoreflex/reflex_tmp_products/muse/muse_scipost_1'
    for file in os.scandir(dict):
        tmp_path = file.name
        hdul_pipe = fits.open(file.path+'/DATACUBE_FINAL.fits')
        obsdate = hdul_pipe[0].header['DATE-OBS']
        print(tmp_path,obsdate)
        if fields[f] == 'CL1449':
            hdul_cube = fits.open(path+'/noise_cal/2nd_run/sky_level_cor1/datacube_skycor_'+obsdate+'.fits')
            #hdu_mask = fits.open(path+'/noise_cal/2nd_run/masks/3Dmask.fits')
        else:
            hdul_cube = fits.open(path+'/noise_cal/3rd_run/sky_level_cor1/datacube_skycor_'+obsdate+'.fits')
            #hdu_mask = fits.open(path+'/noise_cal/3rd_run/masks/3Dmask.fits')
        posang = hdul_cube[0].header['HIERARCH ESO INS DROT POSANG']
        cube = hdul_cube[1].data
        #mask = hdu_mask[0].data
        print('original:',cube.shape)
        inser1 = np.empty((cube.shape[0],add,cube.shape[2]))
        inser1[:] = np.nan
        cube1 = np.concatenate((inser1,cube,inser1),axis=1)
        #mask1 = np.concatenate((inser1,mask,inser1),axis=1)
        inser2 = np.empty((cube1.shape[0],cube1.shape[1],add))
        inser2[:] = np.nan
        cube2 = np.concatenate((inser2,cube1,inser2),axis=2)
        #mask2 = np.concatenate((inser2,mask1,inser2),axis=2)
        print('expanded:',cube2.shape)
        cube_rot = ndimage.rotate(cube2,posang, axes=(1,2), reshape=False,cval=np.nan,order=0)
        #mask_rot = ndimage.rotate(mask2,posang, axes=(1,2), reshape=False,cval=np.nan,order=0)
        i_del = 0
        for i in range(cube_rot.shape[1]):
            j = i - i_del
            if np.all(np.isnan(cube_rot[:,j,:])) == True:
                cube_rot = np.delete(cube_rot,j,1)
                #mask_rot = np.delete(mask_rot,j,1)
                i_del += 1
        a_del = 0
        for a in range(cube_rot.shape[2]):
            b = a - a_del
            if np.all(np.isnan(cube_rot[:,:,b])) == True:
                cube_rot = np.delete(cube_rot,b,2)
                #mask_rot = np.delete(mask_rot,b,2)
                a_del += 1
        print('rotated&cut:',cube_rot.shape)
        hdu_rot = fits.PrimaryHDU(data=cube_rot,header=hdul_cube[1].header)
        hdu_rot.writeto(fields[f]+'/cube_rot_'+obsdate+'.fits',overwrite=True)
        #hdu_rot_mask = fits.PrimaryHDU(data=mask_rot,header=hdul_cube[1].header)
        #hdu_rot_mask.writeto(fields[f]+'/mask_rot/mask_rot_'+obsdate+'.fits',overwrite=True)
        rotated = Cube(fields[f]+'/cube_rot_'+obsdate+'.fits')
        ima = rotated.median(axis=0)
        ima.write(fields[f]+'/ima_rot_'+obsdate+'.fits')
        hdul_pipe.close()
        hdul_cube.close()
        #hdu_mask.close()
    print('-------------------------------------')

