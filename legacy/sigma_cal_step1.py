# ~ min
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.stats import norm
from colorama import Fore, Back, Style
from glob import glob
import os

## TO FIX
# hard-coded paths when inputing and outputing files
# repeated patterns: 
#   1. the similar code structure for row 44-48 and 91-96
#.  2. the same code for row 62-73 and 78-89
# plotting mixed with computation
# one single big for loop

cubes = glob('../4rotateback/revrot*')
hdu_mask = fits.open('../0masks/3Dmask.fits')
mask = hdu_mask[0].data
for f in np.arange(len(cubes)):
    print(Back.GREEN)
    print(Style.RESET_ALL)
    hdu_cube = fits.open(cubes[f])
    obsdate = os.path.basename(cubes[f])[7:-5]
    print(obsdate)
    cube = hdu_cube[0].data
    hdu_var = fits.open('../5remove_stripes/var_stripes_masked_' + obsdate + '.fits')
    var = hdu_var[0].data
    hdu_stripes = fits.open('../0masks/stripe_masks/stripe_mask_'+obsdate+'.fits')
    stripes = hdu_stripes[0].data
    # loop over layers ------------------------------------
    MU = []
    RATIO = [] # mu/(sig/√Npix)
    Scaling = []
    sigma_before = []
    sigma_after = []
    for i in np.arange(len(cube[:,0,0])):
        maskpix = (mask[i,:,:]==0) & (stripes==0)
        sky = cube[i,:,:][maskpix]
        nanpix = np.isnan(sky)==False
        sky = sky[nanpix]
        skyvar = var[i,:,:][maskpix]
        skyvar = skyvar[nanpix]
        SNR_all = sky / np.sqrt(skyvar)
        SNR_trial = SNR_all[abs(SNR_all) < 3]
        if len(SNR_trial) == 0:
            MU.append(np.nan)
            RATIO.append(np.nan)
            Scaling.append(np.nan)
            sigma_before.append(np.nan)
            sigma_after.append(np.nan)
            continue
        mu, sig = norm.fit(SNR_trial)
        sigma_before.append(sig)
        print(f,'th cube:', i, 'trial 1:',mu,sig)
        if round(sig,1) != 1.:
            for l in np.arange(10):
                SNR_trial = SNR_all[(SNR_all < mu+3*sig) & (SNR_all > mu-3*sig)]
                mu1, sig1 = norm.fit(SNR_trial)
                scaling = sig1**2
                print(f,'th cube:', i, 'trial',l+2,':','mu:',mu1,'sig:',sig1,'scaling:',scaling)
                if round(sig1,3) == round(sig,3):
                    var[i,:,:] = var[i,:,:] * scaling
                    Scaling.append(scaling)
                    #sky_fin = cube[i,:,:][mask[i,:,:]==0 and stripes==0]
                    #sky_fin = sky_fin[np.isnan(sky_fin)==False]
                    skyvar_fin = var[i,:,:][maskpix]
                    skyvar_fin = skyvar_fin[nanpix]
                    SNR_fin = sky / np.sqrt(skyvar_fin)
                    SNR_new = SNR_fin[abs(SNR_fin) < 3]
                    mu_after, sig_after = norm.fit(SNR_new)
                    sigma_after.append(sig_after)
                    MU.append(mu_after)
                    RATIO.append(mu_after / sig_after / np.sqrt(len(SNR_new)))
                    break
                sig = sig1
                mu = mu1
                if l == 9:
                    var[i,:,:] = var[i,:,:] * scaling
                    Scaling.append(scaling)
                    #sky_fin = cube[i,:,:][mask[i,:,:]==0 and stripes==0]
                    #sky_fin = sky_fin[np.isnan(sky_fin)==False]
                    skyvar_fin = var[i,:,:][maskpix]
                    skyvar_fin = skyvar_fin[nanpix]
                    SNR_fin = sky / np.sqrt(skyvar_fin)
                    SNR_new = SNR_fin[abs(SNR_fin) < 3]
                    mu_after, sig_after = norm.fit(SNR_new)
                    sigma_after.append(sig_after)
                    MU.append(mu_after)
                    RATIO.append(mu_after / sig_after / np.sqrt(len(SNR_new)))
                    print(f,'th cube:', i, 'not converged','mu:',mu,'sig:',sig,'scaling:',scaling)
        else:
            print(f,'th cube:', i, 'no calibration:', mu,sig)
            MU.append(mu)
            RATIO.append(mu / sig / np.sqrt(len(SNR_trial)))
            scaling = sig**2
            var[i,:,:] = var[i,:,:] * scaling
            Scaling.append(scaling)
            sigma_after.append(sig)
    # save the var cube ------------------------------------
    hdu_var_cal = fits.PrimaryHDU(data=var,header=hdu_var[0].header)
    hdu_var_cal.writeto('varcube_cal/varcube_step1_' + obsdate + '.fits', overwrite=True)
    np.savetxt('sigma_before/sigma_before_' + obsdate + '.csv', sigma_before, delimiter=',')
    np.savetxt('sigma_after/sigma_after_' + obsdate + '.csv', sigma_after, delimiter=',')
    np.savetxt('scaling/scaling_' + obsdate + '.csv', Scaling, delimiter=',')
    np.savetxt('mu/mu_' + obsdate + '.csv', MU, delimiter=',')
    np.savetxt('ratio/mu_err_ratio_' + obsdate + '.csv', RATIO, delimiter=',')
    # plots --------------------------------------------------
    wave = np.arange(4600,4600+cube.shape[0])
    plt.figure(figsize=(10,5))
    plt.scatter(wave,sigma_before,s=2)
    plt.axhline(y=1.0,ls='--',c='r')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\sigma$ before calibration')
    plt.savefig('sigma_before/sigma_before_' + obsdate + '.png')
    plt.close()
    plt.figure(figsize=(10,5))
    plt.scatter(wave,sigma_after,s=2)
    plt.axhline(y=1.0,ls='--',c='r')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\sigma$ after step 1 calibration')
    plt.savefig('sigma_after/sigma_after_' + obsdate + '.png')
    plt.close()
    plt.figure(figsize=(10,5))
    plt.scatter(wave,Scaling,s=2)
    plt.axhline(y=1.0,ls='--',c='r')
    plt.xlabel(r'$\lambda$')
    plt.ylabel('scaling factor')
    plt.savefig('scaling/scaling_' + obsdate + '.png')
    plt.close()
    plt.figure(figsize=(10,5))
    plt.scatter(wave,MU,s=2)
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\mu$')
    plt.savefig('mu/mu_' + obsdate + '.png')
    plt.close()
    plt.figure(figsize=(10,5))
    plt.scatter(wave,RATIO,s=2)
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\mu$/($\sigma$/√Npix)')
    plt.savefig('ratio/ratio_' + obsdate + '.png')
    plt.close()
    #exit()
