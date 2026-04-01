#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# barbidule! pour "table ronde optique astro-fr/discord mars 2026"
# exemple de script de calcul de dOTF à partir de deux series d'images brutes

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import cv2


def shift_add_pond(im, count, i_mx):
    
    tmp = np.zeros((nbpx,nbpx), dtype=np.float64)
    j = 0.
    for i, fn in enumerate(im):
        
        img = np.fromfile(fn, count = count, dtype=img_type) * 1.
        img2 = np.reshape(img,(w,h))
        arg_max = np.argmax(img2)
        maxi = img2.flat[int(arg_max)] * 1.
    
        if maxi < i_mx:
            
            j += 1.
            pos = np.unravel_index(int(arg_max), img2.shape)
            tmp += (img2[pos[0]-nbpx//2:pos[0]+nbpx//2,
                          pos[1]-nbpx//2:pos[1]+nbpx//2] * maxi / i_mx)
    
    print(i, j, np.round(j * 100. / (i*1.), 2))
    if j!=0:
        tmp /= float(j)
    
    return tmp


if __name__ == '__main__':
    
    base = '/home/voyageur/dotf/' # repertoire contenant les images
    
    lmbd = 650 # longueur d'onde de la source
    fnm = "essai12b.wft" # nom du fichier .wft
    img_type = np.uint16    # daA3840-45um 4k raw 12 bits / 2 bytes
    bpp = 2  # bytes per pixel
    i_mx = 4095.  # pixel saturation
    im1 = glob(base+'Basler_daA3840-45um__40702814__20260215_191526432_*.raw')
    im2 = glob(base+'Basler_daA3840-45um__40702814__20260215_191317774_*.raw')
    
    # image size roi, swap h and w <-- np.fromfile & raw image
    w = 1000
    h = 1000
    count = w * h * bpp

    nbpx = 256  # puissance de 2, taille de ROI avant calcul OTF
    chkbd = -2. * (np.indices((nbpx,nbpx)).sum(axis=0) % 2 - 0.5)
    
    psf1 = shift_add_pond(im1, count, i_mx) # empilement des images
    psf2 = shift_add_pond(im2, count, i_mx) # avec ponderation par le flux
        
    psf1 -= np.median(psf1)
    psf2 -= np.median(psf2)

    otf1 = np.fft.ifft2(psf1 * chkbd) * chkbd / np.sum(psf1)
    otf2 = np.fft.ifft2(psf2 * chkbd) * chkbd / np.sum(psf2)
   
    tmp = np.log(np.abs(otf1) + 1.)
    tmp[nbpx//2,nbpx//2] = np.min(tmp)
    # iok = tmp > np.sqrt(np.mean(tmp**2))
    iok = tmp > np.mean(tmp)
    plt.matshow(iok.astype(np.float16), fignum=1)
    
    k = np.min(np.abs(otf1[iok]/otf2[iok]))
    otf2 *= k
    
    dotf = otf1 - otf2
    phase = np.angle(dotf)
    plt.matshow(phase, fignum=2)
    
    ph = np.log(np.abs(np.imag(dotf)/np.max(np.imag(dotf)))+1.)
    windowName = 'Selection wft ("enter" si ok)'
    roi = np.asarray(cv2.selectROI(windowName, ph, False, False))
    cv2.destroyWindow(windowName)
        
    crop = phase.copy()[roi[1]:roi[1]+roi[3],roi[0]:roi[0]+roi[2]]
    crp_shp = crop.shape
    mn_cs = np.min(crop.shape)
    crop = crop[crp_shp[0]//2-mn_cs//2:crp_shp[0]//2+mn_cs//2,
                crp_shp[1]//2-mn_cs//2:crp_shp[1]//2+mn_cs//2]
    ell = str(mn_cs - 0.5)
    
    plt.matshow(crop, fignum=3)
    
    crop *= lmbd / (2. * np.pi * 1000.)
    
    hdr = ' '.join(map(str, np.sort(crop.shape)))
    np.savetxt(fnm, np.transpose(crop.flatten()),
               fmt="%f", comments=' ',header=hdr)
    f = open(fnm, "a")
    f.write('ellipse '+ell+' '+ell+' '+ell+' '+ell+'\n' )
    f.write("DIAM 135\n")
    f.write('ROC 1440\n')    
    f.write('lambda '+str(lmbd)+'\n')
    f.close()
    
