import sys, string, os, math, shutil, time
import numpy as np
import glob
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import shgo
from scipy.optimize import minimize
from scipy import optimize
from scipy import special
from scipy.stats import norm
from scipy.stats import chisquare
import matplotlib.pyplot as plt
import scipy.integrate as integrate

##   Procedure to convolve a spectrum with instrumental, rotational,
##   and/or macroturbulence profiles from FASTWIND code
##
##   Routine based on convol.pro from S. Simon-Diaz (iacob-broad)
##   and Fastwind/idlproc/convol.pro from J. Puls.
#    CONVOL,sx,sy,cx,cy,px,py,resol=resol,vsini=vsini,beta=beta, '
#    sx,sy  -- Espectro 
#    cx,cy  -- Espectro convolucionado
#    px,py  -- Perfil de convolucion
#    resol  -- spectral resolving power (perfil instrumental)
#    vsini  -- Velocidad de rotacion proyectada (km/s)
#    beta   -- limb darkening coeficient (default: 1.5 --> epsilon=0.6)
#    vmacro -- Macroturbulencia radial (km/s)
#    res0   -- Minima resolucion en x en A
#--------------------------------------------------------

def conv_spec(sx,sy,px,py):
    n=len(sx);
    if (n != len(px)): sys.exit('ERROR: Different array dimensions');  
    shft= int(n/2.0)
    tty= np.zeros(n) + (sy[0] + sy[-1])/2.0
    sy_t= sy - tty
    cx = sx.real
    dxx = cx[-1] - cx[0]
    dxn= dxx/float(n-1)
    cy= np.roll(np.fft.ifft(np.fft.fft(sy_t,norm="forward")*np.fft.fft(py,norm="forward"),norm="forward"),shft) 
    cy= dxn*float(n)*cy.real
    cy= cy + tty     
    return cx,cy

def convol(sx,sy,vsini,vmacro,vgauss,resol):
    cc          = 299792.458                # VELOCIDAD DE LA LUZ (km/s)
    #res0        = (sx[3]-sx[2])*0.501
    res0        = (sx[1] - sx[0])*0.1
    rsmpl_fac   = 0.7
    sxx         = sx               # ESPECTRO DE ENTRADA (x)
    syy         = sy               # ESPECTRO DE ENTRADA (y)
    n           = len(sxx) 
    nn          = n
    rd          = (sxx[-1]-sxx[0])/(n-1.)   # DISTANCIA ENTRE PUNTOS
    sxm         = (sxx[-1]+sxx[0])/(2.)     # PUNTO MEDIO
    rdd         = rd                                 
    intp        = 0.0
#--------- Identificar vectores x no equisdistantes ---------#
    dsx     = abs(sxx[1:] - sxx[0:nn-1])
    rd_min  = min(dsx)
    rd_max  = max(dsx)
    eps     = (2.0e-03)/5000.
    diff    = abs(rd_max - rd_min)
    not_eqdist = 1 if (diff > (eps*sxm)) else 0
    #------------ Max(DLAM)-----------------
    if (not_eqdist == 1):
        rdd = res0 if (res0 !=0) else max([rd_min*rsmpl_fac,2.7*eps*sxm])
        nn = (sxx[-1]-sxx[0])/rdd
        intp = 1;
    #---------------------------------------    
    xpow = np.log(nn)/np.log(2.0)    
    pow = int(xpow)
    #--------------- if pow-----------------    
    if ((2 < pow <= 15) or (intp ==1)) and ((xpow-pow) != 0.0):
        pow = pow+1
        nn_x = int(2.0**pow)
        if ((float(nn)/float(nn_x)) > rsmpl_fac) and (not_eqdist == 0): 
            nn = nn_x*2
        else:
            nn = nn_x
        #-----------end if pow---------------
        rdd = (sxx[-1] - sxx[0])/(nn-1)
        sxx = np.arange(nn)*rdd + sx[0]
        syy = interp1d(sx,sy,bounds_error=False)(sxx)
        intp = 1
        #       
    px = (np.arange(nn) - int(nn/2))*rdd
    
# *********************************************************
# ---------------rotation profile -------------------------
# *********************************************************
    beta = 1.5
    if vsini>0:
        eps = beta/(1.0+beta)
        lambda0 = (sxx[-1]+sxx[0])/2.0
        normf = cc/(vsini*lambda0)
        xi = normf*px
        index = 0
        while index < len(xi):
            if (abs(xi[index])>1): xi[index]=1.0
            index +=1
        xi_2 = xi*xi
        py = (2.*np.sqrt(1.- xi_2)/np.pi + (1.-xi_2)*beta/2.)*normf /(1.+ (2.*beta/3.))
        cx,cy=conv_spec(sxx,syy,px,py)
        sxx = cx
        syy = cy 

# *********************************************************
# -----radial-tangencial macroturbulence profile ----------
# *********************************************************
    if vmacro>0:
        aat = 1  
        spi = np.sqrt(np.pi)
        lambda0 = (sxx[-1] + sxx[0])/2.
        aat = 1
        mt = vmacro * lambda0 / cc
        cct = 2.*aat/spi/mt
        ccr = 0
        pxmt = abs(px)/mt
        py = cct * (np.exp(-pxmt**2) + spi * pxmt *(special.erf(pxmt) - 1.0))
        cx,cy = conv_spec(sxx,syy,px,py)   
        sxx = cx
        syy = cy
#; ************************************
#; ----- Instrumental profile -------
#; ************************************
    if resol>0:
        lambda0 = (sxx[-1] + sxx[0])/2
        fwhm    = lambda0/resol
        width   = fwhm/2.0/np.sqrt(np.log(2.0))
        sigma   = width
        sqrt_pi = np.sqrt(np.pi)
        py      = np.exp(-(px/sigma)**2)/(sqrt_pi*sigma)
        cx,cy   = conv_spec(sxx,syy,px,py)
        sxx = cx
        syy = cy
    return sxx,syy 

# *********************************************************
# -------------- Gaussian profile ----------------
# *********************************************************    
    if vgauss>0:
        lambda0 = (sxx[-1]+sxx[0])/2.
        sigma = vgauss*lambda0/cc
        sqrt_pi = np.sqrt(np.pi)
        pvg0 = vgauss*lambda0/cc
        py   = np.exp(-((px-pvg0)/sigma**2))/(sqrt_pi*sigma)
        cx,cy = conv_spec(sxx,syy,px,py)   
        sxx = cx
        syy = cy
                      