import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.fft import fft, fftfreq
from scipy.optimize import shgo,minimize,curve_fit
from scipy.signal import find_peaks

def find_line_value(line_name):
    with open("iacob_broad_lines.dat", "r") as file:
        for line in file:
            if line.startswith(line_name):
                line_split = line.split()
                return float(line_split[1])
    print("Line not found in the list")
    return None

def read_asc(file_path):
    data = np.loadtxt(file_path)
    wave = data[:, 0]
    flux = data[:, 1]
    return wave, flux

def continuum(w, f, xc):
    cont_window = 50 
    window = (w > xc - cont_window/2) & (w < xc + cont_window/2)
    peaks, _ = find_peaks(f[window], prominence=0.01)
    peaks = peaks[peaks != np.argmin(f[window])]
    continuum = np.mean(f[window][peaks])
    return continuum

def ewgauss(delta, ew, lam0, resol):
        nn = len(delta)
        ycont = 1
        delta_min = np.min(delta)
        delta_max = np.max(delta)
        paso = (delta_max - delta_min) / nn
        delta = np.arange(nn) * paso + delta_min
        fwhm = lam0 / resol
        b = np.log(2) * (2 / fwhm)**2
        A = ew * np.sqrt(b / np.pi)
        return ycont - A * np.exp(-b * (delta - lam0)**2) 

def fit_gaussian(w, f , ew, xc ,resol):
   # Estimación inicial de los parámetros
    ew_0 = ew
    lam0_0 = xc
    resol_0 = resol
    popt, pcov = curve_fit(gaussian, w, f, p0=[ew_0, lam0_0, resol_0])
    return popt

def gaussian(x, A, lam0, sigma):
    return A*np.exp(-(x - lam0)**2/(2*sigma**2)) + 1

def optimal_win(w, f, xc):
    mask = (w > xc - 5) & (w < xc + 5) 
    w_fit = w[mask]
    f_fit = f[mask]
    # Parámetros iniciales para la función gaussiana
    amp = 0.025  # Amplitud fija
    x0 = xc
    sigma = np.std(w_fit) 
    popt, pcov = curve_fit(gaussian, w_fit, f_fit, p0=[amp, x0, sigma])
    #====================== Optimal window ======================#
    sigma_opt = popt[2]
    windows = []
    if sigma_opt >= 0:
        k_list = [2.0, 2.5, 3.0]
        print("Usando k_list (+) para la elaboración de ventanas")
    else:
        k_list = [-3, -2.5, -2.0, 2.0, 2.5, 3.0]
        print("Usando k_list (-) para la elaboración de ventanas")
    for k in k_list:
        w_min = xc - k * sigma_opt
        w_max = xc + k * sigma_opt
        window_mask = (w >= w_min) & (w <= w_max)
        w_window = w[window_mask]
        f_window = f[window_mask]
        if len(w_window) > 0 and len(f_window) > 0:
            windows.append((w_window, f_window))
    w_opt1 , f_opt1 = windows[0]
    w_opt2 , f_opt2 = windows[1]
    w_opt3 , f_opt3 = windows[2]
    if sigma_opt >= 0:
        w_opt,f_opt = w_opt3 , f_opt3
    else:
        w_opt,f_opt = w_opt1 , f_opt1
    return popt, w_opt, f_opt

def reduce(cx, cy, xc, porcentaje, num):
    wx = np.linspace(cx[0], cx[-1], num)  # num: Rango de valores x para la interpolación
    # Interpolación de los valores y
    fx = np.interp(wx, cx, cy)
    # Establecer punto medio de la línea
    i1 = np.argmax(fx[wx<=xc]) 
    i2 = np.argmax(fx[wx>=xc])
    i3 = np.argmin(fx)
    wx1 = wx[wx<=xc][i1]
    wx2 = wx[wx>=xc][i2]
    fx1 = fx[wx<=xc][i1] 
    fx2 = fx[wx>=xc][i2]
    d1 = fx1 - min(fx) 
    d2 = fx2 - min(fx)
    d1_new = d1*(porcentaje)/100 
    #d2_new = d2*(porcentaje)/100
    #d1_new = 0
    d2_new = 0
    fx1_new = d1_new + min(fx) 
    fx2_new = d2_new + min(fx)
    i1_new = np.argmin(np.abs(fx[:i3] - fx1_new))
    i2_new = len(fx) - np.argmin(np.flipud(np.abs(fx[i3+1:]-fx2_new)))-1
    fx1_new = fx[i1_new] 
    fx2_new = fx[i2_new]
    wx_new = wx[i1_new:i2_new +1]
    fx_new = fx[i1_new:i2_new +1]
    return wx_new , fx_new

def reducce(wx, fx, xc, porcentaje, gauss):
    i1 = np.argmax(fx[wx<=xc]) 
    i2 = np.argmax(fx[wx>=xc])
    i3 = np.argmin(fx)
    wx1 = wx[wx<=xc][i1]
    wx2 = wx[wx>=xc][i2]
    fx1 = fx[wx<=xc][i1] 
    fx2 = fx[wx>=xc][i2]
    d1 = fx1 - min(fx) 
    d2 = fx2 - min(fx)
    d1_new = d1*(porcentaje)/100 
    if gauss>0 :
        d2_new = d2*(porcentaje)/100
    else:
        d2_new = d2*(porcentaje + 10)/100
    fx1_new = d1_new + min(fx) 
    fx2_new = d2_new + min(fx)
    i1_new = np.argmin(np.abs(fx[:i3] - fx1_new))
    i2_new = len(fx) - np.argmin(np.flipud(np.abs(fx[i3+1:]-fx2_new)))-1
    fx1_new = fx[i1_new] 
    fx2_new = fx[i2_new]
    wx_new = wx[i1_new:i2_new + 1]
    fx_new = fx[i1_new:i2_new + 1]
    return wx_new , fx_new

#==================== Fourier transform===================================#

sigma_FT0 =      0.660
ycont     =          1
perc      =         10 # % of the line sampling that is smoothed
ampli     =         40  # Extension factor for the line continuum
cc        =     2.99e5   # Light speed km/s
bb        =          0

def win_cos(nx):
    wx = np.ones(nx)
    nxw = int(np.round(perc*nx/100))
    wi = 0.5*(1 - np.cos(np.pi*np.arange(nxw)/nxw))
    wx[0:nxw] = wi
    wx[(nx-nxw):nx] = wi[:nxw][::-1]
    return wx

def smooth_edge(wx, fx):
    fx = fx - bb
    nx = len(wx)
    correc = win_cos(nx)
    fx = fx*correc
    fx = fx + bb
    return fx

def extend_cont(wx, fx):
    nx = len(wx)
    plus = ampli * nx
    if plus != 0:
        paso = wx[1] - wx[0]
        wextens1 = (np.arange(plus) + 1) * paso + wx[nx - 1]
        wextens2 = -(-np.arange(plus) + plus) * paso + wx[0]
        wx = np.concatenate([wextens2, wx, wextens1])
        fextens1 = np.ones(plus) + bb
        fextens2 = np.ones(plus) + bb
        fx = np.concatenate([fextens1, fx, fextens2])
    return wx, fx

def vfft(wx, fx, xc):
    # Apply edge smoothing to the spectrum
    fx = smooth_edge(wx, fx)
    # Extend the continuum at both edges
    wf , ff = extend_cont(wx, fx)
    nx = len(wf)
    t = (wf[1] - wf[0])
    ds = 1 / (t * nx)
    max_sigma = 0.5 / t
    # Compute the Fourier transform of the spectrum
    res = abs(fft(ff - ycont))
    ew = res[0] * nx * t
    nr = len(res) - 1
    res = res / res[0]
    factor = xc / cc / sigma_FT0
    base   = (np.arange(nr-1)+1)*ds
    xrmax = 0.8 / t
    xrmin = base[np.argmin(res[0:nr-1] < 0.9)]
    xrmax = 1.0 / (xrmax * factor)
    xrmin = 1.0 / (xrmin * factor)
    max_sigma = 1.0 / (max_sigma * factor)
    xfou_obs = 1./(base*factor)
    with np.errstate(divide='ignore'):
        yfou_obs = np.log10(res[0:nr-1])
    return xfou_obs, yfou_obs, max_sigma

