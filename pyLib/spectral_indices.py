# coding: utf-8


import pandas as pd
import numpy as np

def landsat_nbr(df, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_swir2 = df['SWIR2'].astype(t)

    np.seterr(invalid='ignore')
    index_data = (data_nir - data_swir2) / (data_nir + data_swir2) / out_scale_factor
    return index_data


def landsat_bai(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_red = df['RED'].astype(t)

    nir_convergence = 0.06
    red_convergence = 0.1

    # The actual calculation
    data_red = data_red.astype(t) * in_scale_factor
    data_nir = data_nir.astype(t) * in_scale_factor

    data_nir = data_nir - nir_convergence
    data_red = data_red - red_convergence

    np.seterr(invalid='ignore')

    # data_bai = 1.0/((data4 - nir_convergence)**2 + (data3 - red_convergence)**2)
    data_bai = 1.0 / (data_nir ** 2 + data_red ** 2) / out_scale_factor

    return data_bai


def landsat_ndvi(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_red = df['RED'].astype(t)

    data_red = data_red.astype(t)*in_scale_factor
    data_nir = data_nir.astype(t)*in_scale_factor
    
    data_out = (data_nir-data_red)/(data_nir+data_red)
    data_out /= out_scale_factor

    return data_out


def landsat_evi(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_red = df['RED'].astype(t)
    data_blue = df['BLUE'].astype(t)
    
    G = 2.5
    C1 = 6
    C2 = 7.5
    L = 1

    data_red = data_red.astype(t)*in_scale_factor
    data_nir = data_nir.astype(t)*in_scale_factor
    data_blue = data_blue.astype(t)*in_scale_factor
    
    data_out = G * ((data_nir - data_red) / (data_nir + C1 * data_red - C2 * data_blue + L))
    data_out /= out_scale_factor

    return data_out


def landsat_gemi(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_red = df['RED'].astype(t)

    data_red = data_red.astype(t)*in_scale_factor
    data_nir = data_nir.astype(t)*in_scale_factor

    a = (2*(data_nir**2-data_red**2)+1.5*data_nir+0.5*data_red)/(data_nir+data_red+0.5)
    data_out = a*(1-0.25*a) - (data_red-0.125)/(1-data_red)
    data_out /= out_scale_factor

    return data_out


def landsat_mnbr(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_red = df['RED'].astype(t)
    data_swir1 = df['SWIR1'].astype(t)
    data_swir2 = df['SWIR2'].astype(t)

    data_red = data_red.astype(t)*in_scale_factor
    data_nir = data_nir.astype(t)*in_scale_factor
    data_swir1 = data_swir1.astype(t)*in_scale_factor
    data_swir2 = data_swir2.astype(t)*in_scale_factor

    data_out = (data_red+data_nir-data_swir1-data_swir2)/(data_red+data_nir+data_swir1+data_swir2)
    data_out /= out_scale_factor

    return data_out


def landsat_siri(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_swir1 = df['SWIR1'].astype(t)
    data_swir2 = df['SWIR2'].astype(t)
    
    data_swir1 = data_swir1.astype(t)*in_scale_factor
    data_swir2 = data_swir2.astype(t)*in_scale_factor

    data_out = (data_swir1-data_swir2)/(data_swir1+data_swir2)
    data_out /= out_scale_factor

    return data_out


def landsat_mirbi(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_swir1 = df['SWIR1'].astype(t)
    data_swir2 = df['SWIR2'].astype(t)

    data_swir1 = data_swir1.astype(t)*in_scale_factor
    data_swir2 = data_swir2.astype(t)*in_scale_factor

    data_out = (10*data_swir2-9.8*data_swir1+2)
    data_out /= out_scale_factor

    return data_out


def landsat_nbr2(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_swir1 = df['SWIR1'].astype(t)
    data_swir2 = df['SWIR2'].astype(t)
    
    data_swir1 = data_swir1.astype(t)*in_scale_factor
    data_swir2 = data_swir2.astype(t)*in_scale_factor

    data_out = (data_swir1-data_swir2)/(data_swir1+data_swir2)
    data_out /= out_scale_factor

    return data_out


def landsat_nbr3(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_swir2 = df['SWIR2'].astype(t)

    data_nir = data_nir.astype(t)*in_scale_factor
    data_swir2 = data_swir2.astype(t)*in_scale_factor

    data_out = (data_nir-data_swir2)/(data_nir+data_swir2)
    data_out /= out_scale_factor

    return data_out


def landsat_ndmi(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    data_nir = df['NIR'].astype(t)
    data_swir1 = df['SWIR1'].astype(t)

    data_nir = data_nir.astype(t)*in_scale_factor
    data_swir1 = data_swir1.astype(t)*in_scale_factor

    data_out = (data_nir-data_swir1)/(data_nir+data_swir1)
    data_out /= out_scale_factor

    return data_out



def landsat_nlbi(df, in_scale_factor=0.0001, out_scale_factor=1.0):
    t = np.float32
    # The actual calculation
    
    data_nir = df['NIR'].astype(t)
    data_red = df['RED'].astype(t)
    data_swir1 = df['SWIR1'].astype(t)
    data_swir2 = df['SWIR2'].astype(t)

    data_red = data_red.astype(t)*in_scale_factor
    data_nir = data_nir.astype(t)*in_scale_factor
    data_swir1 = data_swir1.astype(t)*in_scale_factor
    data_swir2 = data_swir2.astype(t)*in_scale_factor

    data_mnbr = (data_red+data_nir-data_swir1-data_swir2)/(data_red+data_nir+data_swir1+data_swir2)

    nir_convergence = 0.06
    red_convergence = 0.1

    # The actual calculation
    data_red = data_red.astype(t) * in_scale_factor
    data_nir = data_nir.astype(t) * in_scale_factor

    data_nir = data_nir - nir_convergence
    data_red = data_red - red_convergence

    np.seterr(invalid='ignore')
    data_bai = 1.0 / (data_nir ** 2 + data_red ** 2)
    
    # nlbi
    M = -0.42
    M = 0
    k = 0.2
    t = 110
    L = 0.36
    data_nlbi = (data_mnbr-M)-L/(1.0+np.exp(-k*(data_bai-t)))
    data_nlbi /= out_scale_factor

    return data_nlbi