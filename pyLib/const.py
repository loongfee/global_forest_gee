# coding: utf-8


import pandas as pd
import numpy as np

grid_name_dict = {}
grid_name_dict['AUST'] = ['10S_110E', '10S_120E', '10S_130E', '10S_140E',\
                                                  '20S_110E', '20S_120E', '20S_130E', '20S_140E', '20S_150E',\
                                                  '30S_110E', '30S_120E', '30S_130E', '30S_140E', '30S_150E',\
                                                                                                                  '40S_140E',\
                                                  '30S_170E', '40S_160E', '40S_170E', '50S_160E',]

grid_name_dict['EQAS'] = ['10S_110E', '10S_120E', '10S_130E', '10S_140E',\
                                                  '20S_110E', '20S_120E', '20S_130E', '20S_140E', '20S_150E',\
                                                  '30S_110E', '30S_120E', '30S_130E', '30S_140E', '30S_150E',\
                                                                                                                  '40S_140E',\
                                                  '30S_170E', '40S_160E', '40S_170E', '50S_160E',]


grid_name_dict['NHAF'] = ['20N_020W', '20N_010W', '20N_000E', '20N_010E', '20N_020E', '20N_030E', '20N_040E', '20N_050E',\
                                                   '10N_020W', '10N_010W', '10N_000E', '10N_010E', '10N_020E', '10N_030E', '10N_040E', '10N_050E',]

grid_name_dict['SHAF'] = ['00N_000E', '00N_010E', '00N_020E', '00N_030E', '00N_040E', \
                                                  '10S_010E', '10S_020E', '10S_030E', '10S_040E', '10S_050E',  \
                                                  '20S_010E', '20S_020E', '20S_030E', '20S_040E', \
                                                  '30S_010E', '30S_020E', '30S_030E', ]

grid_name_dict['CEAS'] = ['50N_010W', '50N_000E', '50N_010E', '50N_020E', '50N_030E', '50N_040E', '50N_050E', '50N_060E',\
                                                  '50N_070E', '50N_080E', '50N_090E', '50N_100E', '50N_110E', '50N_120E', '50N_130E', '50N_140E', \
                                                  '40N_050E', '40N_060E', '40N_070E', '40N_080E', '40N_090E', '40N_100E', '40N_110E', '40N_120E',\
                                                  '40N_130E', '40N_140E', \
                                                  '30N_080E', '30N_090E', '30N_100E', '30N_110E', '30N_120E',]

#grid_name_dict['CHINA'] = ['50N_060E', '50N_070E', '50N_080E', '50N_090E', '50N_100E', '50N_110E',\
#                                                  '40N_050E', '40N_060E', '40N_070E', '40N_080E', '40N_090E', '40N_100E', '40N_110E', '40N_120E', '40N_130E', '40N_140E',\
#                                                  '30N_060E', '30N_070E', '30N_080E', '30N_090E', '30N_100E', '30N_110E', '30N_120E',\
#                                                  '20N_070E', '20N_080E', '20N_090E', '20N_100E', '20N_110E', '20N_120E',\
#                                                  '10N_070E', '10N_080E', '10N_090E', '10N_100E', '10N_110E', '10N_120E',]

grid_name_dict['TENA'] = ['50N_130W', '50N_120W', '50N_110W', '50N_100W', '50N_090W', '50N_080W', '50N_070W', '50N_060W',\
                                                  '40N_130W', '40N_120W', '40N_110W', '40N_100W', '40N_090W', '40N_080W',]

grid_name_dict['CEAM'] = ['30N_120W', '30N_110W','30N_100W','30N_090W','30N_080W',\
                                                  '20N_110W','20N_100W','20N_090W','20N_080W','20N_070W',\
                                                  '20N_090W','20N_080W',]

grid_name_dict['SHSA'] = ['00N_090W', '00N_080W', '00N_070W', '00N_060W', '00N_050W', '00N_040W',\
                                                                         '10S_080W', '10S_070W', '10S_060W', '10S_050W', '10S_040W',\
                                                                         '20S_080W', '20S_070W', '20S_060W', '20S_050W', \
                                                                         '30S_080W', '30S_070W', '30S_060W', \
                                                                         '40S_080W', '40S_070W', \
                                                                         '50S_080W', '50S_070W', '50S_060W', ]

#grid_name_dict['SEAS'] = ['00N_090E', '00N_100E', '00N_110E', '00N_120E', '00N_130E', '00N_140E', '00N_150E', '00N_160E',]
grid_name_dict['SEAS'] = ['30N_090E', '20N_070E', '20N_080E', '20N_090E', '20N_100E', '20N_110E', \
                           '10N_070E', '10N_080E', '10N_090E', '10N_100E', '10N_110E', '10N_120E',\
                         '00N_090E', '00N_100E', '00N_110E',  '00N_120E', '00N_130E', '00N_140E', '00N_150E', '00N_160E',]

grid_name_dict['NCHINA'] = ['60N_100E', '60N_110E', '60N_120E', '60N_130E', \
                                                        '50N_100E', '50N_110E', '50N_120E', '50N_130E']

grid_name_dict['WESTA'] = [ '60N_010E', '60N_020E', '60N_030E', '60N_040E', '60N_050E', '60N_060E', '60N_070E', '60N_080E', '60N_090E', \
                                                        '50N_010E', '50N_020E', '50N_030E', '50N_040E', '50N_050E', '50N_060E', '50N_070E', '50N_080E']

grid_name_dict['BOAS'] = [ '60N_100E', '60N_110E', '60N_120E', '60N_130E', '60N_140E', '60N_150E', '60N_160E', '60N_170E', \
                                                         '70N_100E', '70N_110E', '70N_120E', '70N_130E', '70N_140E', '70N_150E', '70N_160E', '70N_170E']



grid_name_dict['CHINA'] = ['50N_070E', '50N_080E', '50N_090E', '50N_100E', '50N_110E', '50N_120E', '50N_130E',\
                                                  '40N_070E', '40N_080E', '40N_090E', '40N_100E', '40N_110E', '40N_120E',\
                                                  '30N_080E', '30N_090E', '30N_100E', '30N_110E', '30N_120E',\
                                                  '20N_100E', '20N_110E',\
                                                  '60N_110E', '60N_120E',]

grid_name_dict['TEST'] = ['10N_020E']

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