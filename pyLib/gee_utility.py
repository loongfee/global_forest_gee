
# coding: utf-8

# In[ ]:


import ee
from IPython.display import Image
import datetime
import wget
import re
try:
    import numpy as np
except ImportError:
    print("You need numpy installed")
    sys.exit(-1)
try:
    from osgeo import gdal
    from osgeo import ogr
    from osgeo import osr
    from osgeo.gdalnumeric import *
    from osgeo.gdalconst import *
except ImportError:
    print("You need GDAL installed")
    sys.exit(-1)

def addL8Clear(image):
    return image.addBands(
            image.select('BQA').eq(2720).Or(
                image.select('BQA').eq(2724)).Or(
                image.select('BQA').eq(2728)).Or(
                image.select('BQA').eq(2732)).rename('clear'))
#    return image.addBands(
#            image.select('BQA').gte(2720).And(
#                image.select('BQA').lte(2732)).rename('clear'))

def addL5Clear(image):
    return image.addBands(
            image.select('BQA').eq(672).Or(
                image.select('BQA').eq(676)).Or(
                image.select('BQA').eq(680)).Or(
                image.select('BQA').eq(684)).rename('clear'))
#    return image.addBands(
#            image.select('BQA').gte(672).And(
#                image.select('BQA').lte(684)).rename('clear'))

def addS2Clear(image):
    return image.addBands(
            image.select('QA60').eq(0).And((\
                                            image.select('B2').eq(0)\
                                            .And(image.select('B3').eq(0))\
                                            .And(image.select('B4').eq(0))).Not()).rename('clear'))

def maskCloud(image):
    return image.mask(image.select('clear'))

def maskNoCloud(image):
    return image.mask(image.select('clear').Not())


def get_landsat_toa_collection(dataset_name_list):
    band_names = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
    for i, dataset_name in enumerate(dataset_name_list):
        str_list = dataset_name.split('/')
        if re.match(r'^LANDSAT/LC08.*', dataset_name):
            blue_band, green_band, red_band, nir_band, swir1_band, swir2_band =  'B2',  'B3',  'B4',  'B5',  'B6', 'B7'
            collection = ee.ImageCollection(dataset_name).map(addL8Clear)
        else:
            blue_band, green_band, red_band, nir_band, swir1_band, swir2_band =  'B1',  'B2',  'B3',  'B4',  'B5', 'B7'
            collection = ee.ImageCollection(dataset_name).map(addL5Clear)
        collection = collection.select([blue_band, green_band, red_band, nir_band, swir1_band, swir2_band, 'clear'],
                                       band_names + ['clear'])
        if i == 0:
            out_collection = collection
        else:
            out_collection = out_collection.merge(collection)
    return ee.ImageCollection(out_collection)


def get_sentinel_toa_collection():
    blue_band, green_band, red_band, nir_band, swir1_band, swir2_band =  'B2',  'B3',  'B4',  'B8',  'B11', 'B12'
    band_names = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
    collection = ee.ImageCollection('COPERNICUS/S2').map(addS2Clear)
    collection = collection.select([blue_band, green_band, red_band, nir_band, swir1_band, swir2_band, 'clear'],
                                    band_names + ['clear'])
    return collection


def NBR(image):
        return image.expression( 'float(NIR - SWIR2) / (NIR + SWIR2)', 
                            { 'NIR': image.select( 'NIR'), 'SWIR2': image.select('SWIR2')}).rename('NBR')  

    
def BAI(image):
    return image.expression('1.0/((0.06-NIR)**2 + (0.1-RED)**2)', 
                            { 'NIR': image.select( 'NIR'), 'RED': image.select('RED')}).rename('BAI')


def SWIR2(image):
    return image.expression('SWIR2', 
                            { 'SWIR2': image.select('SWIR2')}).rename('SWIR2')


def NBR2(image):
    return image.expression('(NIR-SWIR2)/(NIR+SWIR2)', 
                            { 'NIR': image.select( 'NIR'), 'SWIR2': image.select('SWIR2')}).rename('NBR2')

def NDMI(image):
    return image.expression('(NIR-SWIR1)/(NIR+SWIR1)', 
                            { 'NIR': image.select( 'NIR'), 'SWIR1': image.select('SWIR1')}).rename('NDMI')

def NBR3(image):
    return image.expression('(SWIR1-SWIR2)/(SWIR1+SWIR2)', 
                            { 'SWIR1': image.select('SWIR1'), 'SWIR2': image.select('SWIR2')}).rename('NBR3')

def MNBR(image):
    return image.expression('(RED+NIR-SWIR1-SWIR2)/(RED+NIR+SWIR1+SWIR2)', 
                            { 'NIR': image.select( 'NIR'), 'RED': image.select('RED'), 
                             'SWIR1': image.select('SWIR1'), 'SWIR2': image.select('SWIR2')}).rename('MNBR')

def NLBI(image):
    return image.expression('((RED+NIR-SWIR1-SWIR2)/(RED+NIR+SWIR1+SWIR2)-M)-L/(1.0+exp(-k*(1.0/((0.06-NIR)**2 + (0.1-RED)**2)-t)))', 
                            { 'NIR': image.select( 'NIR'), 'RED': image.select('RED'), 
                             'SWIR1': image.select('SWIR1'), 'SWIR2': image.select('SWIR2'),
                             'M':0, 'k': 0.1, 't': 120, 'L': 0.35}).rename('NLBI')

def NLBI(image):
    return image.expression('((RED+NIR-SWIR1-SWIR2)/(RED+NIR+SWIR1+SWIR2)-M)-L/(1.0+exp(-k*(1.0/((0.06-NIR)**2 + (0.1-RED)**2)-t)))', 
                            { 'NIR': image.select( 'NIR'), 'RED': image.select('RED'), 
                             'SWIR1': image.select('SWIR1'), 'SWIR2': image.select('SWIR2'),
                             'M':0, 'k': 0.1, 't': 120, 'L': 0.35}).rename('NLBI')

def NDVI(image):
    return image.expression('(NIR-RED)/(RED+NIR)', 
                            { 'NIR': image.select( 'NIR'), 'RED': image.select('RED')}).rename('NDVI')

def EVI(image):
    return image.expression('2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', 
                            { 'NIR': image.select( 'NIR'), 'RED': image.select('RED'), 'BLUE': image.select('BLUE')}).rename('EVI')

def GEMI(image):
    image = image.addBands(image.expression('(2*(data_nir**2-data_red**2)+1.5*data_nir+0.5*data_red)/(data_nir+data_red+0.5)', 
                                            { 'data_nir': image.select('NIR'), 'data_red': image.select('RED')}).rename('GEMI_a'))
    return image.expression('a*(1-0.25*a) - (data_red-0.125)/(1-data_red)', 
                            { 'a': image.select('GEMI_a'), 'GEMI_a': image.select('NIR'), 'data_red': image.select('RED')}).rename('GEMI')


def maskCloud(image):
    return image.mask(image.select('clear'))
    
def maskNoCloud(image):
    return image.mask(image.select('clear').Not())

def df_to_featureCollection(df, properties):
    feature_list = []
    for _, row in df.iterrows():
        feature = {'type':'Feature',
                   'properties':{},
                   'geometry':{'type':'Point',
                               'coordinates':[]}}
        feature['geometry']['coordinates'] = [0.0, 0.0]
        for prop in properties:
            feature['properties'][prop] = row[prop]
        feature_list += [ee.Feature(feature)]
    num = min(len(feature_list), 3000)
    return ee.FeatureCollection(feature_list[0:num])

def shapefile_to_json(shapfile):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(shapfile)
    layer = ds.GetLayer()
    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    for feature in layer:
        geomcol.AddGeometry(feature.GetGeometryRef())
    geojson = geomcol.ExportToJson()
    layer.ResetReading()
    return geojson

def read_point_collection(point_shapfile):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(point_shapfile)
    layer = ds.GetLayer()
    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    lnglat_list = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        for i in range(0, geom.GetPointCount()):
            # GetPoint returns a tuple not a Geometry
            pt = geom.GetPoint(i)
            lnglat_list += [[pt[0], pt[1]]]

    layer.ResetReading()
    return ee.Geometry.MultiPoint(lnglat_list)
    #return ee.Feature(ee.Geometry.MultiPoint(lnglat_list))
    
def read_polygon_collection(polygon_shapfile):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(polygon_shapfile)
    
    layer = ds.GetLayer()
    feature_list = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        geom = geom.GetGeometryRef(0)
        lnglat_list = []
        for i in range(0, geom.GetPointCount()):
            # GetPoint returns a tuple not a Geometry
            pt = geom.GetPoint(i)
            lnglat_list += [[pt[0], pt[1]]]
        #lnglat_list += [lnglat_list[0]]
        #print(len(lnglat_list))
        feature_list += [ee.Feature(ee.Geometry.Polygon(lnglat_list))]
    features = ee.FeatureCollection(feature_list)

    layer.ResetReading()
    return features