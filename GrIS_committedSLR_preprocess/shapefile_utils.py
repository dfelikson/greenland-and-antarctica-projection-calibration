
import os
import ogr
import numpy as np

import matplotlib.path as mplPath

def createBuffer(inputfn, outputBufferfn, bufferDist):
   inputds = ogr.Open(inputfn)
   inputlyr = inputds.GetLayer()

   shpdriver = ogr.GetDriverByName('ESRI Shapefile')
   if os.path.exists(outputBufferfn):
      shpdriver.DeleteDataSource(outputBufferfn)
   outputBufferds = shpdriver.CreateDataSource(outputBufferfn)
   
   bufferlyr = outputBufferds.CreateLayer(outputBufferfn, geom_type=ogr.wkbPolygon)
   featureDefn = bufferlyr.GetLayerDefn()
   
   inputlyrDefn = inputlyr.GetLayerDefn()
   for i in range(inputlyrDefn.GetFieldCount()):
      bufferlyr.CreateField ( inputlyrDefn.GetFieldDefn(i) )
    
   for feature in inputlyr:
      ingeom = feature.GetGeometryRef()
      geomBuffer = ingeom.Buffer(bufferDist)
      
      outFeature = ogr.Feature(featureDefn)
      outFeature.SetGeometry(geomBuffer)
      
      for i in range(inputlyrDefn.GetFieldCount()): 
         fieldName = inputlyrDefn.GetFieldDefn(i).GetName()
         outFeature.SetField(fieldName, feature.GetField(fieldName))
      bufferlyr.CreateFeature(outFeature)
      outFeature = None

def read_shapefile_points(shapefile):
   if not os.path.isfile(shapefile):
      print('ERROR: shapefile not found: ' + shapefile)
      return None, None
   if shapefile.endswith('.shp'):
      ds = ogr.Open(shapefile)
      lyr = ds.GetLayer()
      feat = lyr.GetFeature(0)
      geom = feat.GetGeometryRef()
      
      x = np.array([])
      y = np.array([])
      
      # Old way:
      #points = geom.GetPoints()
      #for i, pointxy in enumerate(points):
      #   xtmp, ytmp = pointxy
      #   x  = np.append(x,  xtmp)
      #   y  = np.append(y,  ytmp)

      ring = geom.GetGeometryRef(0)
      points = mplPath.Path(ring.GetPoints())

      if points is None:
         print('ERROR: no points found in shapefile')
         return None, None

      x, y = zip(*points.to_polygons()[0])

   elif shapefile.endswith('.csv'):
      data = np.genfromtxt(shapefile, delimiter=',')
      data = data[~np.isnan(data).any(axis=1)]
      x = data[:,xcol]
      y = data[:,ycol]
     
   return x, y

def checkPointInShapefile(x, y, polygon_shapefile, field_name=None):
   field_values_containing_point = list()

   # Read polygon vertices
   ds = ogr.Open(polygon_shapefile)
   lyr = ds.GetLayer()
   for feat in lyr:
      geom = feat.GetGeometryRef()
      ring = geom.GetGeometryRef(0)
      #bbPath = mplPath.Path(np.array([[ring.GetPoint(0)[0], ring.GetPoint(0)[1]],
      #                                [ring.GetPoint(1)[0], ring.GetPoint(1)[1]],
      #                                [ring.GetPoint(2)[0], ring.GetPoint(2)[1]],
      #                                [ring.GetPoint(3)[0], ring.GetPoint(3)[1]]]))
      bbPath = mplPath.Path(ring.GetPoints())

      if np.any(bbPath.contains_point( (x, y) )):
         field_values_containing_point.append(feat.GetField(field_name))

   return field_values_containing_point

