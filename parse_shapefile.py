#!/usr/bin/env python
"""
Title
Created on Wed May 01 12:54:46 2013

@author: u12161
"""

import os
import shapefile
import logging
from files import flModDate

try:
    import ogr
except ImportError:
    from osgeo import ogr

try:
    import osr
except ImportError:
    from osgeo import osr

LOG = logging.getLogger(__name__)

def getProjection(shape_file):
    """
    Get the projection of a shape file
    
    :param str shape_file: Name of a valid shape file

    :returns: :class:`osgeo.osr.SpatialReference` of the shape file

    """

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataset = driver.Open(shape_file)
    layer = dataset.GetLayer()
    spatial_ref = layer.GetSpatialRef()
    
    return spatial_ref

def writeProjectionFile(spatial_ref, filename):
    """
    Write a .prj file containing info on the projection of a
    shape file

    :param spatial_ref: :class:`osgeo.osr.SpatialReference` instance
    :param filename: File name to write to (without extension)

    """

    spatial_ref.MorphToESRI()
    file = open(filename + ".prj", 'w')
    file.write(spatial_ref.ExportToWkt())
    file.close()
    return


def parseShapefile(shape_file):
    """
    Read shapes and records from a shape file and return them.

    :param str shape_file: Shape file to parse.

    :returns: list of :class:`shapefile._Shape` instances (one for each
              feature), a list of lists containing the records for each 
              feature and list of field descriptions.

    """

    LOG.info('Parsing shapefile {0}'.format(os.path.abspath(shape_file)))

    moddate = flModDate(shape_file)
    LOG.debug('Last modified: {0}'.format(moddate))

    # read shapefile
    sf = shapefile.Reader(shape_file)
    shapes = sf.shapes()
    records = sf.records()
    fields = sf.fields

    # remove first header field - no corresponding value in "records"
    fields = fields[1:]

    LOG.debug("{0} features in {1}".format(len(shapes), shape_file))
    LOG.debug("{0} fields in {1}".format(len(fields) - 1, shape_file))
    return shapes, fields, records
    
def writeShapefile(output_file, fields, shapes, records):
    """
    Write the records and shapes to an output shape file.

    :param str output_file: Path for output file.
    :param list fields: List of field definitions.
    :param list shapes: List of :class:`shapefile._Shape` instances.
    :param list records: List of records for each feature.

    """
    LOG.info("Writing data to {0}.shp".format(os.path.abspath(output_file)))
    shp_writer = shapefile.Writer(shapefile.POLYGON)

    # Set the fields:
    for field in fields:
        shp_writer.field(*field)

    # Add shapes and records:
    for rec, shp in zip(records, shapes):
        if len(shp.parts) > 1:
            start = 0
            new_parts = []
            for part in shp.parts[1:]:
                new_parts.append(list(shp.points[start:part-1]))
                start = part
            new_parts.append(list(shp.points[part:-1]))
            shp_writer.poly(parts=new_parts)
        else:
            shp_writer.poly(parts=[shp.points])
        shp_writer.record(*rec)

    try:
        shp_writer.save(output_file)
    except shapefile.ShapefileException:
        LOG.exception("Failed to write to {0}".
                         format(os.path.abspath(output_file)))
        raise
    except AssertionError:
        LOG.exception("Problem in length and precision of fields")
        raise

    return
    
    
