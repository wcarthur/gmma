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

LOG = logging.getLogger(__name__)

def parseShapefile(shape_file):
    """
    Read shapes and records from a shape file and return them
    
    """
    LOG.info('Parsing shapefile {0}'.format(os.path.abspath(shape_file)))

    moddate = flModDate(shape_file)
    LOG.info('Last modified: {0}'.format(moddate))

    # read shapefile
    sf = shapefile.Reader(shape_file)
    shapes = sf.shapes()
    records = sf.records()
    fields = sf.fields

    # remove first header field - no corresponding value in "records"
    fields = fields[1:]
    
    return shapes, fields, records
    
def writeShapefile(output_file, fields, shapes, records):
    """
    Write the records and shapes to an output shape file
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
    
    
