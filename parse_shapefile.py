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
import pdb


logger = logging.getLogger(__name__)

def parse_shapefile(shape_file):
    """
    Read shapes and records from a shape file and return them
    
    """
    logger.info('Parsing shapefile {0}'.format(os.path.abspath(shape_file)))

    moddate = flModDate(shape_file)
    logger.info( 'Last modified: {0}'.format(moddate) )

    # read shapefile
    shp_reader = shapefile.Reader(shape_file)
    shapes = shp_reader.shapes()
    records = shp_reader.records()
    fields = shp_reader.fields

    # remove first header field - no corresponding value in "records"
    #fields_append = fields # keep this for writing new shp
    fields = fields[1:]
    #fieldnames = [fields[i][0] for i in range(len(fields))]
    
    return shapes, fields, records
    
def write_shapefile(output_file, fields, shapes, records):
    """
    Write the records and shapes to an output shape file
    """
    logger.info("Writing data to {0}.shp".format(os.path.abspath(output_file)))
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
        logger.exception("Failed to write to {0}".
                         format(os.path.abspath(output_file)))
        raise
    except AssertionError:
        logger.exception("Problem in length and precision of fields")
        raise

    return
    
    