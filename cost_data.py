#!/usr/bin/env python
"""
 Title: cost_data.py - access building construction cost data
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2013-05-01
 Description: Read in building construction cost data from csv file and produce
              a lookup table (dict of dicts) that holds the information for
              future processing.

 Version: 0.1
 Id: $Id$

 (C) Commonwealth of Australia (Geoscience Australia) 2012
 This product is released under the Creative Commons Attribution 3.0
 Australia Licence

 http://creativecommons.org/licenses/by/3.0/au/legalcode
 """

import os
import csv
from time import ctime
import numpy as np
import logging

LOG = logging.getLogger(__name__)

def parseCostFile(cost_file):
    """
    Read in building construction cost data from csv file and produce
    a lookup table (dict of dicts) that holds the information for
    future processing.
    """
    LOG.info('Reading cost data from %s' % os.path.abspath(cost_file))
    try:
        sinfo = os.stat(cost_file)
    except (IOError, WindowsError):
        LOG.exception( "Failed to open cost file: %s" % cost_file)
        raise IOError("Failed to open cost file: %s" % cost_file)
    else:
        moddate = ctime(sinfo.st_mtime)
        LOG.info('Last modified %s' % moddate)

    building_cost_data = dict()

    with open(cost_file,'rb') as cost_file_handle:
        reader = csv.reader(cost_file_handle, quotechar='"')

        for row in reader:
            if reader.line_num == 1:
                # First line contains headers.  The first two columns will
                # hold the L4 and L5 land use description. Remaining
                # columns are data:
                fields = row[2:]
            else:
                # Reading the data:
                lu4 = row[0]
                lu5 = row[1]
                data = row[2:]
                entry = dict([[fields[i], int(data[i])]
                              for i in range(len(data))
                              if data[i] != ''])

                if building_cost_data.has_key(lu4):
                    building_cost_data[lu4][lu5] = entry
                else:
                    building_cost_data[lu4] = dict()
                    building_cost_data[lu4][lu5] = entry

    return building_cost_data

def calculateValue(floor_area, landuse4, landuse5, building_type,
                    building_costs):
    """
    Calculate the value of a building type, based on it's floor area,
    building construction type and land-use classification (level 4
    and 5). Data stored in the building_costs dict is the cost per
    square metre for construction, so we need to multiply by the floor
    area to obtain the total value for the building type.

    The building type is a string comprised of the construction type
    concatenated with the building height (e.g. L_1, L_2, M, H, V, E
    or S). Splitting on the underscore gives the building type.

    The landuse values correspond to the level 4 and 5 land use
    classifications used in the GMMA exposure database. See the
    documentation for that data set for information on the
    classifications used.

    For those records where there is no matching construction cost
    information (for the given land use classes), the value of that
    building type is set to zero.

    :param floor_area: Floor area values for the given building type, across
                       all polygons.
    :type  floor_area: :class:`numpy.ndarray`
    :param list landuse4: Level 4 land use class for the polygons.
    :param list landuse5: Level 5 land use class for the polygons.
    :param str building_type: Construction type, including the building height
                              class.
    :param dict building_costs: A dict that holds the building costs (per
                                sq. m value) for the full range of building
                                types and (valid) combinations of level 4 & 5
                                land use class.
                                
    """

    LOG.debug("Calculating value of building type {0}".format(building_type))
    # Key part of the building type is the construction method and
    # building height:
    bldg_type = building_type.split('_', 1)[0]

    value = np.empty(len(floor_area))

    for i, [lu4, lu5, fla] in enumerate(zip(landuse4, landuse5, floor_area)):
        # Hard-coded fix for coverage issues (detailed by G. Davies
        # 2013-03-12)
        if (lu4 == 'Informal Settlements') and (lu5 == 'Residential'):
            lu4 = 'Formal Settlements'
        elif (lu4 == 'Heavy Industry') and (lu5 == 'Services'):
            lu4 = 'Government'
        elif (lu4 == 'Formal Settlements') and (lu5 == 'Yet to be defined'):
            lu5 = 'Mixed Residential and Small Commercial'
        elif (lu4 == 'Heavy Industry') and (lu5 == 'Office'):
            lu5 = 'Manufacturing'
        elif (lu4 == 'Cultural') and (lu5 == 'Exhibitions'):
            lu4 = 'Leisure'
        elif (lu4 == 'Vacant Areas'):
            lu5 = 'Yet to be defined'
        elif (lu4 == 'Natural Areas'):
            lu5 = 'National Parks'


        try:
            if building_costs[lu4][lu5].has_key(bldg_type):
                # Pull the cost/sq. metre data from the building_costs dict:
                value[i] = building_costs[lu4][lu5][bldg_type] * fla
            else:
                value[i] = 0.0
        except KeyError:
            LOG.warn(("No building cost information for {0} buildings "
                      "in {1}-{2} class land parcels "
                      "(parcel ID# {3})").format(bldg_type, lu4, lu5, i))
            value[i] = 0.0
    return value

