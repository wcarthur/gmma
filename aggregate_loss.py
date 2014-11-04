#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Title: aggregate_loss.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: Wed May 01 16:17:09 2013
 Description: Aggregate the fields on the basis of a given aggregation field.

 Users need to set the return periods in the array at the bottom of this file.

 Input:
 This script reads a series of command line arguments to establish the
 parameters:
 -h, --help: displays this help message;
 -z, --zonefield: Name of the field in the input shape file
         that will be used to aggregate results. This field needs to be
         populated;
 -f, --featureid: Fieldname of a unique identifier for each feature in the
         input dataset. Assumed to be an integer.
 -o, --output: Path to folder for storing the output - defaults to current
         working directory;
 -s, --shapefile: Path (including extension) of the shape file holding the
                 zone features to process;

 Output:
 Data are aggregated by barangay (suburb) and municipality and values are
 written to CSV files for the damage, costs and damaged floor-area-equivalent.
 The script will also plot probability-loss curves for each municipality.

 To change the appearance of the plots, make changes to the plot_loss_output
 function in AggUtils.py


 Usage:
 1) Ensure the shapefile to be used has fields representing the damage, costs
    and damaged floor-area-equivalent in each feature, for a number of return
    periods.
 2) Update the 'returnPeriods' variable at the bottom of this file to contain
    the same values as you have damage fields in the shape file. e.g. if you
    have fields named 'V5', 'V10', 'V20' and 'V25', then the 'returnPeriods'
    variable should be [5,10,20,25].
 3) Run the script:
    e.g. at a command prompt, enter the following:
    C:\\Python26\\python aggregate_loss.py -f "OID1" -o ..\\Output
        -s ..\\Output\\annualised_loss.shp -z "AGGPOLY_ID"

 Version: 0.1 2013-05-24


 (C) Commonwealth of Australia (Geoscience Australia) 2012
 This product is released under the Creative Commons Attribution 3.0
 Australia Licence

 http://creativecommons.org/licenses/by/3.0/au/legalcode

 Id: $Id$

"""

import os
import sys
import getopt
import logging
from os.path import join as pjoin

try:
    import numpy as np
    from files import flConfigFile, flStartLog
    from datetime import datetime
    from parse_shapefile import parse_shapefile
    import AggUtils

except ImportError as error:
    print "Cannot import all required modules"
    print "Import error: {0}".format( error )
    sys.exit(-1)

# Turn off runtime warnings:
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

logger = logging.getLogger()

__eps__ = 1.0e-6
__version__ = "0.1"

def ShowSyntax(exit_code=0):
    """Display how to use this script and exit with the given exit
    code"""
    if sys.stdout.isatty():
        print sys.argv[0]
        print __doc__

    sys.exit(exit_code)

def parse_args(argv):
    """
    Parse the command line arguments
    """
    # Default values:
    shape_file = None
    output_path = os.getcwd()
    featureid = None
    zonefield = None

    try:
        opts, args = getopt.getopt(argv, 'hf:o:s:z:',
                                   ['help', 'featureid=', 'outputpath=',
                                    'shapefile=', 'zonefield='])
    except getopt.GetoptError:
        ShowSyntax(-1)
    else:
        for opt, arg in opts:
            if opt in ( "-h", "--help" ):
                ShowSyntax(0)
            elif opt in ( "-f", "--featureid"):
                featureid = arg
            elif opt in ( "-o", "--outputpath" ):
                output_path = arg
            elif opt in ( "-s", "--shapefile" ):
                shape_file = arg
            elif opt in ( "-z", "--zonefield" ):
                zonefield = arg

        return shape_file, output_path, featureid, zonefield

def aggregate_loss(records, fields, features, return_periods, zoneid,
                   output_path):
    """
    Aggregate loss for all return periods and annualised loss across
    some field (e.g. suburb boundaries)
    """

    lguoutput, bgyoutput = AggUtils.aggregate(records, fields, zoneid,
                                              features, return_periods)

    lgu_loss_file = pjoin(os.path.abspath(output_path), 'lgu_loss.csv')
    lgu_cost_file = pjoin(os.path.abspath(output_path), 'lgu_cost.csv')
    lgu_dmg_file = pjoin(os.path.abspath(output_path), 'lgu_dmg.csv')

    bgy_loss_file = pjoin(os.path.abspath(output_path), 'bgy_loss.csv')
    bgy_cost_file = pjoin(os.path.abspath(output_path), 'bgy_cost.csv')
    bgy_dmg_file = pjoin(os.path.abspath(output_path), 'bgy_dmg.csv')


    AggUtils.write_loss_output(lgu_loss_file, lguoutput, return_periods)
    AggUtils.write_loss_output(bgy_loss_file, bgyoutput, return_periods, False)

    AggUtils.write_cost_output(lgu_cost_file, lguoutput, return_periods)
    AggUtils.write_cost_output(bgy_cost_file, bgyoutput, return_periods, False)

    AggUtils.write_dmg_output(lgu_dmg_file, lguoutput, return_periods)
    AggUtils.write_dmg_output(bgy_dmg_file, bgyoutput, return_periods, False)

    return

def main(argv=None, return_periods=None):
    """
    Main section of the script - process command line args and call
    other functions to perform the aggregation
    """

    shape_file, output_path, featureid, zonefield = parse_args(argv)

    # Name of the output shape file, with no extension:
    if not os.path.isdir(os.path.abspath(output_path)):
        try:
            os.makedirs(os.path.abspath(output_path))
        except:
            print "Cannot create output path: {0}".format(output_path)
            raise

    # Load the exposure file:
    shapes, fields, records = parse_shapefile(shape_file)
    features = AggUtils.load_zones_from_records(records, fields, featureid,
                                                zonefield)
    #features = AggUtils.load_zones(zonefile)
    aggregate_loss(records, fields, features, return_periods, featureid,
                   output_path)

if __name__ == '__main__':

    if len(sys.argv) == 1:
        ShowSyntax()

    __STARTTIME__ = datetime.now()

    logger = flStartLog(log_file=flConfigFile('.log'), log_level='INFO',
                        verbose=True, datestamp=True)

    # Set the return periods here:
    return_periods = np.array([2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                              75, 100, 150, 200, 250, 300, 350, 400, 450,
                              500, 1000, 2000, 2500, 5000, 10000])

    main(sys.argv[1:], return_periods)
