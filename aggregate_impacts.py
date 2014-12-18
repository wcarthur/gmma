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
import argparse
import logging
from os.path import join as pjoin, abspath, isdir

try:
    import numpy as np
    from files import flConfigFile, flStartLog
    from datetime import datetime
    from parse_shapefile import parseShapefile
    import AggUtils

except ImportError as error:
    print "Cannot import all required modules"
    print "Import error: {0}".format( error )
    sys.exit(-1)

# Turn off runtime warnings:
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

LOG = logging.getLogger()

EPSILON = 1.0e-6
__version__ = "0.1"

def showSyntax(exit_code=0):
    """Display how to use this script and exit with the given exit
    code"""
    if sys.stdout.isatty():
        print sys.argv[0]
        print __doc__

    sys.exit(exit_code)

def aggregateLoss(records, fields, features, zoneid,
                  output_path):
    """
    Aggregate loss for all return periods and annualised loss across
    some field (e.g. suburb boundaries)
    """

    lguoutput, bgyoutput = AggUtils.aggregate_events(records, fields, zoneid,
                                              features)

    lgu_loss_file = pjoin(abspath(output_path), 'lgu_impact.csv')
    bgy_loss_file = pjoin(abspath(output_path), 'bgy_impact.csv')

    AggUtils.writeEventOutput(lgu_loss_file, lguoutput)
    AggUtils.writeEventOutput(bgy_loss_file, bgyoutput)

    return

def main():
    """
    Main section of the script - process command line args and call
    other functions to perform the aggregation
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--featureid', 
                        help=("Fieldname of a unique identifier for each "
                              "feature in the input dataset. Assumed to "
                              "be an integer."))
    parser.add_argument('-o', '--output',
                        help=("Path to folder for storing the output "
                              "- defaults to current working directory"))
    parser.add_argument('-s','--shapefile',
                        help=("Path (including extension) of the shape "
                              "file holding the zone features to process"))
    parser.add_argument('-z', '--zonefield',
                        help=("Name of the field in the input shape file "
                              "that will be used to aggregate results. "
                              "This field needs to be populated"))


    args = parser.parse_args()
    shape_file = args.shapefile
    output_path = args.output
    featureid = args.featureid
    zonefield = args.zonefield

    # Name of the output shape file, with no extension:
    if not isdir(abspath(output_path)):
        try:
            os.makedirs(abspath(output_path))
        except:
            print "Cannot create output path: {0}".format(output_path)
            raise

    # Load the exposure file:
    shapes, fields, records = parseShapefile(shape_file)
    features = AggUtils.loadZonesFromRecords(records, fields, featureid,
                                             zonefield)

    aggregateLoss(records, fields, features, featureid, output_path)

if __name__ == '__main__':

    if len(sys.argv) == 1:
        showSyntax()

    __STARTTIME__ = datetime.now()

    LOG = flStartLog(log_file=flConfigFile('.log'), log_level='INFO',
                        verbose=True, datestamp=True)

    main()
