#!/usr/bin/env python
"""
 Title: calculate_event_loss.py - parse exposure shape file for GMMA
        project
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2012-11-02
 Description: Process a shape file that holds proportional building type
 information and hazard information (return period wind speeds) for land
 parcels to calculate the loss for each return period, the annualised
 (integrated) loss and the total loss across the extent of the shape file.

 Dependencies: shapefile.py (http://code.google.com/p/pyshp/)
               scipy
               numpy
               matplotlib

 Users need to set the return periods in the array at the bottom of this file.

 Input:
 This script reads a series of command line arguments to establish the
 parameters:
 -h, --help: displays this help message;
 -a, --aggregatefield: Name of the field in the input shape file
                       that will be used to aggregate results. This
                       field needs to be populated;
 -c, --costs: Path to a csv format file containing the cost data for building
              types and land-use groupings;
 -o, --output: Path to folder for storing the output;
 -s, --shapefile: Path (including extension) of the shape file holding the
                 zone features to process (e.g. land-use parcels, meshblocks,
                 etc);
 -v, --vulnerability: Path to a csv format file containing the mean, sigma and
                      scale values for building vulnerability curves.

 Output:
 The input zonal feature dataset (shapefile) is copied to the output folder and
 updated with two new fields for each return period, plus new fields for
 annualised loss and costs.

 The new fields are named 'cost<year>' and 'dmg<year>' where <year>
 corresponds to each of the specified return periods (see below). The
 annualised costs are stored in 'ann_costs' and the annualised losses
 (representing loss as a fraction of replacement cost) in 'ann_loss'.

 Usage:
 1) Ensure the shapefile to be used has fields representing the
    wind hazard in each feature, for a number of return periods.
 2) Update the 'returnPeriods' variable at the bottom of this
    file to contain the same values as you have hazard fields in
    he shape file. e.g. if you have fields named 'V5', 'V10',
    'V20' and 'V25', then the  'returnPeriods' variable should
    be [5,10,20,25].
 3) Run the script:
    e.g. at a command prompt, enter the following:
    C:\\Python26\\python parse_gmma_exposure_shpfile_costs.py
       -c ..\\Exposure\\Building_types_by_land_use.csv
       -o ..\\Output
       -s ..\\Exposure\\GMMA_Exposure_Compilation_Final_Feb2013.shp
       -v ..\\Vulnerability\\wind_vuln_final.csv


 Version: 0.1 2012-11-02
          0.2 2012-11-27 - correct summation of damage for each parcel
          0.3 2012-11-28 - changed to use command line arguments
          0.4 2013-01-22 - calculate proportion of damage for different
                           building height classes (i.e. low, medium or
                           high-rise buildings)
          0.5 2013-03-07 - added parser for construction cost data
          0.6 2013-05-03 - calculate fractional loss based on costs, not
                           weighted by floor area.


 (C) Commonwealth of Australia (Geoscience Australia) 2012
 This product is released under the Creative Commons Attribution 3.0
 Australia Licence

 http://creativecommons.org/licenses/by/3.0/au/legalcode

"""

__version__ = "0.7"

import os
import sys
import csv
import argparse
import logging
import time
from datetime import datetime
from functools import wraps
import traceback
import pdb

from os.path import join as pjoin, isdir, abspath
try:
    import numpy as np
    from files import flConfigFile, flModDate, flStartLog
    from scipy.integrate import simps
    import matplotlib
    matplotlib.use('Agg', warn=False)

    from matplotlib import pyplot
    from get_records import getField
    from cost_data import parseCostFile, calculateValue
    from parse_shapefile import writeShapefile, parseShapefile
    from damage import damage, adjustDamageCurves, adjustFragilityCurves
    from probability import probability
    from AvDict import AvDict
except ImportError as error:
    print "Cannot import all required modules"
    print "Import error: {0}".format( error )
    sys.exit(-1)

# Turn off runtime warnings:
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

LOG = logging.getLogger()

def timer(func):
    """A simple timing decorator"""
    
    @wraps(func)
    def wrap(*args, **kwargs):
        t1 = time.time()
        res = func(*args, **kwargs)

        tottime = time.time() - t1
        msg = "%02d:%02d:%02d " % \
          reduce(lambda ll, b : divmod(ll[0], b) + ll[1:],
                        [(tottime,), 60, 60])

        LOG.info("Time for {0}: {1}".format(func.func_name, msg) )
        return res

    return wrap


def showSyntax(exit_code=0):
    """Display how to use this script and exit with the given exit
    code"""
    if sys.stdout.isatty():
        print sys.argv[0]
        print __doc__

    sys.exit(exit_code)

def totalDamage(building_types, fields, records, building_costs):
    """
    Calculate the total damage for each feature, for each return
    period.  Total damage is the weighted sum of damage to each
    building type, weighted by the fraction of floor area for that
    building type.
    """

    LOG.info("Calculating total damage to building stock")
    LOG.debug("Calculating total value of building stock in each feature")

    nrecords = len(records)
    n_bldg_type = len(building_types.keys())

    lu4 = getField('L4_USE', fields, records, str)
    lu5 = getField('L5_USE', fields, records, str)
    vintage = getField('ERA_CONST', fields, records, str)
    area_sqm = getField('AREA_SQM', fields, records)

    values = np.empty((n_bldg_type, nrecords))

    for i, bld_type in enumerate(building_types.keys()):
        flarea = getField(bld_type, fields, records)
        values[i, :] = calculateValue(flarea, lu4, lu5, bld_type,
                                        building_costs)

    totvalue = np.sum(values, axis=0)
    fields.append(['bldg_value', "N", 19, 0])

    for i, record in enumerate(records):
        record.append(totvalue[i])

    LOG.debug("Appending required fields for loss metrics")
    
    fields.append(['loss', "N", 9, 6])
    fields.append(['cost', "N", 10, 0])
    fields.append(['dmgf', "N", 10, 6])

    # Create a mask to allow us to modify the vulnerability for
    # buildings of different ages. The key time is 1992. This applies
    # to the MWS, CWS, CHB and C1-L building types, where there were
    # changes to the construction materials used.
    vmask = np.array(vintage)=='Post-1992'

    LOG.info("Processing wind speed")
    wind_speed = getField('vmax', fields, records)

    costs = np.empty((len(building_types.keys()), nrecords))
    dmg_flarea = np.empty((len(building_types.keys()), nrecords))

    for i, bld_type in enumerate(building_types.keys()):
        # Calculate corresponding value of built assets:
        flarea = getField(bld_type, fields, records)
        value = calculateValue(flarea, lu4, lu5, bld_type, building_costs)

        # Calculate the damage (as a fraction of replacement cost):
        mu = building_types[bld_type]['mu']*np.ones(len(wind_speed))
        sigma = building_types[bld_type]['sigma']*np.ones(len(wind_speed))
        scale = building_types[bld_type]['scale']*np.ones(len(wind_speed))

        mu, sigma, scale = adjustDamageCurves(bld_type, vmask, mu, sigma, scale)

        d = damage(wind_speed, mu, sigma, scale)

        costs[i, :] = d*value
        dmg_flarea[i, :] = 100.*d*flarea/area_sqm

    totcost = np.sum(costs, axis=0)
    totdmg = totcost/totvalue
    np.putmask(totdmg, totvalue==0, 0)

    # Store the damaged floor area equivalent in units of hectare/sq km.
    totdmgfl = np.sum(dmg_flarea, axis=0)

    del costs
    del dmg_flarea

    for i, record in enumerate(records):
        record.extend([totdmg[i], totcost[i], totdmgfl[i]])

    return fields, records

def damageState(building_types, fields, records, dmgstate):
    """
    Calculate probability of damage states for all building types.

    :param dict building_types: Dict of building type data, including the
                                parameters for the fragility curves. 
    :param list fields: List of field names from the input file.
    :param list records: List of lists of records from the input file.
    :param str state: Damage state to calculate. Must be one of
                      'slight', 'moderate', 'extensive' or 'complete'.

    """
    if dmgstate not in ['slight', 'moderate', 'extensive', 'complete']:
        raise KeyError("Invalid damage state requested: {0}".format(dmgstate))

    nrecords = len(records)
    n_bldg_type = len(building_types.keys())
    vintage = getField('ERA_CONST', fields, records, str)
    area_sqm = getField('AREA_SQM', fields, records)
    
    if dmgstate == 'slight':
        state = 'sl'
    elif dmgstate == 'moderate':
        state = 'mod'
    elif dmgstate == 'extensive':
        state = 'ext'
    elif dmgstate == 'complete':
        state = 'c'
        
    values = np.empty((n_bldg_type, nrecords))

    vmask = np.array(vintage)=='Post-1992'
    wind_speed = getField('vmax', fields, records)
    prob_state = np.empty((len(building_types.keys()), nrecords))
    LOG.debug("Appending required fields for damage state metrics")

    for i, bld_type in enumerate(building_types.keys()):
        LOG.debug("New field name is: {0}".format('_'.join([bld_type, state])))
        fieldname = '_'.join([bld_type, state])
        fields.append([fieldname, "N", 10, 6])
        flarea = getField(bld_type, fields, records)

        # Calculate the probability of being in a damage state:
        mu = building_types[bld_type][state+'_mu']*np.ones(len(wind_speed))
        sigma = building_types[bld_type][state+'_sd']*np.ones(len(wind_speed))
        scale = building_types[bld_type][state+'_scale']*np.ones(len(wind_speed))

        mu, sigma, scale = adjustFragilityCurves(bld_type, vmask, mu,
                                              sigma, scale, dmgstate)

        prob_state[i,:] = damage(wind_speed, mu, sigma, scale)
        
    for i, record in enumerate(records):
        record.extend(*[prob_state[:,i]])

    return fields, records

def calculatePopulation(fields, records, vulnerability_file):
    """
    Calculate the population affected statistic for the event.
    "Population affected" is defined to be the population in building types
    that have 70% or greater probability of being in a moderate damage state,
    or 30% or greater probability of being in an extensive damage state.

    The population in a land parcel is proportionally distributed amongst all
    building types in the parcel.

    :param list fields: List of fields in the dataset.
    :param list records: List of lists of records from the input file.
    :param str vulnerability_file: Filename of csv-format file that holds the
                                   parameters for the vulnerability curves
                                   for all buildings.
    """

    LOG.info("Calculating affected population")
    building_types = parseVulnerability(vulnerability_file)
    nrecords = len(records)
    n_bldg_type = len(building_types.keys())
    flarea = getField('FLAREA_SUM', fields, records)
    pop = getField('POP_EST', fields, records)
    pop_affected = np.zeros((n_bldg_type, nrecords))
    
    for i, bld_type in enumerate(building_types.keys()):
        mod_dmg_name = bld_type + '_mod'
        ext_dmg_name = bld_type + '_ext'
        mod_dmg_prob = getField(mod_dmg_name, fields, records)
        ext_dmg_prob = getField(ext_dmg_name, fields, records)

        # Damage thresholds: >70% probability of moderate damage; 
        #                    >30% probability of extensive damage.
        mod_idx = np.where(mod_dmg_prob >= 0.7)
        ext_idx = np.where(ext_dmg_prob >= 0.3)

        bld_flarea = getField(bld_type, fields, records)
        xx = bld_flarea * pop / flarea

        pop_affected[i, mod_idx] = xx[mod_idx]
        pop_affected[i, ext_idx] = xx[ext_idx]

    population_affected = np.sum(pop_affected, axis=0)

    np.putmask(population_affected, population_affected > pop, pop)

    fields.append(['POP_AFFECT', "N", 10, 2])
    for i, record in enumerate(records):
        record.extend([population_affected[i]])

    return fields, records
    
def processDamage(fields, records, vulnerability_file, cost_file):
    """
    Process the data to calculate loss information:
    Input arguments:
        fields - list of fields contained in the input exposure shape file
        records - list of a list of records for each feature in the exposure
            shape file
        vulnerability_file - filename of a csv-format file that holds the
            parameters for the vulnerability curves for all building types
        cost_file - filename for a csv-format file that holds the construction
            cost information for the combination of land use categories and
            building types
        return_periods - array of return period values for which to calculate
            losses.

    Output:
        Returns updated field and record objects wth additional fields and
        values corresponding to the return period losses and damage that
        are calculated.
    """

    building_types = parseVulnerability(vulnerability_file)
    building_costs = parseCostFile(cost_file)

    output_fields, \
        output_records = totalDamage(building_types, fields,
                                     records, building_costs)
    return output_fields, output_records

def processFragility(fields, records, vulnerability_file, state):
    building_types = parseVulnerability(vulnerability_file)
    output_fields, \
      output_records = damageState(building_types, fields, records, state)

    return output_fields, output_records

def parseVulnerability(vuln_file):
    """
    Return a dict of building classes, with each value another dict
    containing the alpha and beta values for the vulnerability model
    for that class of building.
    """

    LOG.info('Reading vulnerability parameters from {0}'.format(
                    abspath(vuln_file)))
    moddate = flModDate(vuln_file)
    LOG.info('Last modified: {0}'.format(moddate))

    up_vuln_type = []
    mu = []
    sigma = []
    scale = []
    params = AvDict()
    with open(vuln_file, 'rb') as csvfile:
        reader = csv.reader(csvfile, quotechar='"')
        for row in reader:
            if reader.line_num == 1:
                pass
            else:
                try:
                    params[str(row[0])] = {'mu':float(row[1]),
                                           'sigma':float(row[2]),
                                           'scale':float(row[3]),
                                           'sl_mu':float(row[4]),
                                           'sl_sd':float(row[5]),
                                           'sl_scale':float(row[6]),
                                           'mod_mu':float(row[7]),
                                           'mod_sd':float(row[8]),
                                           'mod_scale':float(row[9]),
                                           'ext_mu':float(row[10]),
                                           'ext_sd':float(row[11]),
                                           'ext_scale':float(row[12]),
                                           'c_mu':float(row[13]),
                                           'c_sd':float(row[14]),
                                           'c_scale':float(row[15])}

                except ValueError:
                    # Catch the case of missing values:
                    params[str(row[0])] ={'mu':10000.0,
                                           'sigma':0.1,
                                           'scale':0.0,
                                           'sl_mu':10000.,
                                           'sl_sd':0.1,
                                           'sl_scale':0.0,
                                           'mod_mu':10000.,
                                           'mod_sd':0.1,
                                           'mod_scale':0.0,
                                           'ext_mu':10000.,
                                           'ext_sd':0.1,
                                           'ext_scale':0.0,
                                           'c_mu':10000.,
                                           'c_sd':0.1,
                                           'c_scale':0.0}
                    
    LOG.debug("Vulnerability information")
    LOG.debug("Class : mu : sigma : scale ")
    for key, value in params.items():
        LOG.debug(("{0:<6}: {1:<5.3f} : {2:<5.2f} : {3:<6.4f} : "
                   "{4:<5.3f} : {5:<5.2f} : {6:<6.4f} : "
                   "{7:<5.3f} : {8:<5.2f} : {9:<6.4f} : "
                   "{10:<5.3f} : {11:<5.2f} : {12:<6.4f} : "
                   "{13:<5.3f} : {14:<5.2f} : {15:<6.4f}" ).format(
                       key, value['mu'], value['sigma'], value['scale'],
                       value['sl_mu'], value['sl_sd'], value['sl_scale'],
                       value['mod_mu'], value['mod_sd'], value['mod_scale'],
                       value['ext_mu'], value['ext_sd'], value['ext_scale'],
                       value['c_mu'], value['c_sd'], value['c_scale']))

    return params

def calculateAverageLoss(records, fields, output_folder):
    """
    Calculate average loss (for return periods and annualised loss)
    for the region.
    """

    fh = open(pjoin(output_folder,'annual_loss.csv'), 'w')

    fh.write( "Return period, average loss, cost\n" )

    loss = getField('loss', fields, records)
    cost = getField('cost', fields, records)
    values = getField('bldg_value', fields, records)
    avg_loss = np.sum(loss*cost)/np.sum(cost)

    fh.write("Average loss: {0:f}\n".format(avg_loss))
    fh.write("Cost: P{0:,} (total value: P{1:,})".
             format(int(np.sum(cost)), int(np.sum(values))))
    fh.close()

    return

def plotResults(avg_loss, tot_cost, return_periods, output_folder):
    """
    Plot the results as PML curves or similar
    """

    output_path = pjoin(output_folder, 'plots')
    
    return_periods = np.concatenate([[1], return_periods])
    probs = probability(return_periods)

    avg_loss = np.concatenate([[0], avg_loss])
    tot_cost = np.concatenate([[0], tot_cost])

    LOG.info("Plotting probability-loss curve")
    fig1 = pyplot.figure()
    ax1 = fig1.add_subplot(111)
    ax1.fill_between(avg_loss*100., probs, 0.0, lw=2, alpha=0.5)
    ax1.set_yscale('log', nonposy='clip')
    ax1.set_ylim(0.0001, 1.0)
    ax1.grid(True, which='both', linestyle=':', color='k')
    pyplot.xlabel("Damage ratio (% of reconstruction cost)")
    pyplot.ylabel("Annual probability")
    pyplot.title("Probability-loss curve")
    pyplot.savefig(pjoin(output_path, "probability-loss.png"))
    LOG.info("Plotting return period-loss curve")
    pyplot.figure(2)
    pyplot.semilogx(return_periods, avg_loss*100., lw=2, marker='+', ms=10)
    pyplot.xlabel("Return period (years)")
    pyplot.ylabel("Damage ratio (% of reconstruction cost)")
    pyplot.grid(True, which='both', linestyle=':', color='k')
    pyplot.title("Return period loss curve")
    pyplot.savefig(pjoin(output_path, "rp-loss-curve.png"))
    LOG.info("Plotting probability-cost curve")
    fig3 = pyplot.figure(3)
    ax3 = fig3.add_subplot(111)
    ax3.fill_between(tot_cost/np.power(10., 9), probs, 0.0, lw=2, alpha=0.5)
    ax3.set_yscale('log', nonposy='clip')
    ax3.set_ylim(0.0001, 1.0)
    ax3.grid(True, which='both', linestyle=':', color='k')
    pyplot.xlabel("Loss (billion Peso)")
    pyplot.ylabel("Annual probability")
    pyplot.title("Probability-loss curve")
    pyplot.savefig(pjoin(output_path, "probability-cost.png"))

    return

#@timer
def main():
    """
    Main section of the script - process command line arguments and call
    other functions to process the data
    """

    flStartLog(log_file=flConfigFile('.log'), log_level='DEBUG',
               verbose=True, datestamp=True)
    LOG.info("Parsing command line arguments")
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--costs',
                        help=('csv format file containing the cost data '
                              'for building types and land-use groupings'))
    parser.add_argument('-o', '--outputpath',
                        help=('Path to folder for storing the output' ))
    parser.add_argument('-s', '--shapefile',
                        help=('Path (including extension) of the shape '
                              'file holding the zone features to process '
                              '(e.g. land-use parcels, meshblocks, etc.)'))
    parser.add_argument('-v', '--vulnerability',
                        help=('csv format file containing the mean, sigma '
                              'and scale values for building '
                              'vulnerability curves'))
    
    args = parser.parse_args()
    cost_file = args.costs
    output_path = args.outputpath
    shape_file = args.shapefile
    vulnerability_file = args.vulnerability

    # Name of the output shape file, with no extension:
    if not isdir(abspath(output_path)):
        try:
            os.makedirs(abspath(output_path))
        except:
            print "Cannot create output path: {0}".format(output_path)
            raise
            
    if not isdir(abspath(pjoin(output_path, 'plots'))):
        try:
            os.makedirs(pjoin(abspath(output_path), 'plots'))
        except:
            print "Cannot create output path: {0}".format(pjoin(output_path,
                                                                'plots'))
            raise    

    # Add a datestring to the output file name so we know when it was created
    curdate = datetime.now()
    curdatestr = curdate.strftime('%Y%m%d%H%M')
    
    output_file = pjoin(abspath(output_path), 
                        "event_loss_{0}".format(curdatestr))
    

    # Load the exposure file
    shapes, fields, records = parseShapefile(shape_file)

    output_fields, output_records = processDamage(fields, records, vulnerability_file,
                                                  cost_file)

    # Write out the data to another shapefile
    writeShapefile(output_file, output_fields, shapes, output_records)

    damage_states = ['slight', 'moderate', 'extensive', 'complete']
    for state in damage_states:
        output_fields, output_records = processFragility(fields, records,
                                                         vulnerability_file,
                                                         state)

    output_fields, output_records = calculatePopulation(output_fields,
                                                        output_records,
                                                        vulnerability_file)
    output_file = pjoin(abspath(output_path), 
                        "event_damage_states_{0}".format(curdatestr))
    writeShapefile(output_file, output_fields, shapes, output_records)
    # Calculate average loss across the entire region contained in the
    # input shapefile:
    calculateAverageLoss(output_records, output_fields, output_path)

#*******************************************************************************
# NOW START MAIN CODE
if __name__ == "__main__":

    if len(sys.argv) == 1:
        showSyntax()
            
    try:
        main()
    except Exception: #pylint: disable-msg=W0703
        tblines = traceback.format_exc().splitlines()
        for line in tblines:
            LOG.critical(line.lstrip())
        sys.exit(-1)

    LOG.info("Completed {0}".format(sys.argv[0]))




