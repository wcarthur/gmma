#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Title: AggUtils.py - aggregation utilities for GMMA RAP
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: Thu May 02 13:09:26 2013
 Description:

 Version: 0.1
 Id:

"""

import os
import csv
import logging
from os.path import join as pjoin, abspath

import numpy as np
import matplotlib
matplotlib.use('Agg', warn=False)

from matplotlib import pyplot
from get_records import getField
from probability import probability
from AvDict import AvDict

LOG = logging.getLogger(__name__)

def extractZoneValues(featureid, features, records):
    """
    Extract values from a set of records based on a matched feature id field.
    The 'featureid' variable contains the unique identifier for each feature in
    the dataset.
    'features' is a dict keyed by all the available feature id numbers (which
    should contain all those in 'featureid'). This is loaded from the complete
    list of feature id numbers, using load_zones_from_records()
    'records' is the array of values and should be the same length as
    'featureid'

    """
    regions = {}
    subregions = {}

    for feature, record in zip(featureid, records):
        region = features[int(feature)]['region']
        subregion = features[int(feature)]['subregion']
        sub_key = '{0} - {1}'.format(region, subregion)
        if regions.has_key(region):
            regions[region] = np.append(regions[region], record)
        else:
            regions[region] = np.array([record])

        if subregions.has_key(sub_key):
            subregions[sub_key] = np.append(subregions[sub_key], record)
        else:
            subregions[sub_key] = np.array([record])

    return regions, subregions

def getRegions(zonefile):
    """
    Obtain a list of all unique regions and sub regions contained in the input
    zonefile
    """
    LOG.debug("Extracting list of unique regions from {0}".format(
                                              abspath(zonefile)))

    regions = ['UNDEFINED']
    subregions = ['UNDEFINED']
    with open(os.path.abspath(zonefile), 'r') as filehandle:
        reader = csv.reader(filehandle, quotechar='"')
        for row in reader:
            if reader.line_num == 1:
                pass
            else:
                [oid, zone] = row
                try:
                    region = zone.split(' - ', 1)[0]
                except ValueError:
                    pass

                if region not in regions:
                    regions.append(region)
                # Names of subregions may be shared across multiple regions
                # so we capture the combination of region - subregion name
                if zone not in subregions:
                    subregions.append(zone)

    LOG.debug("Found {0:d} unique regions".format(len(regions)))
    LOG.debug("Found {0:d} unique subregions".format(len(subregions)))

    return regions, subregions

def loadZonesFromRecords(records, fields, featureid, zoneid):
    """
    Obtain a list of all unique regions and sub regions contained in the input
    shapefile
    """
    LOG.debug("Extracting list of unique regions from records")

    output = dict()

    featurerecs = getField(featureid, fields, records, dtype=int)
    zonerecs = getField(zoneid, fields, records, dtype=str)

    for feature, zone in zip(featurerecs, zonerecs):
        try:
            region, subregion = zone.split(' - ', 1)
        except ValueError:
            region = zone
            subregion = zone
        finally:
            output[feature] = {'region':region,
                               'subregion':subregion}

    return output


def aggregate_events(records, fields, featureid, features):
    """
    'records' - all records contained in the shapefile
    'fields' - list of fields held in the shape file
    'zoneid' - name of the unique feature id field in the records
    'features' is a dict keyed by all the available feature id numbers (which
    should contain all those in 'featureid'). This is loaded from the complete
    list of feature id numbers, using load_zones_from_records()

    """
    LOG.info("Calculating aggregated losses")
    # Extract from teh shape file the required values
    oid = getField(featureid, fields, records, dtype=int)
    value = getField('bldg_value', fields, records)
    flarea = getField('FLAREA_SUM', fields, records)
    polygon_area = getField('AREA_SQM', fields, records)
    aac = getField('cost', fields, records)
    popn = getField('POP_EST', fields, records)
    popa = getField('POP_AFFECT', fields, records, dtype=float)

    # Aggregate the values to region and subregion:
    r_value, sr_value = extractZoneValues(oid, features, value)
    r_aac, sr_aac = extractZoneValues(oid, features, aac)
    r_area, sr_area = extractZoneValues(oid, features, polygon_area)
    r_flarea, sr_flarea = extractZoneValues(oid, features, flarea)
    r_popn, sr_popn = extractZoneValues(oid, features, popn)
    r_popa, sr_popa = extractZoneValues(oid, features, popa)

    output = AvDict()
    sroutput = AvDict()

    for region in r_value.keys():
        region_value = np.sum(r_value[region])
        # Floor area of the region in hectares:
        region_flarea = np.sum(r_flarea[region])/np.power(10., 4)
        # Total area of the region in square kilometres:
        region_area = np.sum(r_area[region])/np.power(10., 6)

        # Total population in the region:
        region_pop = np.sum(r_popn[region])

        output[region]['VALUE'] = region_value
        output[region]['flarea_sum'] = region_flarea
        output[region]['area_sum'] = region_area
        output[region]['POP_EST'] = region_pop

        output[region]['cost'] = np.sum(r_aac[region])
        output[region]['loss'] = output[region]['cost']/region_value
        output[region]['dmg'] = output[region]['loss'] * \
                                            region_flarea/region_area
        output[region]['POP_AFFECT'] = np.nansum(r_popa[region])

    for region in sr_value.keys():
        sregion_value = np.sum(sr_value[region])
        # Floor area of the region in hectares:
        sregion_flarea = np.sum(sr_flarea[region])/np.power(10., 4)
        # Total area of the region in square kilometres:
        sregion_area = np.sum(sr_area[region])/np.power(10., 6)

        # Total population in the region:
        sregion_pop = np.sum(sr_popn[region])

        sroutput[region]['VALUE'] = sregion_value
        sroutput[region]['flarea_sum'] = sregion_flarea
        sroutput[region]['area_sum'] = sregion_area
        sroutput[region]['POP_EST'] = sregion_pop

        sroutput[region]['cost'] = np.sum(sr_aac[region])
        sroutput[region]['loss'] = sroutput[region]['cost']/sregion_value
        sroutput[region]['dmg'] = sroutput[region]['loss'] * \
                                            sregion_flarea/sregion_area
        sroutput[region]['POP_AFFECT'] = np.nansum(sr_popa[region])

    return output, sroutput

def aggregate(records, fields, featureid, features, return_periods):
    """
    'records' - all records contained in the shapefile
    'fields' - list of fields held in the shape file
    'zoneid' - name of the unique feature id field in the records
    'features' is a dict keyed by all the available feature id numbers (which
    should contain all those in 'featureid'). This is loaded from the complete
    list of feature id numbers, using load_zones_from_records()
    'return_periods' - array of return period values for which to calculate
        aggregated losses.
    """
    LOG.info("Calculating aggregated losses")

    # Extract from teh shape file the required values
    oid = getField(featureid, fields, records, dtype=int)
    value = getField('bldg_value', fields, records)
    flarea = getField('FLAREA_SUM', fields, records)
    polygon_area = getField('AREA_SQM', fields, records)
    aac = getField('ann_cost', fields, records)

    # Aggregate the values to region and subregion:
    r_value, sr_value = extractZoneValues(oid, features, value)
    r_aac, sr_aac = extractZoneValues(oid, features, aac)
    r_area, sr_area = extractZoneValues(oid, features, polygon_area)
    r_flarea, sr_flarea = extractZoneValues(oid, features, flarea)
    output = AvDict()
    sroutput = AvDict()

    for region in r_value.keys():
        region_value = np.sum(r_value[region])
        # Floor area of the region in hectares:
        region_flarea = np.sum(r_flarea[region])/np.power(10., 4)
        # Total area of the region in square kilometres:
        region_area = np.sum(r_area[region])/np.power(10., 6)

        output[region]['VALUE'] = region_value
        output[region]['flarea_sum'] = region_flarea
        output[region]['area_sum'] = region_area

        output[region]['ann_cost'] = np.sum(r_aac[region])
        output[region]['ann_loss'] = output[region]['ann_cost']/region_value
        output[region]['ann_dmg'] = output[region]['ann_loss'] * \
                                            region_flarea/region_area

    LOG.debug("Calculating zonal totals for each return period")
    for ret_per in return_periods:
        dmg_key = 'dmg' + str(int(ret_per))
        cost_key = 'cost' + str(int(ret_per))
        flarea_key = 'flarea' + str(int(ret_per))

        cost = getField(cost_key, fields, records)
        r_cost, sr_cost = extractZoneValues(oid, features, cost)

        for region in r_value.keys():
            cost_sum = np.sum(r_cost[region])
            dmg_sum = cost_sum/np.sum(r_value[region])
            output[region][dmg_key] = dmg_sum
            output[region][flarea_key] = dmg_sum * \
                                         output[region]['flarea_sum'] / \
                                         output[region]['area_sum']

            output[region][cost_key] = cost_sum

    for region in sr_value.keys():
        region_value = np.sum(sr_value[region])
        # Floor area of the region in hectares:
        region_flarea = np.sum(sr_flarea[region])/np.power(10., 4)
        # Total area of the region in square kilometres:
        region_area = np.sum(sr_area[region])/np.power(10., 6)

        sroutput[region]['VALUE'] = region_value
        sroutput[region]['ann_cost'] = np.sum(sr_aac[region])
        sroutput[region]['ann_loss'] = sroutput[region]['ann_cost'] / \
                                        region_value
        sroutput[region]['ann_dmg'] = sroutput[region]['ann_loss'] * \
                                                region_flarea/region_area
        sroutput[region]['flarea_sum'] = region_flarea
        sroutput[region]['area_sum'] = region_area

    for ret_per in return_periods:
        dmg_key = 'dmg' + str(int(ret_per))
        flarea_key = 'flarea' + str(int(ret_per))
        cost_key = 'cost' + str(int(ret_per))
        cost = getField(cost_key, fields, records)
        r_cost, sr_cost = extractZoneValues(oid, features, cost)

        for region in sr_value.keys():
            cost_sum = np.sum(sr_cost[region])
            dmg_sum = cost_sum/np.sum(sr_value[region])
            sroutput[region][dmg_key] = dmg_sum
            sroutput[region][flarea_key] = dmg_sum * \
                                        sroutput[region]['flarea_sum'] / \
                                        sroutput[region]['area_sum']
            sroutput[region][cost_key] = cost_sum

    return output, sroutput

def plotLossOutput(zonename, loss, probs, output_path):
    """
    Plot loss results against probability
    """
    fig = pyplot.figure()
    axis = fig.add_subplot(111)
    if zonename == '':
        filename = pjoin(os.path.abspath(output_path),
                                'unassigned.prob-loss.png')
    else:
        filename = pjoin(os.path.abspath(output_path),
                                "%s.prob-loss.png" % zonename)
    axis.fill_between(loss*100., probs, 0.0, lw=2, alpha=0.5)
    axis.grid(True, which='both', linestyle=':', color='k')
    axis.set_yscale('log', nonposy='clip')
    axis.set_ylim(0.0001, 1.0)
    pyplot.xlabel("Loss ratio\n(% of reconstruction cost)")
    pyplot.ylabel("Annual probability")
    pyplot.title("Probability-loss curve for %s" % zonename)
    pyplot.savefig(filename)

def plotCostOutput(zonename, costs, probs, output_path):
    """
    Plot cost results against probability
    """
    fig = pyplot.figure()
    axis = fig.add_subplot(111)
    if zonename == '':
        filename = pjoin(os.path.abspath(output_path),
                                'unassigned.prob-cost.png')
    else:
        filename = pjoin(os.path.abspath(output_path),
                                "%s.prob-cost.png" % zonename)
    axis.fill_between(costs/1000000., probs, 0.0, lw=2, alpha=0.5)
    axis.grid(True, which='both', linestyle=':', color='k')
    axis.set_yscale('log', nonposy='clip')
    axis.set_ylim(0.0001, 1.0)

    pyplot.xlabel("Cost (Million Peso)")
    pyplot.ylabel("Annual probability")
    pyplot.title("Probability-cost curve for %s" % zonename)
    pyplot.savefig(filename)

def plotDmgOutput(zonename, dmg, probs, output_path):
    """
    Plot damaged floor area equivalent results against probability
    """
    fig = pyplot.figure()
    axis = fig.add_subplot(111)
    if zonename == '':
        filename = pjoin(os.path.abspath(output_path),
                                'unassigned.prob-cost.png')
    else:
        filename = pjoin(os.path.abspath(output_path),
                                "%s.prob-dmg.png" % zonename)
    axis.fill_between(dmg, probs, 0.0, lw=2, alpha=0.5)
    axis.grid(True, which='both', linestyle=':', color='k')
    axis.set_yscale('log', nonposy='clip')
    axis.set_ylim(0.0001, 1.0)

    pyplot.xlabel("Damage (ha/km^2)")
    pyplot.ylabel("Annual probability")
    pyplot.title("Probability-damage curve for %s" % zonename)
    pyplot.savefig(filename)

def writeEventOutput(output_file, data):
    """
    Write aggregated impact information for a single event to file
    """
    
    output_fileh = open(output_file, 'w')
    zones = data.keys()
    header = "REGION,POP_EST,AREA_SQM,FLAREA_SUM,TOTAL_VALUE,LOSS,COST,DMGF,POP_AFFECT\n"
    fmt = "%s, %f, %f, %f, %f, %f, %f, %f, %d\n"
    output_fileh.write(header)
    for zone in zones:
        output_fileh.write(fmt % (zone, data[zone]['POP_EST'],
                                        data[zone]['area_sum'],
                                        data[zone]['flarea_sum'],
                                        data[zone]['VALUE'],
                                        data[zone]['loss'],
                                        data[zone]['cost'],
                                        data[zone]['dmg'],
                                        data[zone]['POP_AFFECT']))
    output_fileh.close()

def writeLossOutput(output_file, data, return_periods, plot=True):
    """
    Write aggregated loss information to file and plot probability-loss
    curves for each zone
    """
    LOG.info("Writing loss data to {0}".format(output_file))
    rps = np.concatenate([[1], return_periods])
    probs = probability(rps)

    output_path = os.path.split(output_file)[0]
    output_fileh = open(output_file, 'w')
    zones = data.keys()
    header = "REGION," + \
             ",".join(["RP" + str(int(r)) for r in return_periods]) + \
             ",ANN_LOSS\n"
    output_fileh.write(header)
    for zone in zones:
        row = [data[zone]["dmg"+str(int(r))] for r in return_periods]

        output_fileh.write( "%s, %s, %7.5f\n" %
            (zone, ", ".join(["%7.5f" % val for val in row]),
             data[zone]['ann_loss'] ) )

        row = np.concatenate([[0], row])
        if plot:
            plotLossOutput(zone, row, probs, pjoin(output_path, 'plots'))

    output_fileh.close()

def writeDmgOutput(output_file, data, return_periods, plot=True):
    """
    Write aggregated loss information to file and plot probability-loss
    curves for each zone
    """
    LOG.info("Writing damage data to {0}".format(output_file))
    rps = np.concatenate([[1], return_periods])
    probs = probability(rps)

    output_path = os.path.split(output_file)[0]
    output_fileh = open(output_file, 'w')
    zones = data.keys()
    header = "REGION," + \
             ",".join(["RP" + str(int(r)) for r in return_periods]) + \
             ",ANN_DMG, TOTAL_FLAREA\n"
    output_fileh.write(header)
    for zone in zones:
        row = [data[zone]["flarea"+str(int(r))] for r in return_periods]
        output_fileh.write( "%s, %s, %7.5f, %9.5f\n" %
            (zone, ", ".join(["%7.5f" % val for val in row]),
             data[zone]['ann_dmg'], data[zone]['flarea_sum'] ) )

        row = np.concatenate([[0], row])

        if plot:
            plotDmgOutput(zone, row, probs, pjoin(output_path, 'plots'))

    output_fileh.close()

def writeCostOutput(output_file, data, return_periods, plot=True):
    """
    Write aggregated cost information to file and plot probability-cost
    curves for each zone.
    """
    LOG.info("Writing cost data to {0}".format(output_file))
    rps = np.concatenate([[1], return_periods])
    probs = probability(rps)

    output_path = os.path.split(output_file)[0]
    output_fileh = open(output_file, 'w')
    zones = data.keys()
    header = "REGION," + \
             ",".join(["RP" + str(int(r)) for r in return_periods]) + \
             ",ANN_COST,TOTAL_VALUE\n"
    output_fileh.write(header)

    for zone in zones:
        row = [data[zone]["cost"+str(int(r))] for r in return_periods]

        output_fileh.write( "%s, %s, %12.0f, %13.0f\n" %
           (zone, ", ".join(["%7.5f" % val for val in row]),
            data[zone]['ann_cost'], data[zone]['VALUE'] ) )

        row = np.concatenate([[0], row])
        if plot:
            plotCostOutput(zone, row, probs, pjoin(output_path, 'plots'))

    output_fileh.close()
    return
