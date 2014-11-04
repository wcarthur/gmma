#!/usr/bin/env python
"""
 Title: get_records.py - extract a given record from a series of records
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2013-04-29
 Description: Extracts values from an array of lists that each contain data
              for a specified field.

 Version: 0.1

 (C) Commonwealth of Australia (Geoscience Australia) 2012
 This product is released under the Creative Commons Attribution 3.0
 Australia Licence

 http://creativecommons.org/licenses/by/3.0/au/legalcode
"""
import numpy
import logging

__logger__ = logging.getLogger(__name__)

def get_field(fieldname, fields, records, dtype=float):
    """
    Extract from the records the value of the field corresponding
    to fieldname.
    """
    __logger__.debug("Extracting {0} from records".format(fieldname))
    nrecords = len(records)
    fieldnames = [fields[i][0] for i in range(len(fields))]

    # Check the field is in the records:
    if fieldname not in fieldnames:
        __logger__.warn("No field '{0}' in the list of fieldnames" .
                    format(fieldname))

        __logger__.warn("Unable to proceed with processing")
        raise ValueError

    # Get the index of the required field name:
    idx = fieldnames.index(fieldname)
    if dtype != str:
        # For non-string data, return a numpy array:
        output = numpy.array([records[rec][idx] for rec in xrange(nrecords)],
                              dtype=dtype)

    else:
        # Otherwise, return a list:
        output = [records[rec][idx] for rec in xrange(nrecords)]


    return output

def remove_field(fieldname, fields, records):
    """
    Remove from the records the given field
    """
    __logger__.debug("Extracting {0} from records".format(fieldname))
    nrecords = len(records)
    fieldnames = [fields[i][0] for i in range(len(fields))]

    # Check the field is in the records:
    if fieldname not in fieldnames:
        __logger__.warn("No field '{0}' in the list of fieldnames" .
                    format(fieldname))

        __logger__.warn("Unable to proceed with processing")
        raise ValueError

    # Get the index of the required field name:
    idx = fieldnames.index(fieldname)
    for rec in xrange(nrecords):
        del records[rec][idx]

    del fields[idx]

    return records, fields