#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Title: probability.py - calculate annual probability given return period
 Author: Craig Arthur
 CreationDate: Wed May 01 16:23:36 2013
 Description:

 Version:
 Id:

"""

import logging
import numpy

logger = logging.getLogger(__name__)

def probability(return_period):
    """Return an annual probability given a return period"""
    logger.debug("Calculating annual probability from return periods")
    prob = 1.0 - numpy.exp(-1.0/return_period)
    return prob