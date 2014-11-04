#!/usr/bin/env python
"""
 Title: damage.py - calculate damage using lognormal CDF
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2013-04-29
 Description: Calculate the damage level based on a given array of wind speed
    values and a given mu and sigma value. The mu and sigma control
    the form of the vulnerability function and are specific to
    each building class. This version uses the log-normal cumulative
    probability function to describe the damage level.

    Wind speed values are stored in metres/second, but the vulnerability
    relations are based on km/h, so we convert on the fly here.

    The scale value 's' is used to reduce the total damage for a building
    type to provide an upper limit to damage (e.g. assumed only damage
    is to windows/cladding)

 Version: 0.1

 (C) Commonwealth of Australia (Geoscience Australia) 2012
 This product is released under the Creative Commons Attribution 3.0
 Australia Licence

 http://creativecommons.org/licenses/by/3.0/au/legalcode

"""

from scipy.stats import lognorm
import numpy as np
import logging

__eps__ = 1.0e-6
LOGGER = logging.getLogger(__name__)

def damage(wind_speed, mu, sigma, scale=1.0):
    """
    Calculate the damage level based on a given array of wind speed
    values and a given mu and sigma value. The mu and sigma control
    the form of the vulnerability function and are specific to
    each building class. This version uses the log-normal cumulative
    probability function to describe the damage level.

    Wind speed values are stored in metres/second, but the vulnerability
    relations are based on km/h, so we convert on the fly here.

    The scale value 's' is used to reduce the total damage for a building
    type to provide an upper limit to damage (e.g. assumed only damage
    is to windows/cladding)
    """

    # mu is the scale parameter, sigma is the shape parameter
    # of the log-normal distribution:
    #if mu == 0.0:
    dmg = np.zeros(len(wind_speed))
    dmg = scale*lognorm.cdf(wind_speed*3.6, sigma, scale=mu)
    # Mask 'small' damage values to be zero:
    np.putmask(dmg, dmg < __eps__, 0.0)
    np.putmask(dmg, mu==0.0, 0.0)

    return dmg

def adjust_curves(bld_type, vmask, mu, sigma, scale):
    """
    Change the vulnerability curve for those building types
    where there are different materials used in different
    eras. NOTE: THESE ARE HARD-CODED VALUES. ANY CHANGE TO
    THE VULNERABILITY CURVES WILL REQUIRE CHANGES TO BE MADE
    IN THIS SECTION
    """

    LOGGER.debug( ("Adjusting curves for building types with"
                    "age-dependent vulnerability") )
    if bld_type.startswith('MWS'):
        np.putmask(mu, vmask, 398.)
        np.putmask(sigma, vmask, 0.16)
    elif bld_type.startswith('CHB'):
        np.putmask(mu, vmask, 469.)
        np.putmask(sigma, vmask, 0.44)
    elif bld_type.startswith('CWS'):
        np.putmask(mu, vmask, 398.)
        np.putmask(sigma, vmask, 0.16)
    elif bld_type.startswith('C1_L'):
        np.putmask(mu, vmask, 477.)
        np.putmask(sigma, vmask, 0.46)
    elif bld_type.startswith('C2_L'):
        np.putmask(mu, vmask, 477.)
        np.putmask(sigma, vmask, 0.46)
    elif bld_type.startswith('C4_L'):
        np.putmask(mu, vmask, 477.)
        np.putmask(sigma, vmask, 0.46)
    elif bld_type.startswith('PC1_L'):
        np.putmask(mu, vmask, 477.)
        np.putmask(sigma, vmask, 0.46)
    elif bld_type.startswith('PC2_L'):
        np.putmask(mu, vmask, 477.)
        np.putmask(sigma, vmask, 0.46)

    return mu, sigma, scale
