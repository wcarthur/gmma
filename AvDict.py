#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Title: AvDict.py - autovivication of dicts for python
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: Fri May 03 11:39:12 2013
 Description:

 Version:
 Id:

"""

class AvDict(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
