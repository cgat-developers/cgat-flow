import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
import cpgReport

from cgatReport.Tracker import *
from cgatReport.odict import OrderedDict as odict

##########################################################################
##########################################################################
##########################################################################


class chromatinProfileTracker(TrackerImages):

    """Chromatin profile per gene """
