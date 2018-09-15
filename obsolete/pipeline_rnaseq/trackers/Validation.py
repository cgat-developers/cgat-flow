import os
import sys
import re
import types
import itertools

from cgatReport.Tracker import *
from RnaseqReport import *


class ExonValidationSummary(RnaseqTracker, SingleTableTrackerRows):
    table = "exon_validation"
