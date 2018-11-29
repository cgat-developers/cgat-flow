from PeakcallingReport import *
from cgatReport.Tracker import *


class BamSummary(CallingTracker, SingleTableTrackerRows):
    table = "bam_stats"
