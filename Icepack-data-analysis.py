#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from datetime import datetime, date, time, timedelta
import os
import sys
import cmocean as cmo

from SfcRetrieval import SfcRetrieval
import configparser


from Icepackdata import IcepackData
from Icepackdata import IcepackDatasetup

config = configparser.ConfigParser()
config.read('./namelist_Icepack.ini')

MetaData = IcepackDatasetup(config['Icepack_setup'],config_Labels = config['Diag_infos'])
Data = IcepackData(meta = MetaData)

