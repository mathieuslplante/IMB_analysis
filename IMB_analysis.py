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


from TimeUtil import TimeUtil
from Icepackdata import IcepackData
from Icepackdata import IcepackDatasetup

from SAMSIMBdata import SAMSIMB_DAsetup
from SAMSIMBdata import SAMSIMBdata
from SAMSIMBdata import SAMSIMBheatdata

from EC_WStation_data import ECStation_setup
from EC_WStation_data import ECStationData

from visualisation import visuals as vis

#------------------------------------------------------------
# Loading the IMB data
#------------------------------------------------------------

OutputFolder = 'Outputs/'
visuals = vis()

config_IMB1 = configparser.ConfigParser()
config_IMB1.read('Buoy_data/namelist_IMB.ini')

BuoySetup = SAMSIMB_DAsetup(config=config_IMB1['IMB_setup'])
timeIMB1 = TimeUtil(config = config_IMB1['Time'])
IMB1Data = SAMSIMBdata(ExpSetup=BuoySetup, time = timeIMB1)
IMB1_heat = SAMSIMBheatdata(ExpSetup=BuoySetup, time = timeIMB1)

Rtrvl_IMB1 = SfcRetrieval(config=config_IMB1['Retrieval'],Data = IMB1Data.data, Buoy = BuoySetup)
Rtrvl_IMB1.compute_inferfaces_minimisation(Data= IMB1Data.data,Buoy= BuoySetup)
Rtrvl_IMB1.PrintMassBalance(Data = IMB1Data, OutputFolder = OutputFolder +'IMB_MassBalance')

#End of script
