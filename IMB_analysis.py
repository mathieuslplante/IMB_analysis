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
from DAL_IMBdata import DAL_IMBdata

from EC_WStation_data import ECStation_setup
from EC_WStation_data import ECStationData

from visualisation import visuals as vis

#------------------------------------------------------------
# Loading the IMB data
#------------------------------------------------------------

OutputFolder = 'Outputs/'
visuals = vis()

config_IMB1 = configparser.ConfigParser()
config_IMB1.read('Buoy_data/namelist_IMB_DAL.ini')

BuoySetup = SAMSIMB_DAsetup(config=config_IMB1['IMB_setup'])

timeIMB1 = TimeUtil(config = config_IMB1['Time'])

IMB1Data = DAL_IMBdata(ExpSetup=BuoySetup, time = timeIMB1)

Rtrvl_IMB1 = SfcRetrieval(config=config_IMB1['Retrieval'],Data = IMB1Data.data, Buoy = BuoySetup)
Rtrvl_IMB1.compute_inferfaces_minimisation(Data= IMB1Data.data,Buoy= BuoySetup)
Rtrvl_IMB1.PrintMassBalance(Data = IMB1Data, OutputFolder = OutputFolder +'IMB_MassBalance')

# Figure 3
visuals.figure_temperatures_weather(Data = IMB1Data,
                                Rtrvl = Rtrvl_IMB1,
                                krun = 0,
                                OutputFolder = OutputFolder)

# Figure 5
visuals.figure_thickness_obs(data_buoy = IMB1Data,
                                       rtrvl = Rtrvl_IMB1,
                                       krun = 0,
                                       OutputFolder = OutputFolder)

#End of script
