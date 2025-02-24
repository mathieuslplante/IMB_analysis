
import shutil
import os
import sys
import array
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, date, time, timedelta
import urllib.request as request
from contextlib import closing
import cmocean as cmo
from scipy import stats


import cartopy.crs as ccrs
import cartopy.feature

import shapely.geometry as sgeom
import cmocean as cmo

class DAL_IMBdata:
    def __init__(self,
                 ExpSetup=None, time = None, url = None, OutputFolder=None,StartTime = None,loc = None):

        #-----------------------------------
        # Load the data from the input file
        #-----------------------------------

        data = np.genfromtxt(ExpSetup.DataFile,delimiter=',',dtype=None,names=True)
        header = data.dtype.names

        #--------------------------
        #Find the indice corresponding to required columns
        #--------------------------

        a = len(data[ExpSetup.Sensor0_Header])
        b = len(header)


        indt0 = self.find_unstruct(header,ExpSetup.Sensor0_Header)
        indtlast = self.find_unstruct(header,ExpSetup.SensorLast_Header)
        indtair = self.find_unstruct(header,ExpSetup.Tair_Header)
        indtwater = self.find_unstruct(header,ExpSetup.Tocean_Header)

        date_all = data[ExpSetup.Time_Header][:]
        time_all = data[str(int(ExpSetup.Time_Header)+1)][:]
        reference_date = datetime(2022,1,1,hour=0)

        #--------------------------
        # Get the number of time steps and start row (kinit)
        #--------------------------

        nstep = time.nstep
        kinit = 0
        print(ExpSetup.Time_Header)
        print(data[str(int(ExpSetup.Time_Header)+1)][:])

        date_k = data[ExpSetup.Time_Header][kinit]
        time_k = data[str(int(ExpSetup.Time_Header)+1)][kinit]
        #finding the date of the first data, based on the
        #reference (Jan 01 2022, 0Z)
        ThisDate = self.extract_time(date_data = date_k, time_data = time_k)
        #Get the row corresponding to the start time
        while ThisDate < time.StartDate:
            kinit = kinit+1
            date_k = data[ExpSetup.Time_Header][kinit]
            time_k = data[str(int(ExpSetup.Time_Header)+1)][kinit]
            ThisDate = self.extract_time(date_data = date_k, time_data = time_k)

        self.tstep = ExpSetup.tstep
        print(kinit, nstep)

        #--------------------------
        # Prepare datasets
        #--------------------------

        self.Lat = np.zeros(nstep)
        self.Lat[:] = np.nan
        self.Lon = np.zeros(nstep)
        self.Lon[:] = np.nan
        self.airT = np.zeros(nstep)
        self.airT[:] = np.nan
        self.oceT = np.zeros(nstep)
        self.oceT[:] = np.nan
        self.datelist = np.zeros(nstep)
        self.datelist[:] = np.nan
        self.t = np.zeros(nstep)
        self.t[:] = np.nan
        self.data = np.zeros((nstep,ExpSetup.nsensors+1))
        self.data[:,:] = np.nan

        #initializing memory variables of previous date and hour to impossible values
        prev_date = reference_date - timedelta(days=10000)
        prev_hour = -1


        #--------------------------
        # Fetch hourly data
        #   Note: here, we re-organise the data into 1hour
        #   bins, so that we can make a pcolor. This means
        #   that we are skipping few duplicate transmissions
        #--------------------------
        for kdata in range(kinit,a):

            #Update the date object
            date_k = data[ExpSetup.Time_Header][kdata]
            time_k = data[str(int(ExpSetup.Time_Header)+1)][kdata]
            date_k = self.extract_time(date_data = date_k, time_data = time_k)

            #record the first date as reference t=0
            if kdata == kinit:
                ref_time = date_k

            #Attributing the row data to the right time step by rounding the time to nearest hour
            #This deals with double data, or shuffled data due to bad communications.
            new_date = date_k
            new_hour = (date_k - ref_time).days*24 +round((date_k - ref_time).seconds/(60*60))
            kout = int(new_hour/ExpSetup.tstep)

            #Standard time in days since 01-01-1900.
            #44562.0 is 2022-01-01 at 00.
            new_time = 44562.0 +(date_k - reference_date).days + ((date_k - reference_date).seconds /(24*60*60))
            print(new_time)


            #----------------------------------------------
            # Extract the data
            #----------------------------------------------
            if kout < nstep:
                self.airT[kout]=data[header[indtair]][kdata]
                for indname in range(indt0, indtlast+1):
                    self.data[kout,indname-indt0+1]= data[header[indname]][kdata]
                    self.datelist[kout] =  new_date.strftime("%Y%m%d%H")
                    self.t[kout]= new_time
                    self.oceT[kout]=data[header[indtwater]][kdata]

    def find_unstruct(self,header, name1):
        b = len(header)
        for kname in range(0,b):
            if header[kname]==name1:
                indname = kname
        return indname


    def extract_time(self,date_data = None, time_data = None):
        year_data = int(date_data[0:4])
        month_data = int(date_data[5:7])
        day_data = int(date_data[8:10])
        hour_data = int(time_data[0:2])
        time_k = time(hour_data,0,0)
        date_k = date(year_data,month_data,day_data)
        date_time_k = datetime.combine(date_k,time_k)

        return date_time_k

