
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

from math import floor,ceil,pi


import cartopy.crs as ccrs
import cartopy.feature

import shapely.geometry as sgeom
import cmocean as cmo

class ECStation_setup:
    def __init__(self,config = None):
        self.DataInPath = config['DataInPath']
        self.OutputFolder = config['OutPath']
        self.tstep = int(config['dt'])
        self.Tair_Header = config['Tair_Header'] 
        self.UUA_Header = config['UUA_Header'] 
        self.VVA_Header = config['VVA_Header'] 
        self.Time_Header = config['Time_Header']
        self.Lon_Header = config['Lon_Header'] 
        self.Lat_Header = config['Lat_Header'] 
        self.Tdew_Header = config['Tdew_Header']
        self.Precip_Header= config['Precip_Header']
        self.RH_Header= config['RH_Header']
        self.WDir_Header = config['WDir_Header']
        self.WS_Header = config['WS_Header']
        self.Press_Header = config['Press_Header'] 
        self.DataFile = self.DataInPath

            
class ECStationData:
    def __init__(self, 
                 ExpSetup=None, time = None, url = None, OutputFolder=None,StartTime = None,loc = None):


        #--------------------------
        # Prepare datasets
        #-------------------------- 
        
        self.Lat = np.array([])
        self.Lon = np.array([])
        self.airT = np.array([])
        self.UUA = []
        self.VVA = []
        self.datelist = []
        self.t = []
        self.Press = []
        self.Tdew = []
        self.Precip= np.array([])
        self.RH= []
        self.WDir = np.array([])
        self.WS = np.array([])


        reference_date = datetime(2022,1,1,hour=0)
        #-----------------------------------
        # Load the data from the input file
        #-----------------------------------
        
        n = 0

        for months in range(time.StartMonth,time.EndMonth+1):
            print(months)
            DataFile = self.get_name_files(ExpSetup,mm=months,yy=time.StartYear)
            data_m = np.genfromtxt(DataFile,delimiter=';',dtype=float,names=True)
            data_mstr = np.genfromtxt(DataFile,delimiter=';',dtype=None,names=True)
            header = data_m.dtype.names
            print(header)
            self.airT = np.append(self.airT,data_m[str(ExpSetup.Tair_Header)][:])
            if ExpSetup.UUA_Header != "None":
                self.UUA = np.append(self.UUA,data_m[ExpSetup.UUA_Header][:])
                self.VVA = np.append(self.VVA,data_m[ExpSetup.VVA_Header][:])
            if ExpSetup.Tdew_Header != "None":
                self.Tdew = np.append(self.Tdew,data_m[ExpSetup.Tdew_Header][:])

            if ExpSetup.RH_Header != "None":
                self.RH = np.append(self.RH,data_m[ExpSetup.RH_Header][:])
            if ExpSetup.WDir_Header != "None":
                self.WDir = np.append(self.WDir,data_m[ExpSetup.WDir_Header][:])
            if ExpSetup.WS_Header != "None":
                self.WS = np.append(self.WS,data_m[ExpSetup.WS_Header][:])
            if ExpSetup.Press_Header != "None":
                self.Press = np.append(self.Press,data_m[ExpSetup.Press_Header][:])
            if ExpSetup.Precip_Header != "None":
                self.Precip = np.append(self.Precip,data_m[ExpSetup.Precip_Header][:])
            ExpSetup.Lon_Header = header[0]
            self.Lat = np.append(self.Lat,data_m[ExpSetup.Lat_Header][:])
            self.Lon = np.append(self.Lon,data_m[ExpSetup.Lon_Header][:])

            self.datelist = np.append(self.datelist,data_mstr[ExpSetup.Time_Header][:])
            
        n = len(self.Lat)
        for k in range(0,n):
       #     self.Lat[k] = float(str(self.Lat[k][2:-2],'utf-8'))
       #     self.Lat[k] = float(self.Lat[k])
            print(self.WS[k])
       #     sadfgsa
       #     self.Lon[k] = float(str(self.Lon[k][1:-1],'utf-8'))
       #     self.airT[k] = float(str(self.airT[k][1:-1],'utf-8'))
       #     #self.UUA[k] = float(self.UUA[k])
            #self.VVA[k] = float(self.VVA[k])
       #     self.Press[k] = float(str(self.Press[k][1:-1],'utf-8'))
       #     self.Tdew[k] = float(str(self.Tdew[k][1:-1],'utf-8'))
       #     self.Precip[k]= float(str(self.Precip[k][1:-1],'utf-8'))
       #     #self.RH[k]= float(str(self.RH[k][1:-1],'utf-8'))
       #     try:
       #         self.WDir[k] = float(str(self.WDir[k][1:-1],'utf-8'))
       #     except:
       #         self.WDir[k] = np.nan 
       #     self.WS[k] = float(str(self.WS[k][1:-1],'utf-8'))

      
        if ExpSetup.UUA_Header  != "None":
            if ExpSetup.WS_Header == "None":
                self.WS = (self.UUA**2.0 + self.VVA**2.0)**0.5    
                self.WDir = np.arctan(self.VVA[:]/self.UUA[:])
                self.WDir[self.UUA<0.0] = self.WDir - pi

        #--------------------------
        # Get the time correspinding to each time stamps
        #--------------------------

        self.t = np.arange(n)*np.nan
        for k in range(0,n):
            date_k = self.extract_time(date_data = str(self.datelist[k]))
            self.t[k] = 44562.0 +(date_k - reference_date).days + ((date_k - reference_date).seconds /(24*60*60))

    def get_name_files(self,ExpSetup=None,mm=None, yy=None): 
        
        station = "NL"
        station_num = "8502799"
        freq = "P1H"
        month = str(mm)
        month = month.zfill(2) 
        year = str(yy)
        name_file = "%sen_climate_hourly_%s_%s_%s-%s_%s.csv" % (ExpSetup.DataInPath, station, station_num, month, year, freq)

        return name_file
                    
    def find_unstruct(self,header, name1):
        b = len(header)
        for kname in range(0,b):
            if header[kname]==name1:
                indname = kname
        return indname
    
    def extract_time(self,date_data = None):
        #print(date_data)
        year_data = int(date_data[2:6])
        #print(year_data)
        month_data = int(date_data[7:9])
        #print(month_data)
        day_data = int(date_data[10:12])
        #print(day_data)
        hour_data = int(date_data[13:15])
        #print(day_data)
        time_k = time(hour_data,0,0)
        date_k = date(year_data,month_data,day_data)
        date_time_k = datetime.combine(date_k,time_k)
        #print(date_time_k)

        return date_time_k

    def show_temperatures(self, OutputFolder = None):

        Figure = plt.figure(figsize=[6.0, 8.0])
  
        ax1 = Figure.add_axes([.2, 0.07, 0.75, 0.25])
        plt.plot(self.t,self.airT) 
        plt.ylim(-30.0,10.0)
        plt.ylabel('Sfc Temperature ($^{\circ}$)')
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        

        ax2 = Figure.add_axes([0.2, 0.7, 0.75, 0.25])
        plt.plot(self.t,self.WS) 
        #plt.plot(self.t,self.WDir) 
        plt.ylabel('Winds')
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        
        ax3 = Figure.add_axes([0.2, 0.4, 0.75, 0.25])
        plt.plot(self.t,self.Precip) 
        plt.plot(self.t,self.Press/1000.0) 
        plt.ylabel('Pressure')
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])

        Figure.savefig('%sTemperature_and_sonars.png' % OutputFolder)
        plt.close(Figure)

      

    def make_time_labels(self,StartDate,EndDate,time_ax = None):
        """ Function to plot the different model run outputs.
        """
        MM = ['Jan. 01','Feb. 01',
        'Mar. 01','Apr. 01','May. 01','Jun. 01','Jul. 01',
        'Aug. 01','Sep. 01','Oct. 01','Nov. 01','Dec. 01',]
        reference_date = datetime(2022,1,1,hour=0)
        start_date = datetime(StartDate[0],StartDate[1],StartDate[2])
        lim1 = 44562.0 +(start_date  - reference_date).days + ((start_date  - reference_date).seconds /(24*60*60))
        end_date = datetime(EndDate[0],EndDate[1],EndDate[2])
        lim2 = 44562.0 +(end_date  - reference_date).days + ((end_date  - reference_date).seconds /(24*60*60))
        nmonths = int(ceil((end_date - start_date).days)/30)

        DateLabels = []
        DateTicks = []
        Month_tick = int(StartDate[1])-1
		
        for months in range(0, nmonths+2):
            k = int((Month_tick + months) - 12*(floor((Month_tick + months)/12)))
            DateLabels.append(MM[k])
            Yeark = StartDate[0] + floor((Month_tick + months)/12)
			
            time_1 = datetime(Yeark,k+1,1)
            time_2 = datetime(Yeark,k+1,16)
            addedtime = 44562.0 +(time_1 - reference_date).days + ((time_1 - reference_date).seconds /(24*60*60))
            DateTicks.append(addedtime)
            addedtime = 44562.0 +(time_2 - reference_date).days + ((time_2 - reference_date).seconds /(24*60*60))
            DateTicks.append(addedtime)
            DateLabels.append(' ')
        
        if time_ax is not None:
            t0 = time_ax[0]
            tlast = time_ax[len(time_ax)-1]
            DeltaTime = tlast-t0
            DeltaPixel = len(time_ax)-1
            itt = 0
            for tics in DateTicks:
                new_tic = ((tics-t0)/DeltaTime)*DeltaPixel
                DateTicks[itt] = new_tic
                itt = itt+1                
            lim1 = ((lim1-t0)/DeltaTime)*DeltaPixel
            lim2 = ((lim2-t0)/DeltaTime)*DeltaPixel
        plt.xticks(DateTicks,DateLabels)
        plt.xlim(lim1,lim2)
        
        return (lim1,lim2)
       

    def show_buoy_location(self, OutputFolder = None):

        cmap = cmo.cm.balance
        cmap.set_bad('grey')

        proj = ccrs.NorthPolarStereo(central_longitude=315)
        extent = [298.5, 298.15, 56.4,56.6]

        land_10m = cartopy.feature.GSHHSFeature(scale='h', levels=[1],facecolor=cartopy.feature.COLORS['land'])    
    
        Figure1 = plt.figure(figsize=[4.0,4.0])
        
        ax =  Figure1.add_axes([0.0, 0.0, 1.0, 1.0],projection=proj)
        
        ax.set_extent(extent)

        ax.add_feature(cartopy.feature.OCEAN)

        ax.add_feature(land_10m)

        # mark a known place to help us geo-locate ourselves
        ax.plot(self.Lon[:], self.Lat[:], '*r', markersize=5, transform=ccrs.Geodetic())
        

        ax.plot(-61.696888, 56.541683, 'ok', markersize=5, transform=ccrs.Geodetic())
        ax.text(-61.752, 56.525, 'Nain',transform=ccrs.Geodetic())

        Figure1.savefig("%sBuoys_location.png" % OutputFolder)





