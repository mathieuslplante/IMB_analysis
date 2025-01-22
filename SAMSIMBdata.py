
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

class SAMSIMB_DAsetup:
    def __init__(self,config = None):
        self.DataFile = config['DataInPath']
        self.OutputFolder = config['OutPath']
        self.Type = config['Type']
        self.Plateau = config['Plateau']
        self.nsensors = int(config['nsensors']) 
        self.tstep = int(config['dt'])
        self.Tair_Header = config['Tair_Header'] 
        self.Tocean_Header = config['Tocean_Header'] 
        self.Sensor0_Header = config['Sensor0_Header']
        self.SensorLast_Header = config['SensorLast_Header'] 
        self.Time_Header = config['Time_Header']
        self.Lon_Header = config['Lon_Header'] 
        self.Lat_Header = config['Lat_Header'] 


            

class SAMSIMBdata:
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
        
        indt149 = self.find_unstruct(header,'T149')
        indt90 = self.find_unstruct(header,'T90')
        
        date_all = data[ExpSetup.Time_Header][:]
        
        reference_date = datetime(2022,1,1,hour=0)

        #--------------------------
        # Get the number of time steps and start row (kinit)
        #--------------------------

        #If there is no time specifications, we read the entire data .csv files
        #and find total number of data
        if time == None:
            nrows = len(date_all)
            kinit = 0
            
            #If we specified a start time, we need to
            # find the corresponding first row:
            if StartTime != None:
                
                date_k = data[ExpSetup.Time_Header][kinit]
                ThisDate = self.extract_time(date_data = date_k)
                
                #loop until we find the row corresponding to the start time
                while ThisDate < StartTime:
                    kinit = kinit+1
                    date_k = data[ExpSetup.Time_Header][kinit]
                    ThisDate = self.extract_time(date_data = date_k)
                    
            date_init = date_all[kinit]
            date_init = self.extract_time(date_data = date_init)
            
            date_last = date_all[nrows-1]
            date_last = self.extract_time(date_data = date_last)
            
            nstep = int(int((date_last - date_init).days)*24/ExpSetup.tstep + int((date_last - date_init).seconds/ExpSetup.tstep*60*60 +1))
            
            
        #If time specifications are included, we fetch the specific time interval 
        else: 
            nstep = time.nstep
            kinit = 0
            date_k = data[ExpSetup.Time_Header][kinit]
            #finding the date of the first data, based on the 
            #reference (Jan 01 2022, 0Z)
            ThisDate = self.extract_time(date_data = date_k)
            #Get the row corresponding to the start time
            while ThisDate < time.StartDate:
                kinit = kinit+1
                date_k = data[ExpSetup.Time_Header][kinit]
                ThisDate = self.extract_time(date_data = date_k)
                
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
            date_k = self.extract_time(date_data = date_k)
            
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
                
                #Filter out invalid data:
                if ((data['T10'][kdata] > -60) & 
                            (data['T10'][kdata] < 30) & 
                            (data['T60'][kdata] < -1.0) &
                            (data['MsgType'][kdata]==0)):  
                    #Get the first part of the thermistor string data
                    if data['TotalMsgPart'][kdata] == 1:
                        print(data[header[indtair]][kdata])
                        self.airT[kout]=data[header[indtair]][kdata]
                        print(self.airT[kout])
                        for indname in range(indt0, indt149+1):
                            self.data[kout,indname-indt0+1]= data[header[indname]][kdata]
                            #Skip a row in the PI buoy case
                            if (loc=='PI1'):
                                if (indname >= indt0+46):
                                    self.data[kout,indname-indt0]= data[header[indname]][kdata]
                    
                    #Second part of the thermistor string data                
                    if data['TotalMsgPart'][kdata] == 2:
                        for indname in range(indt0, indt90+1):
                            self.data[kout,indname-indt0+151] = data[header[indname]][kdata]
          
                        self.datelist[kout] =  new_date.strftime("%Y%m%d%H")
                        self.Lat[kout]= data[ExpSetup.Lat_Header][kdata]
                        self.Lon[kout]= data[ExpSetup.Lon_Header][kdata]
                        self.t[kout]= new_time
                        self.oceT[kout]=data[header[indtwater]][kdata]

                    
                        if (np.amin(self.data[kout,2:ExpSetup.nsensors+2]) < -60.0):
                            self.data[kout,:] = np.nan
                        if (np.amax(self.data[kout,2:ExpSetup.nsensors+2]) > 30):
                            self.data[kout,:] = np.nan
                            
                        print(self.t[kout],self.airT[kout],self.oceT[kout])
                    
            
    def find_unstruct(self,header, name1):
        b = len(header)
        for kname in range(0,b):
            if header[kname]==name1:
                indname = kname
        return indname
    
    def extract_time(self,date_data = None):
        year_data = int(date_data[0:4])
        month_data = int(date_data[5:7])
        day_data = int(date_data[8:10])
        hour_data = int(date_data[11:13])
        time_k = time(hour_data,0,0)
        date_k = date(year_data,month_data,day_data)
        date_time_k = datetime.combine(date_k,time_k)
        
        return date_time_k

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







   
class SAMSIMBheatdata:
    def __init__(self, 
                 ExpSetup=None, time = None,StartTime = None):

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
        date_all = data[ExpSetup.Time_Header][:]
        reference_date = datetime(2022,1,1,hour=0)

        #--------------------------
        # Get the start and end rows for the analysis
        #--------------------------

        #If we read the entire data csv
        if time == None:
            nrows = len(date_all)
            kinit = 0
            
            #If we specified a start time, we need to
            # find the corresponding first row:
            if StartTime != None:
                
                date_k = data[ExpSetup.Time_Header][kinit]
                ThisDate = self.extract_time(date_data = date_k)
                
                #loop until we find the row corresponding to the start time
                while ThisDate < StartTime:
                    kinit = kinit+1
                    date_k = data[ExpSetup.Time_Header][kinit]
                    ThisDate = self.extract_time(date_data = date_k)
                    
            date_init = date_all[kinit]
            date_init = self.extract_time(date_data = date_init)
            
            date_last = date_all[nrows-1]
            date_last = self.extract_time(date_data = date_last)
            
            nstep = int((date_last - date_init).days)
            
            print(date_init)
            print(date_last)
            
        #If we want a specific time interval 
        #           (requires time object)
        else: 
            nstep = time.ndays
            kinit = 0
            date_k = data[ExpSetup.Time_Header][kinit]
            #finding the date of the first data, based on the 
            #reference (Jan 01 F2022, 0Z)
            ThisDate = self.extract_time(date_data = date_k)
            #Get the row corresponding to the start time
            while ThisDate < time.StartDate:
                kinit = kinit+1
                date_k = data[ExpSetup.Time_Header][kinit]
                ThisDate = self.extract_time(date_data = date_k)
                
        print(kinit, nstep)
        
        #--------------------------
        # Prepare datasets and initializing parameters
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
        self.data = np.zeros((nstep,ExpSetup.nsensors+1,2))
        self.data[:,:,:] = np.nan
        
        #initializing memory variables of previous date and hour to impossible values
        prev_date = reference_date - timedelta(days=10000)
        prev_hour = -1
        
        smooth_t = 1    #width of smooting in time
        smooth_z = 2    #range of smooting in z
        
        #Make vertical vector corresponding to the sensor positions
        z = np.zeros((ExpSetup.nsensors+1,1))
        for k in range(1, ExpSetup.nsensors+1):
            z[k] = k*-2.0
        zf = len(z)

        ah,bh,ch = self.data.shape
        heat_ratio = np.zeros((ah,bh))*np.nan   #data array for the heating cycles
        first_der = np.zeros((ah,bh))*np.nan    #data array of the 1st order derivative in Z
        second_der_h = np.zeros((ah,bh))*np.nan #data array of the 2nd order derivative in Z
    


        #--------------------------
        # Fetch hourly data
        #   Note: here, we re-organise the data into 1hour
        #   bins, so that we can make a pcolor. This means 
        #   that we are skipping few duplicate transmissions
        #-------------------------- 
        
        for kdata in range(kinit,a):
            
            #Update the data
            date_k = data[ExpSetup.Time_Header][kdata]
            date_k = self.extract_time(date_data = date_k)
            
            #record the first date
            if kdata == kinit:
                ref_time = date_k
            
            #Calculate the row indice
            new_date = date_k
            new_hour = (date_k - ref_time).days
            kout = new_hour
            print(new_date,new_hour,kout)
            
            new_time = 44562.0 + ((date_k - reference_date).seconds /(24*60*60))

            #----------------------------------------------
            # Extract data
            #----------------------------------------------
            if kout < nstep:
                
                #Get temperature data
                if data['MsgType'][kdata]==1:
                    self.data[kout,0,0] = kdata
                    for indname in range(indt0, indtlast+1):
                        self.data[kout,indname-indt0+1,0] = data[header[indname]][kdata]
      
                if data['MsgType'][kdata]==2:
                    self.data[kout,0,1] = kdata
                    for indname in range(indt0, indtlast+1):
                        self.data[kout,indname-indt0+1,1] = data[header[indname]][kdata]
          
                        self.datelist[kout] =  new_date.strftime("%Y%m%d%H")
            
                        self.Lat[kout]= data[ExpSetup.Lat_Header][kdata]
                        self.Lon[kout]= data[ExpSetup.Lon_Header][kdata]
                        self.t[kout]= new_time
                        
                        
    

        #--------------------------------------------------------------
        #  Compute the heat cycles
        #--------------------------------------------------------------

        #Get the temperature change between 1min and 2min of heating
        for k1 in range(0, ah):
            for kz in range(1, zf):
                heat_ratio[k1,kz] = (self.data[k1,kz,1]/ self.data[k1,kz,0])
        

        #--------------------------------------------------------------
        #  Smoothing
        #--------------------------------------------------------------
    
    
        #Smoothing the heat ratio in time. Smooth_t is the number of points used for the smooting.
        heat_ratio_sm = heat_ratio.copy()
        for k1 in range(0+smooth_t, ah-smooth_t):
            if (smooth_t > 0):
                for n in range(0, smooth_t):
                    if (n == 0):
                        heat_ratio_sm[k1,:]  = heat_ratio[k1,:]
                        denom = 1.0
                    else:
                        denom = denom+2
                        heat_ratio_sm[k1,:]= heat_ratio_sm[k1,:] + heat_ratio[k1-n,:] + heat_ratio[k1+n,:]
                heat_ratio_sm[k1,:] = heat_ratio_sm[k1,:]/denom
        
        #Smoothing the heat ratio vertically. Smooth_z is the number of points used for the smooting.
        for k1 in range(0, ah-1):
            for kz in range(1+smooth_z, zf-smooth_z):
                if (smooth_z > 0):
                    for n in range(0, smooth_z):
                        if (n == 0):
                            heat_ratio_sm[k1,kz] = heat_ratio_sm[k1,kz]
                            denom = 1.0
                        else:
                            denom = denom+2
                            heat_ratio_sm[k1,kz]= heat_ratio_sm[k1,kz] + heat_ratio[k1,kz-n] + heat_ratio[k1,kz+n]
                    heat_ratio_sm[k1,kz] = heat_ratio_sm[k1,kz]/denom
 
 
        #--------------------------------------------------------------
        #  Compute the derivatives
        #--------------------------------------------------------------

            # Compute the first Z-derivative (center difference)
            for kz in range(2, zf-2):
                first_der[k1,kz] = (heat_ratio_sm[k1,kz+2] - heat_ratio_sm[k1,kz-2])/8.0
        
            # Compute the second Z-derivative (center difference)
            for kz in range(2, zf-2):
                second_der_h[k1,kz] = (first_der[k1,kz+1] - first_der[k1,kz-1])/4.0
        
        #Record data in objects:
        self.heat_ratio = heat_ratio_sm
        self.second_der = second_der_h
        
        
    def find_unstruct(self,header, name1):
        b = len(header)
        for kname in range(0,b):
            if header[kname]==name1:
                indname = kname
        return indname
    

    def extract_time(self,date_data = None):
        year_data = int(date_data[0:4])
        month_data = int(date_data[5:7])
        day_data = int(date_data[8:10])
        hour_data = int(date_data[11:13])
        date_k = date(year_data,month_data,day_data)
        time_k = time(hour_data,0,0)
        date_time_k = datetime.combine(date_k,time_k)
        
        return date_time_k

