
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
from math import floor,ceil


import cartopy.crs as ccrs
import cartopy.feature

import shapely.geometry as sgeom
import cmocean as cmo

class IcepackDatasetup:
	
    def __init__(self,config = None, config_Labels = None,config_rtrvl = None):
		
		#IO
        self.OutputFolder = config['OutPath'] 
        self.DataIn = config['DataInPath']
        #Dimensions
        self.tstep = int(config['dt']) #in hour
        self.nslyr = int(config['nslyr'])
        self.nilyr = int(config['nilyr'])
        self.ncat =int(config['ncat'])
        
        #Diagnostic labels
        self.SICLabel = config_Labels['SICLabel']
        self.HiavgLabel = config_Labels['HiavgLabel']
        self.HsavgLabel= config_Labels['HsavgLabel']
        self.TairLabel = config_Labels['TairLabel']
        self.SWLabel = config_Labels['SWLabel']
        self.LWLabel = config_Labels['LWLabel']
        self.SnowPcPLabel = config_Labels['SnowPcPLabel']
        self.RainPcPLabel = config_Labels['RainPcPLabel']
        self.SavgLabel = config_Labels['SavgLabel'] 
        self.TsfcLabel = config_Labels['TsfcLabel']
        self.LWoutLabel = config_Labels['LWoutLabel']
        self.LWabsLabel = config_Labels['LWabsLabel']
        self.FsensLabel = config_Labels['FsensLabel'] 
        self.FlatentLabel = config_Labels['FlatentLabel']
        self.FsurfLabel = config_Labels['FsurfLabel']
        self.CondSublLabel = config_Labels['CondSublLabel']
        self.TopMeltLabel = config_Labels['TopMeltLabel']
        self.BotMeltLabel = config_Labels['BotMeltLabel']
        self.LatMeltLabel = config_Labels['LatMeltLabel'] 
        self.NewIceLabel = config_Labels['NewIceLabel']
        self.CongLabel = config_Labels['CongLabel'] 
        self.SnowiceLabel = config_Labels['SnowiceLabel'] 
        self.IceVolLabel = config_Labels['IceVolLabel']
        self.IntEnergyChgLabel = config_Labels['IntEnergyChgLabel']
        self.SSTLabel= config_Labels['SSTLabel']
        self.SSSLabel= config_Labels['SSSLabel']
        self.FreezeTempLabel= config_Labels['FreezeTempLabel']
        self.FocnLabel= config_Labels['FocnLabel']
        self.FcondbotLabel= config_Labels['FcondbotLabel']
        self.FbotLabel= config_Labels['FbotLabel']
        self.CumulSnowMeltLabel = config_Labels['CumulSnowMeltLabel']
        self.CumulTopMeltLabel = config_Labels['CumulTopMeltLabel']
        self.CumulBotMeltLabel = config_Labels['CumulBotMeltLabel']
        self.CumulLatMeltLabel = config_Labels['CumulLatMeltLabel']
        self.CumulCatMeltLabel = config_Labels['CumulCatMeltLabel'] 
        self.CumulNewIceLabel = config_Labels['CumulNewIceLabel']
        self.CumulCongLabel = config_Labels['CumulCongLabel'] 
        self.CumulSnowiceLabel = config_Labels['CumulSnowiceLabel']
        self.CumulCatAreaChgLabel = config_Labels['CumulCatAreaChgLabel']
        self.SiceLabel = config_Labels['SiceLabel'] 
        self.TiceLabel = config_Labels['TiceLabel']
        self.QiceLabel = config_Labels['QiceLabel']
        self.TsnowLabel = config_Labels['TsnowLabel'] 
        self.AlbIceLabel  = config_Labels['AlbIceLabel']
        self.AlbSnowLabel = config_Labels['AlbSnowLabel'] 
        self.AlbPndLabel  = config_Labels['AlbPndLabel'] 
        self.LiqFractionLabel = config_Labels['LiqFractionLabel']
        self.LyrPermLabel = config_Labels['LyrPermLabel'] 
        

class IcepackData:
    def __init__(self, 
                 meta = None, time = None):

        #-----------------------------------
        # Load the data from the input file
        #-----------------------------------
        
        self.datafile = open(meta.DataIn,'r').readlines()
        print(meta.DataIn)
        self.ndata, self.dates = self.extract_ndata(self.datafile)
        
        #Making a standard time vector
        self.reference_date = datetime(2022,1,1,hour=0)
        self.datelist = self.dates.copy()
        self.t = self.dates.copy()
        self.Temp= []
        
        
        for kt in range(0,len(self.datelist)):
            self.t[kt] = (44562.0 +  
                   (self.datelist[kt] - self.reference_date).days + 
                   ((self.datelist[kt] - self.reference_date).seconds/(24*60*60)))
        
        #--------------------------
        # Get the start and end rows for the analysis
        #--------------------------

        self.nstep = self.ndata
        self.tstep = meta.tstep      
        self.OutputFolder = meta.OutputFolder
        self.nslyr = meta.nslyr
        self.nilyr = meta.nilyr
        
        #--------------------------
        # retrieve all fields
        #-------------------------- 
                
        self.SIC=self.extract_diagnostic(field=meta.SICLabel,nlayers=1,FieldName = "area")
        self.Hiavg=self.extract_diagnostic(field=meta.HiavgLabel,nlayers=1,FieldName = "thickness")
        self.Hsavg=self.extract_diagnostic(field=meta.HsavgLabel,nlayers=1,FieldName = "snowdepth")
        self.Tair=self.extract_diagnostic(field=meta.TairLabel,nlayers=1,FieldName = "Tair")
        self.Fsurf=self.extract_diagnostic(field=meta.FsurfLabel,nlayers=1,FieldName = "Fsurf")
        self.SW=self.extract_diagnostic(field=meta.SWLabel,nlayers=1,FieldName = "fsw")
        self.LW=self.extract_diagnostic(field=meta.LWLabel,nlayers=1,FieldName = "flw")
        self.SnowPcP=self.extract_diagnostic(field=meta.SnowPcPLabel,nlayers=1,FieldName = "PcPsnow")
        self.RainPcP=self.extract_diagnostic(field=meta.RainPcPLabel,nlayers=1,FieldName = "PcPrain")
        self.Savg=self.extract_diagnostic(field=meta.SavgLabel,nlayers=1,FieldName ="Avg-salinity" )
        self.Tsfc=self.extract_diagnostic(field=meta.TsfcLabel,nlayers=1,FieldName = "Tsfc")
        self.LWout=self.extract_diagnostic(field=meta.LWoutLabel,nlayers=1,FieldName = "flwout")
        self.LWabs=self.extract_diagnostic(field=meta.LWabsLabel,nlayers=1,FieldName = "flwabs")
        self.Fsens=self.extract_diagnostic(field=meta.FsensLabel,nlayers=1,FieldName = "Fsens")
        self.Flatent=self.extract_diagnostic(field=meta.FlatentLabel,nlayers=1,FieldName ="Flatent")
        self.CondSubl=self.extract_diagnostic(field=meta.CondSublLabel,nlayers=1,FieldName ="subl_cond")
        self.TopMelt=self.extract_diagnostic(field=meta.TopMeltLabel,nlayers=1,FieldName ="top-melt")
        self.BotMelt=self.extract_diagnostic(field=meta.BotMeltLabel,nlayers=1,FieldName ="bottom-melt")
        self.LatMelt=self.extract_diagnostic(field=meta.LatMeltLabel,nlayers=1,FieldName = "lateral-melt")
        self.NewIce=self.extract_diagnostic(field=meta.NewIceLabel,nlayers=1,FieldName ="new-ice")
        self.Cong=self.extract_diagnostic(field=meta.CongLabel,nlayers=1,FieldName ="congelation") 
        self.Snowice=self.extract_diagnostic(field=meta.SnowiceLabel,nlayers=1,FieldName ="snowice")
        self.IceVol=self.extract_diagnostic(field=meta.IceVolLabel,nlayers=1,FieldName ="volume")
        self.IntEnergyChg=self.extract_diagnostic(field=meta.IntEnergyChgLabel,nlayers=1,FieldName ="enrgy-chng")
        self.SST = self.extract_diagnostic(field=meta.SSTLabel,nlayers=1,FieldName ="SST")
        self.SSS = self.extract_diagnostic(field=meta.SSSLabel,nlayers=1,FieldName ="SSS")
        self.FreezeTemp = self.extract_diagnostic(field=meta.FreezeTempLabel,nlayers=1,FieldName ="Freezing_temp")
        self.Focn = self.extract_diagnostic(field=meta.FocnLabel,nlayers=1,FieldName ="Focn")
        self.Fcondbot = self.extract_diagnostic(field=meta.FcondbotLabel,nlayers=1,FieldName ="Fcond_bot")
        self.Fbot = self.extract_diagnostic(field=meta.FbotLabel,nlayers=1,FieldName ="Fbot")
        self.CumulSnowMelt=self.extract_diagnostic(field=meta.CumulSnowMeltLabel,nlayers=1,FieldName = "cumul_snmelt")
        self.CumulTopMelt=self.extract_diagnostic(field=meta.CumulTopMeltLabel,nlayers=1,FieldName = "cumul_topmelt")
        self.CumulBotMelt=self.extract_diagnostic(field=meta.CumulBotMeltLabel,nlayers=1,FieldName = "cumul_botmelt")
        self.CumulLatMelt=self.extract_diagnostic(field=meta.CumulLatMeltLabel,nlayers=1,FieldName = "cumul_latmelt")
        self.CumulCatMelt=self.extract_diagnostic(field=meta.CumulCatMeltLabel,nlayers=1,FieldName = "cumul catmelt")
        self.CumulNewIce=self.extract_diagnostic(field=meta.CumulNewIceLabel,nlayers=1,FieldName = "cumul_newice")
        self.CumulCong=self.extract_diagnostic(field=meta.CumulCongLabel,nlayers=1,FieldName = "cumul_cong")
        self.CumulSnowice=self.extract_diagnostic(field=meta.CumulSnowiceLabel,nlayers=1,FieldName = "cumul_snowice")
        self.CumulCatAreaChg=self.extract_diagnostic(field=meta.CumulCatAreaChgLabel,nlayers=1,FieldName = "cumul_cat_area")
        self.CongTotal = self.CumulCong + self.CumulNewIce - self.CumulBotMelt
        self.Sice=self.extract_diagnostic(field=meta.SiceLabel,nlayers=meta.nilyr,FieldName = "Salin_int")
        self.Tice=self.extract_diagnostic(field=meta.TiceLabel,nlayers=meta.nilyr,FieldName = "Tice_int")
        self.Qice=self.extract_diagnostic(field=meta.QiceLabel,nlayers=meta.nilyr,FieldName = "Qice_int")        
        self.Tsnow=self.extract_diagnostic(field=meta.TsnowLabel,nlayers=meta.nslyr,FieldName = "Tsnow_int")

        self.AlbIce=self.extract_diagnostic(field=meta.AlbIceLabel,nlayers=meta.ncat,FieldName = "albice")
        self.AlbSnow=self.extract_diagnostic(field=meta.AlbSnowLabel,nlayers=meta.ncat,FieldName = "albsno" )
        self.AlbPnd=self.extract_diagnostic(field=meta.AlbPndLabel,nlayers=meta.ncat,FieldName = "albpnd")
        
        self.LiqFraction=self.extract_diagnostic(field=meta.LiqFractionLabel,nlayers=meta.nilyr,FieldName = "liquidfrac")
        self.LyrPerm=self.extract_diagnostic(field=meta.LyrPermLabel,nlayers=meta.nilyr,FieldName = "lyrperm")
        
        try:
            self.w=self.extract_diagnostic(field="cat. darcy vel. (m)",nlayers=meta.ncat,FieldName = "Darcy_cat")
        except:
            self.w=self.CatTsfc.copy()*np.nan
        try:
            self.dSdt=self.extract_diagnostic(field="salinity changes",nlayers=meta.nilyr,FieldName = "dSdt")
        except:
            self.dSdt=self.LiqFraction.copy()*np.nan 
        
        self.data = []
        
        #Get the ice mass balance details to print
        DateMax = self.datelist[np.squeeze(self.Hiavg[:]==np.nanmax(self.Hiavg))]
        DateMax = DateMax[0] 
        Congel = self.CumulCong[np.squeeze(self.Hiavg[:]==np.nanmax(self.Hiavg))]
        Congel = Congel[0] 
        NewIce = self.CumulNewIce[np.squeeze(self.Hiavg[:]==np.nanmax(self.Hiavg))]
        NewIce = NewIce[0] 
        BotMelt = self.CumulBotMelt[np.squeeze(self.Hiavg[:]==np.nanmax(self.Hiavg))]
        BotMelt = BotMelt[0] 
        CongTotal= self.CongTotal[np.squeeze(self.Hiavg[:]==np.nanmax(self.Hiavg))]
        CongTotal = CongTotal[0]
        Snice = self.CumulSnowice[np.squeeze(self.Hiavg[:]==np.nanmax(self.Hiavg))]
        Snice = Snice[0] 
        
        #Print the ice mass balance in a text file
        with open('%sPrinted_MassBalance.txt' % self.OutputFolder, 'w') as f:	
				
            f.write('Initial ice thickness and snow depth with date: %s, %s, %s \n' % (self.Hiavg[0],self.Hsavg[0],self.datelist[0])) 
            f.write('Max ice thickness and date: %s, %s \n' % (np.nanmax(self.Hiavg),DateMax)) 
            f.write('Total ice congelation and snow-ice: %s, %s \n' % (CongTotal, Snice)) 
            f.write('Congelation split in growht, melt, newice: %s, %s, %s \n' % (Congel, BotMelt, NewIce)) 
            
            DateMax = self.datelist[np.squeeze(self.Hsavg[:]==np.nanmax(self.Hsavg))]
            DateMax = DateMax[0] 
            f.write('Max snow thickness and date: %s, %s \n' % (np.nanmax(self.Hsavg),DateMax)) 
            hfb = self.Hiavg - (self.Hsavg*330 + self.Hiavg*917)/1026
            for k in range(0,len(hfb)):
                f.write('hi, hs, freeboard is : %s, %s, %s,,,,,%s, %s,              on %s \n' % (self.Hiavg[k],self.Hsavg[k],hfb[k], self.CongTotal[k],self.CumulSnowice[k],self.datelist[k]))


        return

    def make_time_labels(self,StartDate,EndDate):
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
            
        plt.xticks(DateTicks,DateLabels)
        plt.xlim(lim1,lim2)
        
        return (lim1,lim2)
        
            
    def extract_ndata(self,File = None):
		#This code is used to find the number of outputs (timesteps) 
		#and the corresponding dates in the Icepack diagnostic file
        kt = 0
        dates = []
        
        for line in File:
            if "istep1:" in line:
                columns = line.split(' ')
                for i in range(len(columns)-1, 0, -1):
                    if columns[i]=='':
                        del columns[i]
                
                yeark = columns[3][0:4]
                monthk = columns[3][4:6]
                dayk = columns[3][6:8]
                hourk = int(int(columns[5])/3600)
                minutek = int((int(columns[5]) - hourk*3600)/60)
                secondk = int(int(columns[5]) - hourk*3600 - minutek*60)
                time_k = time(hourk,minutek,secondk)
                date_k = date(int(yeark),int(monthk),int(dayk))
                date_time_k = datetime.combine(date_k,time_k)
                dates.append(date_time_k)
                kt = kt+1
        Ndata = kt
        dates = np.transpose(dates)
        
        return(Ndata,dates)
      
    def extract_diagnostic(self,field=None,FieldName=None,nlayers=None):
		#This code is used to find retrieve the data from the icepack diagnostic file
		#and print a simple plot of the data.

                
        diag = np.zeros((self.ndata,nlayers))
        diag[:,:] = np.nan
        
        kt = -1
        for line in self.datafile:
            if "istep1:" in line:
                kt = kt+1
                klayer = 0
                
            elif field in line:
                columns2 = line.split(' ')
                for i in range(len(columns2)-1, 0, -1):
                    if columns2[i]=='':
                        del columns2[i]
                data = columns2[len(columns2)-1]
                if data == '*':
                    diag[kt,klayer] = np.nan
                else:
                    diag[kt,klayer] = data
                klayer=klayer+1


        #plot the data timeseries        
        if FieldName != None:                
            fieldfig = plt.figure(figsize=[6.0, 4.0]) 
            ax3 = fieldfig.add_axes([0.2, 0.2, 0.6, 0.6]) 
            figname = self.OutputFolder + FieldName +'.png'
            for kplot in range(0,nlayers):
                plt.plot(self.t,diag[:,kplot],label='layer %s' % str(kplot))
                plt.xlabel("date")
                plt.ylabel(FieldName)
                plt.grid(True)
                (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,4,16])
            plt.legend(loc='lower right',fontsize=8.0 )
            fieldfig.savefig(figname)
            plt.close(fieldfig) 
   
        return(diag)


    def Interpolate_into_IMB_profile(self,nsensors=None,meta = None):
    #------------------------------------------------
    # Interpolate the temperature into 2 cm layers, to match the obs.
    # - The initial snow-ice interface is set at z=0, 50th sensor
    #------------------------------------------------

        self.data = np.zeros((int(self.ndata),nsensors+2))*np.nan
  
        zf = np.zeros((nsensors,1))
        zf.fill(np.nan)  
        for kz in range(0, nsensors):
            zf[kz] = 100.0-(kz*2.0)

        for k in range(0,int(self.ndata)):
			
	        #-------------------------------------------------     
            #Make a vector with height of each layer and interfaces
            #-------------------------------------------------
            z = np.zeros((meta.nilyr+meta.nslyr+2,1)) *np.nan #position of the layer interfaces

            #Get the thickness of each layer
            hik = self.Hiavg[k]*100
            hsk = self.Hsavg[k]*100
            hilyrk = hik/meta.nilyr
            hslyrk = hsk/meta.nslyr
            
            #The air-snow interface
            z[0] = hsk
            #Snow layers
            for ks in range(1,meta.nslyr+1):
                z[ks] = (hslyrk)*(ks-0.5)
            #Ice layers
            for ki in range(1,meta.nilyr+1):
                z[ki+meta.nslyr] = -(hilyrk)*(ki-0.5)    
            z[meta.nslyr+meta.nilyr+1] = -hik
      
            #Raising everything when there is snow ice, as the top of the ice is above 0.
            z = z + self.CumulSnowice[k]*100.0
    
	        #-------------------------------------------------     
            # Interpolate the temperature in simulated IMB sensor positions
            #    with z=0 being the original ice surface and 50 sensors about the ice (1m) 
            #-------------------------------------------------    
            
            
            #Temperature at snow-ice interface, estimated from the temperature gradient in the ice
            Z_0 = self.CumulSnowice[k]*100.0
            T_0 = self.Tice[k,0] + (self.CumulSnowice[k]*100.0-z[1+meta.nslyr])*(self.Tice[k,1]-self.Tice[k,0])/((z[2+meta.nslyr]-z[1+meta.nslyr]))

            for kz in range(0, nsensors):
				
				#Air temperature at Tsfc above the snow
                if (zf[kz]>= z[0]):
                    self.data[k,kz] = self.Tsfc[k]
                    
                #Just below the air-snow interface 
                elif zf[kz] >= z[1]:
                    self.data[k,kz] = self.Tsfc[k] + (zf[kz]-z[0])*(self.Tsnow[k,0]-self.Tsfc[k])/((z[1]-z[0]))

                #In the snow layers  (not happening if only 1 snow layer)
                elif zf[kz] >= z[meta.nslyr]:
                    for indz in range(2,meta.nslyr+1):
                        ksnow = indz-1
                        if (zf[kz] < z[indz-1]) and (zf[kz] >= z[indz]):
                            self.data[k,kz] = self.Tsnow[k,ksnow-1] + (zf[kz]-z[indz-1])*(self.Tsnow[k,ksnow]-self.Tsnow[k,ksnow-1])/((z[indz]-z[indz-1]))
		
                #Just above the snow-ice interface
                elif zf[kz] >= Z_0:
                    self.data[k,kz] = self.Tsnow[k,meta.nslyr-1] + (zf[kz]-z[meta.nslyr])*(T_0-self.Tsnow[k,meta.nslyr-1])/((Z_0-z[meta.nslyr]))

                #Just below the snow-ice interface 
                elif zf[kz] >= z[meta.nslyr+1]:	
                    self.data[k,kz] = T_0 + (zf[kz]-Z_0)*(self.Tice[k,0]-T_0)/((z[1+meta.nslyr]-Z_0))
                    
                #In the ice layers  (not happening if only 1 ice layer)
                elif zf[kz] >= z[meta.nslyr+meta.nilyr]:	
                    for indz in range(2+meta.nslyr,meta.nslyr+meta.nilyr+1):
                        kice= indz-1-meta.nslyr
                        if (zf[kz] < z[indz-1]) and (zf[kz] >= z[indz]):
                            self.data[k,kz] = self.Tice[k,kice-1] + (zf[kz]-z[indz-1])*(self.Tice[k,kice]-self.Tice[k,kice-1])/((z[indz]-z[indz-1]))
		
				#Just above the ice-ocean interface
                elif (zf[kz] >= z[meta.nslyr+meta.nilyr+1]):
                    self.data[k,kz] = self.Tice[k,meta.nilyr-1] + (zf[kz]-z[meta.nslyr+meta.nilyr])*(self.SST[k]-self.Tice[k,meta.nilyr-1])/((z[meta.nslyr+meta.nilyr+1]-z[meta.nslyr+meta.nilyr]))

				#In the ocean
                else:  
                    self.data[k,kz] = self.SST[k]
							



