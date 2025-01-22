

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import array

from datetime import datetime, date, time, timedelta
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean as cmo
from scipy import stats
from scipy import signal
from math import floor,ceil
from matplotlib.ticker import FormatStrFormatter,MultipleLocator, Locator

# ===========================================================================================================
class visuals:
    def __init__(self):

        self.Fig3Key = 7
        Fig3 = plt.figure(self.Fig3Key,figsize=[6.0, 8.0])
        self.axT0 = np.nan
        self.Fig4Key = 10
        Fig4 = plt.figure(self.Fig4Key,figsize=[10.0, 7.0])

        self.Fig5Key = 15
        Fig5 = plt.figure(self.Fig5Key,figsize=[5.0, 6.0])
        self.axH1 = np.nan
        self.axH2 = np.nan

        self.Fig6Key = 21
        Fig6 = plt.figure(self.Fig6Key,figsize=[10.0, 5.0])

        self.Fig7Key = 25
        Fig7 = plt.figure(self.Fig7Key,figsize=[6.5, 5.0])

        self.Fig8Key = 26
        Fig8 = plt.figure(self.Fig8Key,figsize=[6.5, 5.0])

        self.F11Key = 37
        Fig11 = plt.figure(self.F11Key,figsize=[10.0, 8.0])

        self.Fig12Key = 27
        Fig12 = plt.figure(self.Fig12Key,figsize=[6.5, 5.0])
        
        self.FigSupKey = 47
        FigSup = plt.figure(self.FigSupKey,figsize=[10.0, 7.0])

    def make_y_axis(self, Zice0 = None):
        if Zice0 is not None:
            print(Zice0)
            plt.yticks([Zice0-125, Zice0-100, Zice0-75, Zice0-50,Zice0-25,
                        Zice0,Zice0+25,Zice0+50,Zice0+75],
                        ('-250', '-200', '-150', '-100','-50', '0','50','100','150'))
        else:
            plt.yticks([86, 136, 186, 236],('-300', '-200', '-100', '0'))
        plt.ylim(126, 236)

    def make_time_labels(self,StartDate,EndDate,time_ax = None, monthly = None):
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

        if monthly is None:
            print(lim1,lim2)
            DateTicks = np.arange(lim1,lim2,7)
            flag = 0
            for k in range(0,len(DateTicks)):
                delta = timedelta(weeks=k)
                newdate = start_date + delta
                print(newdate)
                thisstring = '%s%s' % (MM[newdate.month-1][0:5],newdate.day)
                print(newdate.day)
                if flag == 0:
                    DateLabels.append(' ')
                    flag = 1
                else:
                    DateLabels.append(thisstring)
                    flag = 0
                #thisstring = MM[]

        else:

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


#================================================================================================
# FIGURE 3 in Plante et al. 2023
#================================================================================================

    def figure_temperatures_weather(self,Data = None,
                                Data_model = None,
                                Data_station = None,
                                Rtrvl = None,
                                krun = None,
                                OutputFolder = None):


        nameF3 = '%sFigure3_Observed_T_and_weather.png' % OutputFolder
        Fig3 = plt.figure(self.Fig3Key)

        # First panel with recorded air temperature from IMB buoys and GDPS
        #----------------------------------------------------------------

        if krun == 0:
            self.axT0 = Fig3.add_axes([0.15, 0.83, 0.7, 0.125])
            if Data_model is not None:
                plt.plot(Data_model.t+5/24,Data_model.Tair,linestyle='-',color='k', linewidth=1, label = 'GDPS')
            plt.plot(Data.t,Data.airT,linestyle='-',color='b', label = 'SIMBA1')
            plt.plot(Data.t,Data.airT-1000,linestyle='-',color='g',label = 'SIMBA2')
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                ncol=3, borderaxespad=0.,fontsize=7.0 )

        elif krun == 1:
            plt.sca(self.axT0)
            plt.grid()
            plt.plot(Data.t,Data.airT,linestyle='-',color='g',label = 'IMB2')


        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        plt.ylim(-35, 5)
        plt.ylabel(r'$T_{air}$ ($^\circ$C)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.text(-0.15, 1.2,'a)',ha='center', va='center', transform=self.axT0.transAxes,fontsize=12)
        #plt.grid()

        # Second panel with recorded precipitations from Nain airport and GDPS
        #----------------------------------------------------------------

        if krun == 0:
            axP = Fig3.add_axes([0.15, 0.62, 0.7, 0.125] )
            Precip_station = Data_station.Precip[0:-2]+Data_station.Precip[1:-1]+Data_station.Precip[2:]
            plt.plot(Data_station.t[1:-1]-5.0/24.0,Precip_station/3.0,linestyle='-',color='black', label='Pcp. W. station')
            plt.plot(Data_model.t,Data_model.SnowPcP*1000.0/3.0,linestyle='-', color='cornflowerblue',label='snow, GDPS')
            plt.plot(Data_model.t,Data_model.RainPcP*1000.0/3.0,linestyle='-', color='blue',label='rain, GDPS')
            (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
            plt.grid()
            plt.ylabel('Precipitation (mm h$^{-1}$)',fontsize=8)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.text(-0.15, 1.2,'b)',ha='center', va='center', transform=axP.transAxes,fontsize=12)
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                         ncol=2, borderaxespad=0.,fontsize=7.0 )

        # Third and Fourth panels with recorded vertical temperature profiles from IMB1 and IMB2
        #----------------------------------------------------------------

        if krun == 0:
            ax1 = Fig3.add_axes([0.15, 0.32, 0.77, 0.2],)

        elif krun == 1:
            ax1 = Fig3.add_axes([0.15, 0.05, 0.77, 0.2])

        z0 = Rtrvl.nsensors
        xtime = np.arange(0,Rtrvl.Ndata)
        snowbot_init = np.zeros(len(Rtrvl.bottomsnow_minim))
        snowbot_init[:] = int(Rtrvl.snowbot_mode)
        icetop_init = np.zeros(len(Rtrvl.bottomsnow_minim))
        icetop_init[:] = int(Rtrvl.icetop_mode)

        ax1.set_facecolor("gray")
        data_p= np.ma.array(Data.data, mask=Data.data==np.nan)
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')
        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 = -snowbot_init[0]+z0)
        #plt.clim(-20.0,-1.8) # Limits of colorbar
        plt.ylabel('sensor position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)

        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data.t)
        plt.plot(xtime,-Rtrvl.topsnow_minim+z0,linestyle='-',color='b', label='snow surface')
        plt.plot(xtime,-snowbot_init+z0,linestyle='--',color='r', label='plateau')
        plt.plot(xtime,-icetop_init+z0,[-Rtrvl.snowbot_mode+z0,-Rtrvl.snowbot_mode+z0],linestyle='--',color='r',label='')
        plt.plot(xtime,-Rtrvl.icetop_minim+z0,linestyle='-',color='k', label='ice interface')
        plt.plot(xtime,-Rtrvl.icebottom_minim+z0,linestyle='-',color='k', label='')

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )

        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])
        clb.ax.tick_params(labelsize=8)
        plt.text(0.5,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax1.transAxes,fontsize=8)

        if krun == 0:
            plt.text(-0.15, 1.2,'c)',ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'d)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
            ax1.arrow(0.37, 0.33, 0.06, 0.25,transform=ax1.transAxes,color='rebeccapurple',linewidth = 1,head_width=0.01,overhang=0,zorder=100)
        Fig3.savefig(nameF3, dpi=600)

        if krun == 1:
            plt.close(Fig3)

#================================================================================================
# FIGURE 4 in Plante et al. 2023
#================================================================================================

    def figure_RtrvlValidation(self,Data=None,
                                time=None,
                                Rtrvl = None,
                                Data_heat = None,
                                krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                OutputFolder = None):

        nameF4 = '%sFigure4_validatione_dtdzHeat.png' % OutputFolder

        z0 = Rtrvl.nsensors
        xtime = np.arange(0,time.nstep)
        a,b = Data.data.shape
        snowbot_init = np.zeros(len(Rtrvl.bottomsnow_minim))
        snowbot_init[:] = int(Rtrvl.snowbot_mode)
        icetop_init = np.zeros(len(Rtrvl.bottomsnow_minim))
        icetop_init[:] = int(Rtrvl.icetop_mode)

        # Top pannels showing the cooling / warming rates at each sensor
        #----------------------------------------------------------------

        data_f1 = Data.data[int(24/Data.tstep):a,:] #Temperature at time=n, in Celsius
        data_f2 = Data.data[0:a-int(24/Data.tstep),:] #Temperature at time=n+1, in Celsius
        data_dt = (data_f1-data_f2)/24.0 # Rate of change in temperature, in C/day

        Fig4 = plt.figure(self.Fig4Key)
        ax0 = Fig4.add_axes([0.07+krun*0.5, 0.68, 0.38, 0.2])
        ax0.set_facecolor("gray")

        im = plt.imshow(np.flipud(np.transpose(data_dt)),cmap=cmo.cm.balance,aspect='auto')
        self.make_y_axis(Zice0 = -snowbot_init[0]+z0)
        plt.clim(-0.05,0.05) # Limits of colorbar
        plt.ylabel('sensor position (cm)', fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data.t)
        plt.plot(xtime,-Rtrvl.topsnow_minim+z0,linestyle='-',color='b', label='snow surface')
        plt.plot(xtime,-snowbot_init+z0,linestyle='--',color='r', label='plateau')
        plt.plot(xtime,-icetop_init+z0,[-Rtrvl.snowbot_mode+z0,-Rtrvl.snowbot_mode+z0],linestyle='--',color='r',label='')
        plt.plot(xtime,-Rtrvl.icetop_minim+z0,linestyle='-',color='k', label='ice interface')
        plt.plot(xtime,-Rtrvl.icebottom_minim+z0,linestyle='-',color='k', label='')

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, mode="expand", borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax0)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1)#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(-0.5,0.06,r'$\frac{\partial T}{\partial t}$ ($^\circ$C h$^{-1}$)',fontsize=8)

        if krun ==0:
            plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=ax0.transAxes,fontsize=12)
        elif krun== 1:
            plt.text(-0.15, 1.2,'b)', ha='center', va='center', transform=ax0.transAxes,fontsize=12)
            ax0.arrow(0.355, 0.33, 0.06, 0.25,transform=ax0.transAxes,color='rebeccapurple',linewidth = 1,head_width=0.01,overhang=0,zorder=100)

        # Middle pannels showing the vertical temperature gradient at each sensor
        #----------------------------------------------------------------

        first_der= np.ma.array(Rtrvl.first_der, mask=Rtrvl.first_der==np.nan)

        cmap = cmo.cm.balance
        cmap.set_bad('gray')

        ax0z = Fig4.add_axes([0.07+krun*0.5, 0.38, 0.38, 0.2])
        ax0z.set_facecolor("gray")
        im = plt.imshow(np.flipud(np.transpose(-first_der*100)),cmap=cmap,aspect='auto')
        self.make_y_axis(Zice0 = -snowbot_init[0]+z0)
        plt.clim(-50.0,50.0) # Limits of colorbar
        plt.ylabel('sensor position (cm)', fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data.t)

        plt.plot(xtime,-Rtrvl.topsnow_minim+z0,linestyle='-',color='b', label='snow surface')
        plt.plot(xtime,-snowbot_init+z0,linestyle='--',color='r', label='plateau')
        plt.plot(xtime,-icetop_init+z0,[-Rtrvl.snowbot_mode+z0,-Rtrvl.snowbot_mode+z0],linestyle='--',color='r',label='')
        plt.plot(xtime,-Rtrvl.icetop_minim+z0,linestyle='-',color='k', label='ice interface')
        plt.plot(xtime,-Rtrvl.icebottom_minim+z0,linestyle='-',color='k', label='')

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, mode="expand", borderaxespad=0.,fontsize=7.0 )

        divider = make_axes_locatable(ax0z)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1)#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(-0.5,60,r'$\frac{\partial T}{\partial z}$ ($^\circ$C $m^{-1}$)',fontsize=8)

        if krun == 0:
            plt.text(-0.15, 1.2,'c)', ha='center', va='center', transform=ax0z.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'d)', ha='center', va='center', transform=ax0z.transAxes,fontsize=12)

        # Bottom pannels showing the warming during heating cycles at each sensor
        #----------------------------------------------------------------

        heat_post1 = np.zeros((a,b))*np.nan
        heat_post2 = np.zeros((a,b))*np.nan
        ratio_post = np.zeros((a,b))*np.nan

        for k in range(0,a):
            heat_post1[k,:] =Data_heat.data[int(k/int(24/time.tstep)),:,0]
            heat_post2[k,:] =Data_heat.data[int(k/int(24/time.tstep)),:,1]
            ratio_post[k,:] =Data_heat.heat_ratio[int(k/int(24/Data.tstep)),:]



        ax1 = Fig4.add_axes([0.07+krun*0.5, 0.08, 0.38, 0.2])
        ax1.set_facecolor("gray")

        ratio_post= np.ma.array(ratio_post, mask=ratio_post==np.nan)
        cmap = cmo.cm.thermal
        cmap.set_bad('gray')
        im2 =plt.imshow(np.flipud(np.transpose(ratio_post)),cmap=cmap, aspect = 'auto')
        plt.clim(0.5,3) # Limits of colorbar

        plt.plot(xtime,-Rtrvl.topsnow_minim+z0,linestyle='-',color='b', label='snow surface')
        plt.plot(xtime,-snowbot_init+z0,linestyle='--',color='r', label='plateau')
        plt.plot(xtime,-icetop_init+z0,[-Rtrvl.snowbot_mode+z0,-Rtrvl.snowbot_mode+z0],linestyle='--',color='r',label='')
        plt.plot(xtime,-Rtrvl.icetop_minim+z0,linestyle='-',color='k', label='ice interface')
        plt.plot(xtime,-Rtrvl.icebottom_minim+z0,linestyle='-',color='k', label='')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, mode="expand", borderaxespad=0.,fontsize=7.0 )

        self.make_y_axis(Zice0 = -snowbot_init[0]+z0)

        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data.t)
        cmap.set_bad('gray')

        plt.ylabel('sensor position (cm)', fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)

        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im2,cax=cax1)
        clb.ax.tick_params(labelsize=8)
        plt.text(-0.5,3.25,r'$\Delta$T ($^\circ$C)',fontsize=8)

        if krun ==0:
            plt.text(-0.15, 1.2,'e)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'f)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)

        Fig4.savefig(nameF4, dpi=600)

        if krun == 1:
            plt.close(Fig4)

#================================================================================================
# FIGURE 5 in Plante et al. 2023
#================================================================================================

    def figure_thickness_obs(self,
                            data_buoy = None,
                            rtrvl = None,
                            krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                            OutputFolder=None):

        max_buoy = 1.0

        hfb_buoy_long = rtrvl.hi_int - (rtrvl.hs_int*330 + rtrvl.hi_int*917)/1026
        plt.rc('legend',fontsize=12)

        Fig5 = plt.figure(self.Fig5Key)

        # First pannel showing the ice thickness, snow thickness, and freeboard
        #----------------------------------------------------------------

        if krun == 0:
            self.axH1 = Fig5.add_axes([0.15, 0.575, 0.75, 0.3])
            plt.sca(self.axH1)
            plt.plot(data_buoy.t,rtrvl.hi_int,linestyle='-', color='b',linewidth=2)
            plt.plot(data_buoy.t,rtrvl.hs_int,linestyle='-', color='darkgreen',linewidth=2)
            plt.plot(data_buoy.t,hfb_buoy_long,linestyle='-', color='darkorange',linewidth=2)
        elif krun == 1:
            plt.sca(self.axH1)
            plt.plot(data_buoy.t,rtrvl.hi_int,linestyle='-', color='c',linewidth=2)
            plt.plot(data_buoy.t,rtrvl.hs_int,linestyle='-', color='yellowgreen',linewidth=2)
            plt.plot(data_buoy.t,hfb_buoy_long,linestyle='-', color='gold',linewidth=2)

        max_buoy = np.nanmax([max_buoy, np.nanmax(rtrvl.hi_int)])
        plt.grid()

        maxlim = max_buoy*1.05
        plt.ylim(-10.0, maxlim)
        plt.ylabel('thickness (cm)')
        plt.grid(True)
        if krun == 1:
            plt.legend(('$h_\mathrm{i}$ SIMBA1','$h_\mathrm{s}$ SIMBA1','$h_\mathrm{fb}$ SIMBA1',
                        '$h_\mathrm{i}$ SIMBA2','$h_\mathrm{s}$ SIMBA2','$h_\mathrm{fb}$ SIMBA2'),bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=2, mode="expand", borderaxespad=0.,fontsize=10.0 )
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=self.axH1.transAxes,fontsize=14)



        # Second pannel showing the congelation ice and snow-ice formation
        #----------------------------------------------------------------

        if krun == 0:

            self.axH2 = Fig5.add_axes([0.15, 0.075, 0.75, 0.3])
            plt.sca(self.axH2)
            plt.plot(data_buoy.t,rtrvl.snowice,linestyle='-', color='b',linewidth=2)
            plt.plot(data_buoy.t,rtrvl.congelation,linestyle='-', color='darkorange',linewidth=2)

        elif krun == 1:
            plt.sca(self.axH2)
            plt.plot(data_buoy.t,rtrvl.snowice,linestyle='-', color='c',linewidth=2)
            plt.plot(data_buoy.t,rtrvl.congelation,linestyle='-', color='gold',linewidth=2)

        plt.ylim(0.0, 20.0)
        plt.ylabel('cumulative growth (cm)')
        plt.grid(True)

        if krun == 1:
            plt.legend(('SIMBA1 snow-ice','SIMBA1 congelation','SIMBA2 snow-ice','SIMBA2 congelation'),bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=2, mode="expand", borderaxespad=0.,fontsize=10.0 )

        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        plt.text(-0.15, 1.2,'b)', ha='center', va='center', transform=self.axH2.transAxes,fontsize=14)

        nameFig5 = '%sFigure5_Thickness_Obs.png' % OutputFolder
        Fig5.savefig(nameFig5, dpi=600)
        if krun == 1:
            plt.close(Fig5)

#================================================================================================
# FIGURE 6 in Plante et al. 2023
#================================================================================================

    def figure_temperature_intercomp(self,
                                Data_buoy = None,
                                Data_model = None,
                                Data_mushy = None,
                                Rtrvl_buoy = None,
                                krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                OutputFolder = None):

        nameF6 = '%s_Ice_temperature.png' % OutputFolder

        #Top panels with simulated vertical temperature profiles, BL99 simulations
        #----------------------------------------------------------------

        Fig6 = plt.figure(self.Fig6Key)
        if krun == 0:
            ax1 = Fig6.add_axes([0.07, 0.55, 0.39, 0.32])
        elif krun == 1:
            ax1 = Fig6.add_axes([0.57, 0.55, 0.39, 0.32])

        z0 = Rtrvl_buoy.nsensors
        xtime = np.arange(0,Rtrvl_buoy.Ndata)

        ax1.set_facecolor("gray")
        data_p= np.ma.array(Data_model.data, mask=Data_model.data==np.nan)
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')
        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 =-50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data_model.t)

        dailytime = np.arange(len(Data_mushy.t))
        topsnow = -50.0+Data_model.Hsavg*50+z0
        topice = -50.0+Data_model.CumulSnowice*50.0+z0
        botice = -50.0-Data_model.Hiavg*50.0+z0+Data_model.CumulSnowice*50.0

        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax1.transAxes,fontsize=8)
        if krun == 0:
            plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'b)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)


        #Bottom panels with simulated vertical temperature profiles, mushy simulations
        #----------------------------------------------------------------

        if krun == 0:
            ax1 = Fig6.add_axes([0.07, 0.1, 0.39, 0.3])
        elif krun == 1:
            ax1 = Fig6.add_axes([0.57, 0.1, 0.39, 0.3])
        ax1.set_facecolor("gray")
        data_p= np.ma.array(Data_mushy.data, mask=Data_mushy.data==np.nan)
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')

        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 = -50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data_model.t)
        dailytime = np.arange(len(Data_mushy.t))

        topsnow = -50.0+Data_mushy.Hsavg*50+z0
        topice = -50.0+Data_mushy.CumulSnowice*50.0+z0
        botice = -50.0-Data_mushy.Hiavg*50.0+z0+Data_mushy.CumulSnowice*50.0
        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax1.transAxes,fontsize=8)

        if krun == 0:
            plt.text(-0.15, 1.2,'c)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'d)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
            Fig6.savefig(nameF6, dpi=600)
            plt.close(Fig6)
#================================================================================================
# FIGURE 7, 8 and 12 in Plante et al. 2023
#================================================================================================

    def F7_thickness_models(self, data_buoy = None,
                                  rtrvl = None,
                                  data_1 = None,
                                  data_2 = None,
                                  krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                  Figure = None, #This is the figure number (7, 8 or 12)
                                  OutputFolder=None):

        if Figure == 7:
            ExpLabel = "BL99"
            Fig= plt.figure(self.Fig7Key)
        elif Figure == 8:
            ExpLabel = "mushy"
            Fig= plt.figure(self.Fig8Key)
        if Figure == 12:
            ExpLabel = "mushy"
            Fig= plt.figure(self.Fig12Key)

        # Top pannels showing the ice thickness, snow thickness, and freeboard.  Thin lines are added if we have a no-snowice sim to compare with
        #----------------------------------------------------------------

        max_buoy = 1.0
        hfb_buoy_long = rtrvl.hi_int - (rtrvl.hs_int*330 + rtrvl.hi_int*917)/1026
        hfb_1 = data_1.Hiavg - (data_1.Hsavg*330 + data_1.Hiavg*917)/1026
        if data_2 is not None:
            hfb_2 = data_2.Hiavg - (data_2.Hsavg*330 + data_2.Hiavg*917)/1026

        plt.rc('legend',fontsize=10)
        if krun == 0:
            axH1b = Fig.add_axes([0.085, 0.55, 75.0/165.0, 0.3])
        elif krun == 1:
            axH3b = Fig.add_axes([0.64, 0.55, 55.0/165.0, 0.3])

        plt.plot(data_buoy.t,rtrvl.hi_int,linestyle='-', color='k',linewidth=2)
        plt.plot(data_buoy.t,rtrvl.hs_int,linestyle='-', color='dimgray',linewidth=2)
        plt.plot(data_buoy.t,hfb_buoy_long,linestyle='-', color='darkgrey',linewidth=2)


        max_buoy = np.nanmax([max_buoy, np.nanmax(rtrvl.hi_int)])
        plt.grid()

        plt.plot(data_1.t,data_1.Hiavg*100.0,linestyle='-', color='c',linewidth=2)
        plt.plot(data_1.t,data_1.Hsavg*100.0,linestyle='-', color='yellowgreen',linewidth=2)
        plt.plot(data_1.t,hfb_1*100.0,linestyle='-', color='gold',linewidth=2)

        if data_2 is not None:
            plt.plot(data_2.t,data_2.Hiavg*100.0,linestyle='-', color='c',linewidth=1)
            plt.plot(data_2.t,data_2.Hsavg*100.0,linestyle='-', color='yellowgreen',linewidth=1)
            plt.plot(data_2.t,hfb_2*100.0,linestyle='-', color='gold',linewidth=1)

        max_buoy = np.nanmax([max_buoy, np.nanmax(data_1.Hiavg)])
        maxlim = max_buoy*1.05
        plt.ylim(-10.0, 120)
        plt.ylabel('thickness (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.grid(True)
        if krun == 0:
            if data_2 is not None:
#                plt.text(0.1, 0.92,'SIMBA1', ha='center', va='center', transform=axH1b.transAxes,fontsize=10)
                plt.legend(('$h_\mathrm{i}$ SIMBA1','$h_\mathrm{s}$ SIMBA1','$h_\mathrm{fb}$ SIMBA1',
                            '$h_\mathrm{i}$ %s' % (ExpLabel),'$h_\mathrm{s}$ %s' % (ExpLabel),'$h_\mathrm{fb}$ %s' % (ExpLabel),
                            '$h_\mathrm{i}$ no s-i','$h_\mathrm{s}$  no s-i','$h_\mathrm{fb}$ no s-i'),
                bbox_to_anchor=(0.0, 1.04, 1., .102), loc='lower left', ncol=3, mode = 'expand', borderaxespad=0.,fontsize=6.0 )
            else:
                plt.legend(('$h_\mathrm{i}$ SIMBA1','$h_\mathrm{s}$ SIMBA1','$h_\mathrm{fb}$ SIMBA1',
                            '$h_\mathrm{i}$ %s' % (ExpLabel),'$h_\mathrm{s}$ %s' % (ExpLabel),'$h_\mathrm{fb}$ %s' % (ExpLabel)),
                bbox_to_anchor=(0., 1.04, 1., .102), loc='lower left', ncol=2,  borderaxespad=0.,fontsize=7.0 )
            plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=axH1b.transAxes,fontsize=12)
            (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        elif krun == 1:
            if data_2 is not None:
                plt.legend(('$h_\mathrm{i}$ SIMBA2','$h_\mathrm{s}$ SIMBA2','$h_\mathrm{fb}$ SIMBA2',
                            '$h_\mathrm{i}$ %s' % (ExpLabel),'$h_\mathrm{s}$ %s' % (ExpLabel),'$h_\mathrm{fb}$ %s' % (ExpLabel),
                            '$h_\mathrm{i}$ no s-i','$h_\mathrm{s}$  no s-i','$h_\mathrm{fb}$ no s-i'),
                bbox_to_anchor=(0.0, 1.04, 1., .102), loc='lower center', ncol=3, borderaxespad=0.,fontsize=6.0 )
#                plt.text(0.15, 0.92,'SIMBA2', ha='center', va='center', transform=axH3b.transAxes,fontsize=10)
            else:
                plt.legend(('$h_\mathrm{i}$ SIMBA2','$h_\mathrm{s}$ SIMBA2','$h_\mathrm{fb}$ SIMBA2',
                            '$h_\mathrm{i}$ %s' % (ExpLabel),'$h_\mathrm{s}$ %s' % (ExpLabel),'$h_\mathrm{fb}$ %s' % (ExpLabel)),
                bbox_to_anchor=(0., 1.04, 1., .102), loc='lower left', ncol=2,  borderaxespad=0.,fontsize=7.0 )

            plt.text(-0.17, 1.2,'b)', ha='center', va='center', transform=axH3b.transAxes,fontsize=12)
            (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,4,16])

        # Bottom pannels showing the congelation ice and snow-ice formation. Thin lines are added if we have a no-snowice sim to compare with
        #----------------------------------------------------------------

        if krun == 0:
            axH2b = Fig.add_axes([0.085, 0.075, 75.0/165.0, 0.3])
        elif krun == 1:
            axH4b = Fig.add_axes([0.64, 0.075, 55.0/165.0, 0.3])

        plt.plot(data_buoy.t,rtrvl.snowice,linestyle='-', color='k',linewidth=2)
        plt.plot(data_buoy.t,rtrvl.congelation,linestyle='-', color='dimgrey',linewidth=2)

        plt.plot(data_1.t,data_1.CumulSnowice*100.0,linestyle='-', color='c',linewidth=2)
        plt.plot(data_1.t,data_1.CongTotal*100.0,linestyle='-', color='gold',linewidth=2)
        if data_2 is not None:
            plt.plot(data_2.t,data_2.CumulSnowice*100.0,linestyle='-', color='c',linewidth=1)
            plt.plot(data_2.t,data_2.CongTotal*100.0,linestyle='-', color='gold',linewidth=1)

        plt.ylim(0.0, 20.0)
        plt.ylabel('cumulative growth (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.grid(True)

        if (krun == 0):
            if data_2 is not None:
#                plt.text(0.1, 0.92,'SIMBA1', ha='center', va='center', transform=axH2b.transAxes,fontsize=10)
                plt.legend(('$h_\mathrm{si}$ SIMBA1','$h_\mathrm{c}$ SIMBA1',
                            '$h_\mathrm{si}$ %s' % ExpLabel,'$h_\mathrm{c}$ %s' % ExpLabel,'$h_\mathrm{si}$ no s-i',
                            '$h_\mathrm{c}$ no s-i'),bbox_to_anchor=(0., 1.04, 1., .102), loc='lower center', ncol=3, mode = 'expand', borderaxespad=0.,fontsize=6.0 )
            else:
                plt.legend(('$h_\mathrm{si}$ SIMBA1','$h_\mathrm{c}$ SIMBA1',
                            '$h_\mathrm{si}$ %s' % ExpLabel,'$h_\mathrm{c+f}$ %s' % ExpLabel),bbox_to_anchor=(0., 1.04, 1., .102), loc='lower left', ncol=2, borderaxespad=0.,fontsize=7.0 )

            (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
            plt.text(-0.15, 1.2,'c)', ha='center', va='center', transform=axH2b.transAxes,fontsize=12)

        elif (krun == 1):
            if data_2 is not None:
#                plt.text(0.15, 0.92,'SIMBA2', ha='center', va='center', transform=axH4b.transAxes,fontsize=10)
                plt.legend(('$h_\mathrm{si}$ SIMBA2','$h_\mathrm{c}$ SIMBA2',
                            '$h_\mathrm{si}$ %s' % ExpLabel,'$h_\mathrm{c}$ %s' % ExpLabel,'$h_\mathrm{si}$ no s-i',
                            '$h_\mathrm{c}$ no s-i'),bbox_to_anchor=(0., 1.04, 1., .102), loc='lower center', ncol=3, borderaxespad=0.,fontsize=6.0 )
            else:
                plt.legend(('$h_\mathrm{si}$ SIMBA2','$h_\mathrm{c}$ SIMBA2',
                            '$h_\mathrm{si}$ %s' % ExpLabel,'$h_\mathrm{c+f}$ %s' % ExpLabel),bbox_to_anchor=(0., 1.04, 1., .102), loc='lower left', ncol=2, borderaxespad=0.,fontsize=7.0 )

            (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,4,16])
            plt.text(-0.15, 1.2,'d)', ha='center', va='center', transform=axH4b.transAxes,fontsize=12)
            nameF7_8_12 = '%sFigure%s_simThickness.png' % (OutputFolder, Figure)
            Fig.savefig(nameF7_8_12, dpi=600)
            plt.close(Fig)

#================================================================================================
# FIGURE 9 in Plante et al. 2023
#================================================================================================

    def Figure9_temperature_intercomp(self,
                                Data_1 = None,
                                Data_2 = None,
                                Data_3 = None,
                                Rtrvl= None,
                                krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                OutputFolder = None):


        nameH2 = '%sFigure9_Ice_temperature.png' % OutputFolder
        Fig9 = plt.figure(figsize=[5.0, 6.0])

        #-------------------------------------
        # First pannel (Data_1)
        #-------------------------------------

        ax1 = Fig9.add_axes([0.15, 0.69, 0.75, 0.23])
        z0 = Rtrvl.nsensors
        dailytime = np.arange(len(Data_1.t))

        #Make vertical profiles pcolor
        ax1.set_facecolor("gray")
        data_p= np.ma.array(Data_1.data, mask=Data_1.data==np.nan)
        a,b = data_p.shape
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')
        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 = -50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,4,17],time_ax = Data_1.t)

        # Draw interfaces
        topsnow = -50.0+Data_1.Hsavg*50+z0
        topice = -50.0+Data_1.CumulSnowice*50.0+z0
        botice = -50.0-Data_1.Hiavg*50.0+z0+Data_1.CumulSnowice*50.0
        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')

        # Add legends and labels
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax1.transAxes,fontsize=8)
        plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)

        #-------------------------------------
        # Second pannel (Data_2)
        #-------------------------------------

        ax2 = Fig9.add_axes([0.15, 0.37, 0.75, 0.23])
        z0 = Rtrvl.nsensors
        dailytime = np.arange(len(Data_2.t))

        #Make vertical profiles pcolor
        ax2.set_facecolor("gray")
        data_p= np.ma.array(Data_2.data, mask=Data_2.data==np.nan)
        a,b = data_p.shape
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')
        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 = -50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,4,17],time_ax = Data_2.t)

        # Draw interfaces
        topsnow = -50.0+Data_2.Hsavg*50+z0
        topice = -50.0+Data_2.CumulSnowice*50.0+z0
        botice = -50.0-Data_2.Hiavg*50.0+z0+Data_2.CumulSnowice*50.0
        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')

        # Add legends and labels
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax2)
        cax2 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax2,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax2.transAxes,fontsize=8)
        plt.text(-0.15, 1.2,'b)', ha='center', va='center', transform=ax2.transAxes,fontsize=12)

        #-------------------------------------
        # Third pannel (Data_3)
        #-------------------------------------

        ax3 = Fig9.add_axes([0.15, 0.05, 0.75, 0.23])
        z0 = Rtrvl.nsensors
        dailytime = np.arange(len(Data_3.t))

        #Make vertical profiles pcolor
        ax3.set_facecolor("gray")
        data_p= np.ma.array(Data_3.data, mask=Data_3.data==np.nan)
        a,b = data_p.shape
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')
        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 =-50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,4,17],time_ax = Data_3.t)

        # Draw interfaces
        topsnow = -50.0+Data_3.Hsavg*50+z0
        topice = -50.0+Data_3.CumulSnowice*50.0+z0
        botice = -50.0-Data_3.Hiavg*50.0+z0+Data_3.CumulSnowice*50.0
        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')

        # Add legends and labels
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax3)
        cax3 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax3,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax3.transAxes,fontsize=8)
        plt.text(-0.15, 1.2,'c)', ha='center', va='center', transform=ax3.transAxes,fontsize=12)


        Fig9.savefig(nameH2, dpi=600)
        plt.close(Fig9)

#================================================================================================
# FIGURE 10 in Plante et al. 2023
#================================================================================================

    def figure_layer_thermodynamics(self,
                                Data_buoy = None,
                                Data_1 = None,
                                Data_2 = None,
                                Data_3 =  None,
                                N = None, #this is the layer number
                                krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                OutputFolder = None):

        name2 = '%sInternal_layer_%s_Thermo_.png' % (OutputFolder,N)

        Fig10 = plt.figure(figsize=[5.0, 5.0])
        plt.rc('legend',fontsize=12)

        #-----------------------------
        # First pannel: layer temperature
        #-----------------------------

        ax1 = Fig10.add_axes([0.15, 0.7, 0.8, 0.22])
        plt.plot(Data_1.t[:],Data_1.Tice[:,N],linestyle='-', color = 'b', label='no flooding')
        plt.plot(Data_2.t[:],Data_2.Tice[:,N],linestyle='-', color = 'darkgreen', label='$\phi_\mathrm{min}$ = 0.005')
        plt.plot(Data_3.t[:],Data_3.Tice[:,N],linestyle='-', color = 'darkorange', label='Manual')

        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        plt.grid()
        plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc='lower left', ncol=3, borderaxespad=0.,fontsize=7.0 )
        if N == 0:
            plt.ylim(-6.0, 0.0) #For the top layer
        elif N == 6:
            plt.ylim(-2.5, -1.8) #For the bottom layer
        else:
            pass
        plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        plt.ylabel('T$_\mathrm{lyr}$ ($^\circ$C)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)


        #-----------------------------
        # Second pannel: bulk salinity
        #-----------------------------

        ax2 = Fig10.add_axes([0.15, 0.38, 0.8, 0.22])
        plt.plot(Data_1.t[:],Data_1.Sice[:,N],linestyle='-', color = 'b', label='no flooding')
        plt.plot(Data_2.t[:],Data_2.Sice[:,N],linestyle='-', color = 'green', label='$\phi_\mathrm{min}$ = 0.005')
        plt.plot(Data_3.t[:],Data_3.Sice[:,N],linestyle='-', color = 'darkorange', label='Manual')

        plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc='lower left', ncol=3, borderaxespad=0.,fontsize=7.0 )
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        if N == 0:
            plt.ylim(0.0, 15.0)
        elif N == 6:
            plt.ylim(0.0, 20.0)
        else:
            pass
        plt.text(-0.15, 1.2,'b)', ha='center', va='center', transform=ax2.transAxes,fontsize=12)
        plt.grid()
        plt.ylabel('S$_\mathrm{lyr}$ (PSU)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)

        #-----------------------------
        # Third pannel: Brine salinity
        #-----------------------------

        ax3 = Fig10.add_axes([0.15, 0.06, 0.8, 0.22])

        plt.plot(Data_1.t[:],Data_1.Sice[:,N]/Data_1.LiqFraction[:,N],linestyle='-', color = 'b', label='no flooding')
        plt.plot(Data_2.t[:],Data_2.Sice[:,N]/Data_2.LiqFraction[:,N],linestyle='-', color = 'green', label='$\phi_\mathrm{min}$ = 0.005')
        plt.plot(Data_3.t[:],Data_3.Sice[:,N]/Data_3.LiqFraction[:,N],linestyle='-', color = 'darkorange', label='Manual')

        plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc='lower left', ncol=3, borderaxespad=0.,fontsize=7.0 )
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])

        if N == 0:
            plt.ylim(0.0, 125.0)
        elif N == 6:
            plt.ylim(25.0, 45.0)
        else:
            pass

        plt.grid()
        plt.text(-0.15, 1.2,'c)', ha='center', va='center', transform=ax3.transAxes,fontsize=12)
        plt.ylabel('$S_\mathrm{br}$ (PSU)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)


        Fig10.savefig(name2, dpi=600)
        plt.close(Fig10)

#================================================================================================
# FIGURE 11 in Plante et al. 2023
#================================================================================================

    def Figure11_bottomlayer_thermodynamics(self,
                                Data_buoy = None,
                                Data_1 = None,
                                Data_2 = None,
                                Data_3 =  None,
                                krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                OutputFolder = None):


        nameF11 = '%sBottom_layer_Thermo.png' % (OutputFolder)

        Fig11 = plt.figure(self.F11Key)
        plt.rc('legend',fontsize=12)

        #-----------------------------
        # First pannel: layer temperature
        #-----------------------------
        if krun == 0:
            ax1 = Fig11.add_axes([0.08, 0.77, 0.4, 0.16])
            plt.plot(Data_1.t[:],Data_1.Tice[:,6],linestyle='-', color = 'b', label='$\omega$=-5.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_2.t[:],Data_2.Tice[:,6],linestyle='-', color = 'darkgreen', label='$\omega$=-2.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_3.t[:],Data_3.Tice[:,6],linestyle='-', color = 'darkorange', label='$\omega$=-1.0e$^{-9}$ m s$^{-1}$')

        elif krun == 1:
            ax1 = Fig11.add_axes([0.58, 0.77, 0.4, 0.16])

            plt.plot(Data_1.t[:],Data_1.Tice[:,6],linestyle='-', color = 'b', label='$\phi_\mathrm{i}$=0.85')
            plt.plot(Data_2.t[:],Data_2.Tice[:,6],linestyle='-', color = 'darkgreen', label='$\phi_\mathrm{i}$=0.65')
            plt.plot(Data_3.t[:],Data_3.Tice[:,6],linestyle='-', color = 'darkorange', label='$\phi_\mathrm{i}$=0.45')

        plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc='lower left', ncol=3, borderaxespad=0.,fontsize=7.0 )
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        plt.grid()

        if krun == 0:
            plt.ylim(-2.4, -1.4)
            plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        elif krun == 1:
            plt.ylim(-2.4, -1.4)
            plt.text(-0.15, 1.2,'b)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)

        plt.ylabel('T$_\mathrm{lyr}$ ($^\circ$C)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)



        #-----------------------------
        # Second pannel: bulk salinity
        #-----------------------------
        if krun == 0:
            ax2 = Fig11.add_axes([0.08, 0.53, 0.4, 0.16])
            plt.plot(Data_1.t[:],Data_1.Sice[:,6],linestyle='-', color = 'b', label='$\omega$=-5.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_2.t[:],Data_2.Sice[:,6],linestyle='-', color = 'green', label='$\omega$=-2.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_3.t[:],Data_3.Sice[:,6],linestyle='-', color = 'darkorange', label='$\omega$=-1.0e$^{-9}$ m s$^{-1}$')
        elif krun == 1:
            ax2 = Fig11.add_axes([0.58, 0.53, 0.4, 0.16])
            plt.plot(Data_1.t[:],Data_1.Sice[:,6],linestyle='-', color = 'b', label='$\phi_\mathrm{i}$=0.85')
            plt.plot(Data_2.t[:],Data_2.Sice[:,6],linestyle='-', color = 'green', label='$\phi_\mathrm{i}$=0.65')
            plt.plot(Data_3.t[:],Data_3.Sice[:,6],linestyle='-', color = 'darkorange', label='$\phi_\mathrm{i}$=0.45')

        plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc='lower left', ncol=3, borderaxespad=0.,fontsize=7.0 )
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])

        if krun == 0:
            plt.ylim(0, 25)
            plt.text(-0.15, 1.2,'c)', ha='center', va='center', transform=ax2.transAxes,fontsize=12)
        elif krun == 1:
            plt.ylim(0, 25)
            plt.text(-0.15, 1.2,'d)', ha='center', va='center', transform=ax2.transAxes,fontsize=12)

        plt.grid()
        plt.ylabel('S$_\mathrm{lyr}$ (PSU)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)

        #-----------------------------
        # Third pannel: Brine salinity
        #-----------------------------
        if krun == 0:
            ax3 = Fig11.add_axes([0.08, 0.29, 0.4, 0.16])
            plt.plot(Data_1.t[:],Data_1.Sice[:,6]/Data_1.LiqFraction[:,6],linestyle='-', color = 'b', label=          '$\omega$=-5.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_2.t[:],Data_2.Sice[:,6]/Data_2.LiqFraction[:,6],linestyle='-', color = 'green', label=      '$\omega$=-2.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_3.t[:],Data_3.Sice[:,6]/Data_3.LiqFraction[:,6],linestyle='-', color = 'darkorange', label= '$\omega$=-1.0e$^{-9}$ m s$^{-1}$')

        elif krun == 1:
            ax3 = Fig11.add_axes([0.58, 0.29, 0.4, 0.16])

            plt.plot(Data_1.t[:],Data_1.Sice[:,6]/Data_1.LiqFraction[:,6],linestyle='-', color = 'b', label='$\phi_\mathrm{i}$=0.85')
            plt.plot(Data_2.t[:],Data_2.Sice[:,6]/Data_2.LiqFraction[:,6],linestyle='-', color = 'green', label='$\phi_\mathrm{i}$=0.65')
            plt.plot(Data_3.t[:],Data_3.Sice[:,6]/Data_3.LiqFraction[:,6],linestyle='-', color = 'darkorange', label='$\phi_\mathrm{i}$=0.45')

        plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc='lower left', ncol=3, borderaxespad=0.,fontsize=7.0 )
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])

        if krun == 0:
            plt.ylim(25, 45)
            plt.text(-0.15, 1.2,'e)', ha='center', va='center', transform=ax3.transAxes,fontsize=12)
        elif krun == 1:
            plt.ylim(25, 45)
            plt.text(-0.15, 1.2,'f)', ha='center', va='center', transform=ax3.transAxes,fontsize=12)

        plt.grid()
        plt.ylabel('$S_\mathrm{br}$ (PSU)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)

        #-----------------------------
        # Fourth pannel: Brine drainage
        #-----------------------------

        if krun == 0:
            ax4 = Fig11.add_axes([0.08, 0.05, 0.4, 0.16])
            plt.plot(Data_1.t[:], Data_1.dSdt[:,6],linestyle='-', color = 'b', label='$\omega$=-5.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_2.t[:], Data_2.dSdt[:,6],linestyle='-', color = 'green', label='$\omega$=-2.0e$^{-9}$ m s$^{-1}$')
            plt.plot(Data_3.t[:], Data_3.dSdt[:,6],linestyle='-', color = 'darkorange', label='$\omega$=-1.0e$^{-9}$ m s$^{-1}$')

        elif krun == 1:
            ax4 = Fig11.add_axes([0.58, 0.05, 0.4, 0.16])
            plt.plot(Data_1.t[:], Data_1.dSdt[:,6],linestyle='-', color = 'b', label='$\phi_\mathrm{i}$=0.85')
            plt.plot(Data_2.t[:], Data_2.dSdt[:,6],linestyle='-', color = 'green', label='$\phi_\mathrm{i}$=0.65')
            plt.plot(Data_3.t[:], Data_3.dSdt[:,6],linestyle='-', color = 'darkorange', label='$\phi_\mathrm{i}$=0.45')

        plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc='lower left', ncol=3, borderaxespad=0.,fontsize=7.0 )


        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        if krun == 0:
            plt.ylim(-1.5, 0.1)
            plt.text(-0.15, 1.2,'g)', ha='center', va='center', transform=ax4.transAxes,fontsize=12)
        elif krun == 1:
            plt.ylim(-1.5, 0.1)
            plt.text(-0.15, 1.2,'h)', ha='center', va='center', transform=ax4.transAxes,fontsize=12)

        plt.grid()
        plt.ylabel('Brine drainage rate',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        if krun == 1:
            Fig11.savefig(nameF11, dpi=600)
            plt.close(Fig11)


#================================================================================================
# FIGURES in Appendix, in Plante et al. 2023
#================================================================================================

    def figure_congelation_intercomp(self, Data_buoy = None,
                                           Rtrvl = None,
                                           Data_1 = None,
                                           Data_2 = None,
                                           Data_3 = None,
                                           krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                           Figure = None,
                                           OutputFolder=None):

        plt.rc('legend',fontsize=12)

        FigApp = plt.figure(figsize=[6.0,6.0])
        ax1 =FigApp.add_axes([0.15, 0.1, 0.8, 0.6] )

        #Find the congelation already present in Obs at the model t=0. All model data have same start time.
        first_mod_time = Data_1.t[0]
        ntime = np.arange(len(np.squeeze(Data_buoy.t)))
        time_obs = ntime[Data_buoy.t[:]>first_mod_time]
        cong_init = Rtrvl.congelation[int(time_obs[0])-1]
        plt.ylim(0, 20.0)
        plt.ylabel('basal growth (cm)')
        plt.grid()
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5])
        plt.plot(Data_buoy.t,Rtrvl.congelation-cong_init,linestyle='-', color='k',linewidth=2)

        plt.plot(Data_1.t,Data_1.CongTotal*100.0,linestyle='-', color='blue',linewidth=2)
        plt.plot(Data_1.t,Data_1.CumulCong*100.0-Data_1.CumulBotMelt*100.0,linestyle='--', color='b',linewidth=1)
        plt.plot(Data_1.t,Data_1.CumulNewIce*100.0,linestyle='-.', color='b',linewidth=1)
        plt.plot(np.NaN, np.NaN, '-', color='none')
        plt.plot(Data_2.t,Data_2.CongTotal*100.0,linestyle='-', color='green',linewidth=2,zorder=99)
        plt.plot(Data_2.t,Data_2.CumulCong*100.0-Data_2.CumulBotMelt*100.0,linestyle='--', color='green',linewidth=1,zorder=100)
        plt.plot(Data_2.t,Data_2.CumulNewIce*100.0,linestyle='-.', color='green',linewidth=1,zorder=101)

        if Data_3 is not None:
            print("Making bugged figure...")
            plt.plot(np.NaN, np.NaN, '-', color='none')
            plt.plot(Data_3.t,Data_3.CongTotal*100.0,linestyle='-', color='darkorange',linewidth=2,zorder=102)
            plt.plot(Data_3.t,Data_3.CumulCong*100.0-Data_3.CumulBotMelt*100.0,linestyle='--', color='darkorange',linewidth=1,zorder=104)
            plt.plot(Data_3.t,Data_3.CumulNewIce*100.0,linestyle='-.', color='darkorange',linewidth=1,zorder=106)


        if Figure == 3:
            plt.legend(('SIMBA1','Total, $\phi_\mathrm{init}$=0.85','Cong., $\phi_\mathrm{init}$=0.85',
                        'Frazil, $\phi_\mathrm{init}$=0.85','','Total, $\phi_\mathrm{init}$=0.65',
                        'Cong., $\phi_\mathrm{init}$=0.65','Frazil, $\phi_\mathrm{init}$=0.65','',
                        'Total, $\phi_\mathrm{init}$=0.45','Cong., $\phi_\mathrm{init}$=0.45','Frazil, $\phi_\mathrm{init}$=0.45'),
                                 bbox_to_anchor=(0.0, 1.05, 1., .102), loc='lower left',
                                 ncol=3, mode = "expand", borderaxespad=0.,fontsize=9.0 )
        elif Figure == 2:
            plt.legend(('SIMBA1','Total, std','Cong., std','Frazil, std','','Total, modified','Cong., modified','Frazil,  modified'), bbox_to_anchor=(0.0, 1.05, 1., .102), loc='lower left',
                         ncol=2, mode = "expand", borderaxespad=0.,fontsize=9.0 )
        else:
            plt.legend(('SIMBA1','Total, $\phi_\mathrm{init}$=0.85','Cong., $\phi_\mathrm{init}$=0.85','Frazil, $\phi_\mathrm{init}$=0.85','',
                                 'Total, $\phi_\mathrm{init}$=0.65','Cong., $\phi_\mathrm{init}$=0.65','Frazil, $\phi_\mathrm{init}$=0.65'), bbox_to_anchor=(0.0, 1.05, 1., .102), loc='lower left',
                         ncol=2, mode = "expand", borderaxespad=0.,fontsize=9.0 )

        nameFigApp = '%sCongelation_Intercomparison.png' % (OutputFolder)
        FigApp.savefig(nameFigApp)
        plt.close(FigApp)

#================================================================================================
# FIGURE Supplementary in Plante et al. 2023
#================================================================================================

    def figure_temperature_intercomp_supp(self,
                                Data_buoy = None,
                                Data_model = None,
                                Data_mushy = None,
                                Data_3 = None,
                                Rtrvl_buoy = None,
                                krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                                OutputFolder = None):

        nameFSup = '%s_Ice_temperature.png' % OutputFolder

        #Top panels with simulated vertical temperature profiles, BL99 simulations
        #----------------------------------------------------------------

        FigSup = plt.figure(self.FigSupKey)
        ax1 = FigSup.add_axes([0.07+krun*0.5, 0.68, 0.38, 0.2])
        z0 = Rtrvl_buoy.nsensors
        xtime = np.arange(0,Rtrvl_buoy.Ndata)

        ax1.set_facecolor("gray")
        data_p= np.ma.array(Data_model.data, mask=Data_model.data==np.nan)
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')
        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 =-50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data_model.t)

        dailytime = np.arange(len(Data_model.t))
        topsnow = -50.0+Data_model.Hsavg*50+z0
        topice = -50.0+Data_model.CumulSnowice*50.0+z0
        botice = -50.0-Data_model.Hiavg*50.0+z0+Data_model.CumulSnowice*50.0

        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax1.transAxes,fontsize=8)
        if krun == 0:
            plt.text(-0.15, 1.2,'a)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'b)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)


        #Middle panels with simulated vertical temperature profiles, mushy simulations
        #----------------------------------------------------------------

        
        ax1 = FigSup.add_axes([0.07+krun*0.5, 0.38, 0.38, 0.2])
        ax1.set_facecolor("gray")
        data_p= np.ma.array(Data_mushy.data, mask=Data_mushy.data==np.nan)
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')

        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 = -50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data_model.t)
        dailytime = np.arange(len(Data_mushy.t))

        topsnow = -50.0+Data_mushy.Hsavg*50+z0
        topice = -50.0+Data_mushy.CumulSnowice*50.0+z0
        botice = -50.0-Data_mushy.Hiavg*50.0+z0+Data_mushy.CumulSnowice*50.0
        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax1,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax1.transAxes,fontsize=8)

        if krun == 0:
            plt.text(-0.15, 1.2,'c)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'d)', ha='center', va='center', transform=ax1.transAxes,fontsize=12)
            
        #-------------------------------------
        # Third pannel (Data_3)
        #-------------------------------------

        ax3 = FigSup.add_axes([0.07+krun*0.5, 0.08, 0.38, 0.2])
        ax3.set_facecolor("gray")
        data_p= np.ma.array(Data_3.data, mask=Data_3.data==np.nan)
        cmap = mpl.colormaps["Spectral_r"]
        cmap.set_bad('gray')

        im = plt.imshow(np.flipud(np.transpose(data_p)),cmap=cmap,aspect='auto',norm=mpl.colors.TwoSlopeNorm(vmin=-26., vcenter=-5., vmax=-2.))
        self.make_y_axis(Zice0 = -50.0+z0)
        plt.ylabel('vertical position (cm)',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        (lim1,lim2) = self.make_time_labels([2017,2,23],[2017,5,5],time_ax = Data_3.t)
        dailytime = np.arange(len(Data_3.t))

        topsnow = -50.0+Data_3.Hsavg*50+z0
        topice = -50.0+Data_3.CumulSnowice*50.0+z0
        botice = -50.0-Data_3.Hiavg*50.0+z0+Data_3.CumulSnowice*50.0
        plt.plot(dailytime,topsnow,linestyle='-',color='b', label='snow surface')
        plt.plot(dailytime,dailytime*0.0-50.0+z0,linestyle='--',color='r', label='init. ice sfc.')
        plt.plot(dailytime,topice,linestyle='-',color='k', label='ice interface')
        plt.plot(dailytime,botice,linestyle='-',color='k', label='')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=4, borderaxespad=0.,fontsize=7.0 )
        divider = make_axes_locatable(ax3)
        cax3 = divider.append_axes("right", size="5%",pad=0.2)
        clb = plt.colorbar(im,cax=cax3,ticks=[-26.0, -19.0, -12.0, -5.0, -4.0, -3.0,-2.0])#,shrink=0.5
        clb.ax.tick_params(labelsize=8)
        plt.text(0.4,1.07,r'T ($^\circ C$)',ha='center', va='center', transform=cax3.transAxes,fontsize=8)       
           
        if krun == 0:
            plt.text(-0.15, 1.2,'e)', ha='center', va='center', transform=ax3.transAxes,fontsize=12)
        elif krun == 1:
            plt.text(-0.15, 1.2,'f)', ha='center', va='center', transform=ax3.transAxes,fontsize=12) 
            FigSup.savefig(nameFSup, dpi=600)
            plt.close(FigSup)


#================================================================================================
# Computing errors (MIE, RMSE), in Plante et al. 2023
#================================================================================================

    def Compute_Errors(self, Data_buoy = None,
                             Data_model = None,
                             variable = None,
                             reference = None,
                             krun = None, #This is the Exp itterate, 0=IMB1, 1=IMB2
                             OutputFolder=None):
    ##-------------------------------------------------------------------------------------------------
    # Plotting the air temperature in the simulations are in the Buoys
    ###-------------------------------------------------------------------------------------------------


        plt.rc('legend',fontsize=12)
        #Find the time interval in common
        t_min_common= np.max([Data_buoy.t[0], Data_model.t[0]])
        t_max_common = np.nanmin([np.nanmax(Data_buoy.t), np.nanmax(Data_buoy.t)])
        t_length = ( t_max_common-t_min_common )

        #Re-arrange in hourly data
        t_common = np.arange(t_min_common, t_max_common, 1.0/24.0)


        DeltaT_model = Data_model.t[0]-Data_buoy.t[0]
        DeltaT_buoy = Data_buoy.t[0]-Data_model.t[0]
        DeltaT_min = np.min([DeltaT_model,DeltaT_buoy])

        indices_model = np.arange(0,len(Data_model.t[:]))
        indices_buoy = np.arange(0,len(Data_buoy.t[:]))

        #prepare vectors that will be later populated
        Tdiff = t_common.copy()*np.nan
        Tbuoy = Tdiff.copy()
        Tmodel= Tdiff.copy()

        for k in range(0,len(t_common)):
            t = t_common[k]

            #Get data at corresponding point
            ind_model = indices_model[Data_model.t==t]
            if not ind_model: #in this case, we interpolate from neighboring data points
                if t<Data_model.t[0] or t > np.nanmax(Data_model.t[:]):
                    T_model_loc = np.nan
                else:
					#Find the existing time closest to t in the data
                    dist = (Data_model.t[:]-t)**2.0
                    ind_closest = indices_model[dist[:] == np.nanmin(dist)]
                    #Get the indices of the neighbouring points
                    if Data_model.t[ind_closest] > t:
                        ind1m = ind_closest - 1
                        ind2m = ind_closest
                    else:
                        ind2m = ind_closest + 1
                        ind1m = ind_closest

					#Get the data at time t by interpolation
                    T_model_loc = variable[ind1m]+(t-Data_model.t[ind1m]) *(variable[ind2m] -
                                            variable[ind1m]) / (Data_model.t[ind2m]
					                                     - Data_model.t[ind1m])
            else:
                T_model_loc = variable[ind_model]
            #Populate the model data vector
            Tmodel[k] = T_model_loc


            #Repeat for the Observed values
            t = t_common[k]
            if np.isnan(t):
                Tdummy = np.nan
            else:
                ind_buoy = indices_buoy[Data_buoy.t==t]
                if not ind_buoy:
                    if t<Data_buoy.t[0] or t > np.nanmax(Data_buoy.t[:]):
                        T_buoy_loc = np.nan
                    else:
                        dist = (Data_buoy.t[:]-t)**2.0
                        ind_closest = indices_buoy[dist[:] == np.nanmin(dist)]
                        ind_closest = ind_closest[0]
                        if Data_buoy.t[ind_closest] > t:
                            ind1m = ind_closest - 1
                            ind2m = ind_closest
                        else:
                            ind2m = ind_closest + 1
                            ind1m = ind_closest

                        T_buoy_loc = reference[ind1m]+(t-Data_buoy.t[ind1m]) *(reference[ind2m] -
					      reference[ind1m]) / (Data_buoy.t[ind2m]
					                                     - Data_buoy.t[ind1m])
                else:
                    T_buoy_loc = reference[ind_buoy]
                #Populate the Obs data vector
                Tbuoy[k] = T_buoy_loc

        #Calculate the errors
        Tdiff = Tmodel-Tbuoy
        RMSD = (np.nanmean((Tdiff)**2.0))*0.5
        MIE = np.nanmean(Tdiff)

        #Print the errors in a txt document output
        with open('%sPrinted_error.txt' % OutputFolder, 'w') as f:
            f.write('MIE = %s, RMSD = %s\n' % (MIE,RMSD))
        return


class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        print(majorlocs)
        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in np.arange(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)
        print(minorlocs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


