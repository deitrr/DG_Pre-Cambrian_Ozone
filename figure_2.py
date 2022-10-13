import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy.interpolate as si
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import pdb
#import pathlib
import matplotlib.gridspec as gridspec

#this just suppresses a deprecation warning from netCDF4
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#list of simulations that will go into figure 2
simnames =["const_CO2_O2_1e-9",
           "const_CO2_O2_1e-7",
           "const_CO2_O2_1e-4",
           "const_CO2_O2_1e-3",
           "const_CO2_O2_1e-2",
           "ref_const_O3",
           "temp_cont_O2_1e-9",
           "temp_cont_O2_1e-7",
           "temp_cont_O2_1e-4",
           "temp_cont_O2_1e-3",
           "const_CO2_O2_1e-2",
           "ref_const_O3"]
#o2_title = ['Pre-GOE           ',
#            'Start of GOE      ',
#            'After GOE         ',
#            'Proterozoic 1     ',
#            'Proterozoic 2     ',
#            'PAL (constant O3) ',
#            'PAL (variable O3) ',
#            'Start GOE (hiCO2) ']
#fields = ['TS','LWCF','SWCF','CLDTOT','CLDHGH','CLDMED','CLDLOW','TGCLDCWP']

#nominal O2 levels (use as x-axis)
o2s = np.array([4.5545e-9,6.059e-8,7.255e-5,9.082e-4,9.941e-3,0.2112,4.5545e-9,6.059e-8,7.255e-5,9.082e-4,9.941e-3,0.2112])

#create the arrays
Ts = np.zeros(len(simnames))
o3col = np.zeros(len(simnames))
rh500 = np.zeros(len(simnames))
q500 = np.zeros(len(simnames))
q20 = np.zeros(len(simnames))
rh20 = np.zeros(len(simnames))
T20 = np.zeros(len(simnames))
cldlow = np.zeros_like(Ts)
cldhgh = np.zeros_like(Ts)
#cldmed = np.zeros_like(Ts)
#cwp = np.zeros_like(Ts)
icwp = np.zeros_like(Ts)
lcwp = np.zeros_like(Ts)
#cwp_low = np.zeros_like(Ts)
#cwp_med = np.zeros_like(Ts)
#cwp_high = np.zeros_like(Ts)
lwcf = np.zeros_like(Ts)
swcf = np.zeros_like(Ts)
icldt = np.zeros_like(Ts)
#lcldt = np.zeros_like(Ts)
#tcldt = np.zeros_like(Ts)
icldp = np.zeros_like(Ts)
#lcldp = np.zeros_like(Ts)
#tcldp = np.zeros_like(Ts)
lts = np.zeros_like(Ts)
eis = np.zeros_like(Ts)
#eis_frac = np.zeros_like(Ts)
Ttroptrop = np.zeros_like(Ts)
eqTs = np.zeros_like(Ts)
plTs = np.zeros_like(Ts)

#
#R = 287.0
#g = 9.80665

for i in np.arange(len(simnames)):
  print(simnames[i])

  #Global mean file----------------------
  simpath = "CAM_post_process/" + simnames[i]
  file = simpath + "/" + simnames[i] + ".cam.h0.globmean_0031_0060.nc"
  data = nc.Dataset(file,'r')

  #read in data from global mean
  p = data['lev'][:].squeeze()
  T = data['T'][:].squeeze()
  Ts[i] = data['TS'][:].squeeze()
  o3col[i] = data['cb_ozone_c'][:].squeeze()
  rh = data['RELHUM'][:].squeeze()
  q = data['Q'][:].squeeze()
  cldlow[i] = data['CLDLOW'][:].squeeze()
  #cldmed[i] = data['CLDMED'][:].squeeze()
  cldhgh[i] = data['CLDHGH'][:].squeeze()
#  cwp[i] = data['TGCLDCWP'][:].squeeze()
  lwcf[i] = data['LWCF'][:].squeeze()
  swcf[i] = data['SWCF'][:].squeeze()
  icwp[i] = data['TGCLDIWP'][:].squeeze()
  lcwp[i] = data['TGCLDLWP'][:].squeeze()

  #interpolate these to 500 hPa
  rh500[i] = si.interp1d(p,rh)(500.0)
  q500[i] = si.interp1d(p,q)(500.0)

  #interpolate these to 20 hPa
  q20[i] = si.interp1d(p,q)(20.0)
  rh20[i] = si.interp1d(p,rh)(20.0)
  T20[i] = si.interp1d(p,T)(20.0)


  #Global mean cloud tops----------------------
  file_cldt = simpath + "/" + simnames[i] + ".cam.h0.globmean_0031_0060_cldtp.nc"
  data_cldt = nc.Dataset(file_cldt,'r')

  #read in cloud top data
  icldt[i] = data_cldt['ICLDT'][:].squeeze()
  #lcldt[i] = data_cldt['LCLDT'][:].squeeze()
  #tcldt[i] = data_cldt['TCLDT'][:].squeeze()
  icldp[i] = data_cldt['ICLDP'][:].squeeze()
  #lcldp[i] = data_cldt['LCLDP'][:].squeeze()
  #tcldp[i] = data_cldt['TCLDP'][:].squeeze()
  #icldt[i] = data['ICLDT'][:].squeeze()
  #lcldt[i] = data['LCLDT'][:].squeeze()
  #tcldt[i] = data['TCLDT'][:].squeeze()
  #icldp[i] = data['ICLDP'][:].squeeze()
  #lcldp[i] = data['LCLDP'][:].squeeze()
  #tcldp[i] = data['TCLDP'][:].squeeze()



  #p_int = data_tim.variables['ilev'][:].squeeze()*100.0
  #rho = p*100 / (R * T)
  #dz = np.diff(p_int) / (rho*g)

  #iwc = data['IWC'][:].squeeze()
  #lwc = data['LWC'][:].squeeze()
  #cwp_tmp = (iwc + lwc) * dz
  #if 'ICLDTWP' in data.variables:
   # pdb.set_trace()
  #cwp_low[i] = np.sum(cwp_tmp[p>=700])
  #cwp_med[i] = np.sum(cwp_tmp[np.logical_and(p<700,p>=400)])
  #cwp_high[i] = np.sum(cwp_tmp[np.logical_and(p<400,p>=50)])


  #Tropical mean file--------------------------------
  file_trop = simpath + "/" + simnames[i] + ".cam.h0.tropmean_0031_0060.nc"
  data_trop = nc.Dataset(file_trop,'r')

  #Get data I need from tropical mean
  Ttroptrop[i] = si.interp1d(p,data_trop['T'][:].squeeze())(100.0)
  eqTs[i] = data_trop['TS'][:].squeeze()


  #Polar mean files---------------------------------
  file_spol = simpath + "/" + simnames[i] + ".cam.h0.spolmean_0031_0060.nc"
  data_spol = nc.Dataset(file_spol,'r')
  file_npol = simpath + "/" + simnames[i] + ".cam.h0.npolmean_0031_0060.nc"
  data_npol = nc.Dataset(file_npol,'r')
  plTs[i] = 0.5*(data_spol['TS'][:].squeeze()+data_npol['TS'][:].squeeze())

  #file_lts = simpath+ "/" + simnames[i] + ".cam.h0.globmean_0031_0060_lts.nc"
  #data_lts = nc.Dataset(file_lts,'r')
  #lts[i] = data_lts['LTS'][:].squeeze()

  #file_eis = simpath+ "/" + simnames[i] + ".cam.h0.globmean_0031_0060_eis.nc"
  #data_eis = nc.Dataset(file_eis,'r')
  #eis[i] = data_eis['EIS'][:].squeeze()


  #Need a few things from time mean file-------------
  file_tim = simpath + "/" + simnames[i] + ".cam.h0.timmean_0031_0060.nc"
  data_tim = nc.Dataset(file_tim,'r')
  landfrac = data_tim['LANDFRAC'][:].squeeze()
  lats = data_tim['lat'][:].squeeze()


  #LTS file---------------------------------------
  file_lts = simpath+ "/" + simnames[i] + ".cam.h0.timmean_0031_0060_lts.nc"
  data_lts = nc.Dataset(file_lts,'r')
  lts_arr = data_lts['LTS'][:].squeeze()
  #Horizontal average over marine regions only
  lts[i] = np.sum((lts_arr*np.cos(lats*np.pi/180)[:,None]*(1-landfrac)))/np.sum(np.cos(lats*np.pi/180))/144
  #lts[i] = data_lts['LTS'][:].squeeze()


  #EIS file---------------------------------------
  file_eis = simpath+ "/" + simnames[i] + ".cam.h0.timmean_0031_0060_eis.nc"
  data_eis = nc.Dataset(file_eis,'r')
  eis_arr = data_eis['EIS'][:].squeeze()
  #Horizontal average over marine regions and EIS>0 only
  eis[i] = np.sum((eis_arr*np.cos(lats*np.pi/180)[:,None]*(1-landfrac))[eis_arr>0])/np.sum(np.cos(lats*np.pi/180))/144
  #eis[i] = data_eis['EIS'][:].squeeze()
  #pdb.set_trace()

  #eis_frac[i] = np.sum(((1-landfrac)*np.cos(lats[:,None]*np.pi/180))[eis_arr>0])/np.sum((1-landfrac)*np.cos(lats[:,None]*np.pi/180))


#this will be the synthesis plot
plt.rcParams.update({'font.size':6})
fig = plt.figure(figsize=(7.5,5))

#Set up the plot arrangement 
outer_grid = gridspec.GridSpec(5,1,wspace=0.2,hspace=0.3,left=0.07,right=0.98,bottom=0.07,top=0.96,height_ratios=(2,2,2,2,2))
row1 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[0],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
row2 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[1],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
row3 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[2],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
row4 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[3],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))
row5 = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_grid[4],wspace=0.5,hspace=0.1,width_ratios=(1,1,1,1))

#top row
ax = fig.add_subplot(row1[0])
ax.plot(o2s[:6],Ts[:6],color='k',marker='.',linestyle='--',label='Constant CO$_2$',lw=2)
ax.plot(o2s[6:],Ts[6:],color='k',marker='.',linestyle='-',label='Temperature control')
ax.set(xscale='log',ylabel='Surface temperature (K)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row1[1])
ax.plot(o2s[:6],eqTs[:6]-plTs[:6],color='k',marker='.',linestyle='--',label='Constant CO$_2$',lw=2)
ax.plot(o2s[6:],eqTs[6:]-plTs[6:],color='k',marker='.',linestyle='-',label='Temperature control')
ax.set(xscale='log',ylabel='Equator-pole\n surface $\Delta T$ (K)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row1[2])
ax.plot(o2s[6:],o3col[6:],color='k',marker='.',linestyle='-',label='Temperature\n control')
ax.plot(o2s[:6],o3col[:6],color='k',marker='.',linestyle='--',label='Constant CO$_2$',lw=2)
ax.set(xscale='log',yscale='log',ylabel='Ozone column\n(kg m$^{-2}$)')
ax.tick_params(direction='in')
ax.legend(loc='lower right',fontsize=5,handlelength=4.5,frameon=False)
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row1[3])
ax.plot(o2s[:6],Ttroptrop[:6],color='k',marker='.',linestyle='--',label='Constant CO$_2$',lw=2)
ax.plot(o2s[6:],Ttroptrop[6:],color='k',marker='.',linestyle='-',label='Temperature control')
ax.set(xscale='log',ylabel='Tropical temperature\nat 100 hPa (K)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

#2nd row
ax = fig.add_subplot(row2[0])
ax.plot(o2s[:6],T20[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],T20[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='Temperature\nat 20 hPa (K)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row2[1])
ax.plot(o2s[:6],q20[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],q20[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',yscale='log',ylabel='Specfic Humidity \nat 20 hPa (kg kg$^{-1}$)')
ax.tick_params(direction='in',which='major')
ax.tick_params(direction='in',which='minor')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row2[2])
ax.plot(o2s[:6],q500[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],q500[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='Specific Humidity \nat 500 hPa (kg kg$^{-1}$)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row2[3])
ax.plot(o2s[:6],rh500[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],rh500[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='Relative humidity\nat 500 hPa (%)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

#3rd row
#determine optimal dynamic ranges for cloud properties
#this is just so that like quantities have the same dynamic range
dy_crf = np.max((np.ceil(np.abs(np.max(lwcf)-np.min(lwcf))),np.ceil(np.abs(np.max(swcf)-np.min(swcf)))))
dy_frc = 1.2*np.max((np.abs(np.max(cldhgh)-np.min(cldhgh)),np.abs(np.max(cldlow)-np.min(cldlow))))
dy_cwp = 1.2*np.max((np.abs(np.max(icwp)-np.min(icwp)),np.abs(np.max(lcwp)-np.min(lcwp))))

ax = fig.add_subplot(row3[0])
ax.plot(o2s[:6],cldhgh[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],cldhgh[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='High cloud fraction')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_frc),0.5*(ylims[0]+ylims[1]+dy_frc)))

ax = fig.add_subplot(row3[1])
ax.plot(o2s[:6],icwp[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],icwp[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='Ice CWP (kg m$^{-2}$)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_cwp),0.5*(ylims[0]+ylims[1]+dy_cwp)))

ax = fig.add_subplot(row3[2])
ax.plot(o2s[:6],icldt[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],icldt[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='Ice cloud top\ntemperature (K)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row3[3])
ax.plot(o2s[:6],icldp[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],icldp[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='Ice cloud top\nPressure (hPa)')
ax.tick_params(direction='in')
ax.invert_yaxis()
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

#4th row
ax = fig.add_subplot(row4[0])
ax.plot(o2s[:6],cldlow[:6],marker='.',linestyle='--',color='k',lw=2)
ax.plot(o2s[6:],cldlow[6:],marker='.',linestyle='-',color='k')
ax.set(xscale='log',ylabel='Low cloud fraction')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_frc),0.5*(ylims[0]+ylims[1]+dy_frc)))

ax = fig.add_subplot(row4[1])
ax.plot(o2s[:6],lcwp[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],lcwp[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='Liquid CWP (kg m$^{-2}$)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_cwp),0.5*(ylims[0]+ylims[1]+dy_cwp)))

ax = fig.add_subplot(row4[2])
ax.plot(o2s[:6],lts[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],lts[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='LTS (K)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row4[3])
ax.plot(o2s[:6],eis[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],eis[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='EIS (K)',xlabel='Nominal O$_2$ level (mol mol$^{-1}$)')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ax.tick_params(direction='in')

#5th row
ax = fig.add_subplot(row5[0])
ax.plot(o2s[:6],lwcf[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],lwcf[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='LW CRF (W m$^{-2}$)',xlabel='Nominal O$_2$ level (mol mol$^{-1}$)')
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_crf),0.5*(ylims[0]+ylims[1]+dy_crf)))
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])

ax = fig.add_subplot(row5[1])
ax.plot(o2s[:6],swcf[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],swcf[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='SW CRF (W m$^{-2}$)',xlabel='Nominal O$_2$ level (mol mol$^{-1}$)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_crf),0.5*(ylims[0]+ylims[1]+dy_crf)))

ax = fig.add_subplot(row5[2])
ax.plot(o2s[:6],lwcf[:6] + swcf[:6],color='k',marker='.',linestyle='--',lw=2)
ax.plot(o2s[6:],lwcf[6:] + swcf[6:],color='k',marker='.',linestyle='-')
ax.set(xscale='log',ylabel='LW + SW CRF\n(W m$^{-2}$)',xlabel='Nominal O$_2$ level (mol mol$^{-1}$)')
ax.tick_params(direction='in')
ax.set_xticks([1e-8,1e-6,1e-4,1e-2])
ylims = ax.get_ylim()
ax.set_ylim((0.5*(ylims[0]+ylims[1]-dy_crf),0.5*(ylims[0]+ylims[1]+dy_crf)))

plt.savefig('figures/fig01.pdf')
plt.close()
