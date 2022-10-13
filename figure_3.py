import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import pdb

#suppress deprecation warning from netCDF4
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#simnames =["o3_pre_goe_hiCO2","o3_during_goe_hiCO2",
#           "o3_after_goe_hiCO2","o3_proto1_hiCO2","o3_proto2_test",
#           "o3_pre_goe_test","O3_during_goe_test","o3_after_goe_test","o3_proto1_test",
#           "o3_pal_avg_test","bps_control_1.0"]
simnames = ["temp_cont_O2_1e-9",
            "temp_cont_O2_1e-7",
            "temp_cont_O2_1e-4",
            "temp_cont_O2_1e-3",
            "const_CO2_O2_1e-2",
            "const_CO2_O2_1e-9",
            "const_CO2_O2_1e-7",
            "const_CO2_O2_1e-4",
            "const_CO2_O2_1e-3",
            "ref_const_O3",
            "ref_vary_O3"]

#this is used for the legend
o2_title = ['Nominal O$_2=10^{-9}$',
            'Nominal O$_2=10^{-7}$',
            'Nominal O$_2=10^{-4}$',
            'Nominal O$_2=10^{-3}$',
            'Nominal O$_2=10^{-2}$',
            '',
            '',
            '',
            '',
            'Nominal O$_2$ = PAL,\n  constant O$_3$',
            'Nominal O$_2$ = PAL,\n  varying O$_3$']

#fields = ['T']
#cname = ['Temperature']
cmap = plt.cm.seismic_r
colors=[cmap(0.8),cmap(0.6),cmap(0.4),cmap(0.2),cmap(0),cmap(0.8),cmap(0.6),cmap(0.4),cmap(0.2),'k','0.5']
lstyles=['-','-','-','-','-','--','--','--','--','-','-']
lwidths=[3,1.5,2,2,2,3,1.5,2,2,2,1]

#set up figure, axes, insets
fig, axes = plt.subplots(ncols=2,nrows=3,figsize=(7.5,8))
inset00 = axes[0][0].inset_axes([0.65,0.5,0.34,0.47])
inset10 = axes[1][0].inset_axes([0.63,0.5,0.34,0.47])
inset20 = axes[2][0].inset_axes([0.65,0.5,0.34,0.47])

for i in np.arange(len(simnames)):
  print(simnames[i])

  simpath = "CAM_post_process/" + simnames[i]
  file = simpath + "/" + simnames[i] + ".cam.h0.globmean_0031_0060.nc"
  data = nc.Dataset(file,'r')

  lats = data.variables['lat'][:]
  T = np.squeeze(data.variables['T'][:])
  Q = np.squeeze(data.variables['Q'][:])
  rh = np.squeeze(data.variables['RELHUM'][:])

  #compute pressure levels using average surface pressure
  data_ic = nc.Dataset("cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc")
  p = (data_ic['hyam'][:]*data_ic['P0'][:] + data_ic['hybm'][:]*data['PS'][:]).squeeze()/100

  #if simnames[i] == "o3_pal_avg_test":
  #  Tref = T
  #  Qref = Q
  #  rhref = rh
  #if simnames[i] == "bps_control_1.0":
  #  print("dT = (%f,%f)"%(np.max(np.abs((T-Tref))),np.mean((T-Tref))))
  #  print("dQ = (%f,%f)"%(np.max(np.abs((np.log10(Q)-np.log10(Qref)))),np.mean((np.log10(Q)-np.log10(Qref)))))
  #  print("drh = (%f,%f)"%(np.max(np.abs(rh-rhref)),np.mean(rh-rhref)))

  axes[0][0].plot(T,p,c=colors[i],lw=lwidths[i],linestyle=lstyles[i])
  inset00.semilogy(T,p,c=colors[i],lw=lwidths[i]/2.,linestyle=lstyles[i])

  axes[1][0].plot(Q,p,c=colors[i],lw=lwidths[i],linestyle=lstyles[i])
  inset10.semilogy(Q,p,c=colors[i],lw=lwidths[i]/2.,linestyle=lstyles[i])

  axes[2][0].plot(rh,p,c=colors[i],lw=lwidths[i],linestyle=lstyles[i])
  inset20.semilogy(rh,p,c=colors[i],lw=lwidths[i]/2.,linestyle=lstyles[i])

  #data.close()

  #file = simpath + "/" + simnames[i] + ".cam.h0.tropmean_0031_0060.nc"

  #data = nc.Dataset(file,'r')
  #p = data.variables['lev'][:]
  #lats = data.variables['lat'][:]
  cloud = np.squeeze(data.variables['CLOUD'][:])
  iwc = np.squeeze(data.variables['IWC'][:])
  lwc = np.squeeze(data['LWC'][:])

  #if simnames[i] == "o3_pal_avg_test":
  #  CLref = T
  #  IWCref = Q
  #  LWCref = rh
  #if simnames[i] == "bps_control_1.0":
  #  print("dCF = (%f,%f)"%(np.max(np.abs((T-CLref))),np.mean((T-CLref))))
  #  print("dIWC = (%f,%f)"%(np.max(np.abs((np.log10(Q)-np.log10(IWCref)))),np.mean((np.log10(Q)-np.log10(IWCref)))))
  #  print("dLWC = (%f,%f)"%(np.max(np.abs((np.log10(rh)-np.log10(LWCref)))),np.mean((np.log10(rh)-np.log10(LWCref)))))
  #  pdb.set_trace()

  axes[0][1].plot(cloud,p,c=colors[i],lw=lwidths[i],label=o2_title[i],linestyle=lstyles[i])

  axes[1][1].plot(iwc/1e-6,p,c=colors[i],lw=lwidths[i],linestyle=lstyles[i])

  if o2_title[i] != '':
    axes[2][1].plot(lwc/1e-6,p,c=colors[i],lw=lwidths[i],label=o2_title[i],linestyle=lstyles[i])
  else:
    axes[2][1].plot(lwc/1e-6,p,c=colors[i],lw=lwidths[i],label=o2_title[i],linestyle=lstyles[i])

  data.close()

#axis labels, ranges, etc.
axes[0][0].set_yscale('log')
axes[0][0].invert_yaxis()
axes[0][0].set(xlabel='Temperature (K)',ylabel='Pressure (hPa)')
inset00.set(xlim=(250,290),ylim=(500,1050))
inset00.invert_yaxis()
inset00.set_yticklabels([1],minor=True)

axes[0][1].invert_yaxis()
axes[0][1].set_yscale('log')
axes[0][1].set(xlabel='Cloud fraction',ylim=(1050,50))

axes[1][0].set_yscale('log')
axes[1][0].set_xscale('log')
axes[1][0].invert_yaxis()
axes[1][0].set(xlabel='Specific Humidity (kg kg$^{-1}$)',ylabel='Pressure (hPa)')
inset10.set(xlim=(1e-3,1.2e-2),ylim=(500,1050))
inset10.invert_yaxis()
inset10.set_xscale('log')
inset10.set_yticklabels([1],minor=True)

axes[1][1].invert_yaxis()
axes[1][1].set_yscale('log')
axes[1][1].set(xlabel='Ice cloud density (10$^{-6}$ kg m$^{-3}$)',ylim=(1050,50))

axes[2][0].set_yscale('log')
axes[2][0].invert_yaxis()
axes[2][0].set(xlabel='Relative Humidity (%)',ylabel='Pressure (hPa)')
inset20.set(xlim=(45,85),ylim=(500,1050))
inset20.invert_yaxis()
inset20.set_yticklabels([1],minor=True)

axes[2][1].invert_yaxis()
axes[2][1].set_yscale('log')
axes[2][1].set(xlabel='Liquid water density (10$^{-6}$ kg m$^{-3}$)',ylim=(1050,450))
axes[2][1].legend(loc='upper right',fontsize=6,handlelength=3)

plt.tight_layout()
plt.savefig('figures/fig03.pdf')
plt.close()


