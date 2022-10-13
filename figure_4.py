import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
#from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import pdb
#import pathlib

#simnames = ["o3_pal_avg_test","o3_proto2_test","o3_proto1_hiCO2","o3_after_goe_hiCO2","o3_during_goe_hiCO2","o3_pre_goe_hiCO2"]
simnames = ["ref_const_O3",
            "const_CO2_O2_1e-2",
            "temp_cont_O2_1e-3",
            "temp_cont_O2_1e-4",
            "temp_cont_O2_1e-7",
            "temp_cont_O2_1e-9"]

#right side labels
o2_title = ['Nominal O$_2=$ PAL',
            'Nominal O$_2=10^{-2}$',
            'Nominal O$_2=10^{-3}$',
            'Nominal O$_2=10^{-4}$',
            'Nominal O$_2=10^{-7}$',
            'Nominal O$_2=10^{-9}$']

fields = ['T','Q','U','RELHUM']
cname = ['Temperature (K)','log(Humidity [kg kg$^{-1}$])','Zonal velocity (m s$^{-1}$)','Relative Humidity (%)']
cmaps = ['viridis_r','viridis_r','RdBu_r','viridis_r']

#ranges for top row
cranges = [np.arange(180,310,10),np.linspace(np.log10(1e-6),np.log10(0.1),11),np.arange(-55,60,5),np.arange(0,135,10)]

#label format for contours
clabel_fmt = ['%d','%d','%d','%d']

fig = plt.figure(figsize=(7.5,9))

outer_grid = gridspec.GridSpec(2,1,wspace=0.2,hspace=0.1,left=0.1,right=0.95,bottom=0.05,top=0.95,height_ratios=(1,5))
top_grid = gridspec.GridSpecFromSubplotSpec(2,4,subplot_spec=outer_grid[0],wspace=0.1,hspace=0.1,height_ratios=(1,15))
bot_grid = gridspec.GridSpecFromSubplotSpec(6,4,subplot_spec=outer_grid[1],wspace=0.1,hspace=0.1,height_ratios=(1,15,15,15,15,15))

for i in np.arange(len(simnames)):
  simpath = "CAM_post_process/" + simnames[i]
  file = simpath + "/" + simnames[i] + ".cam.h0.zonmean_0031_0060.nc"

  data = nc.Dataset(file,'r')
  p = data['lev'][:]
  lats = data['lat'][:]

  #dictionary to store reference values for differences
  if i == 0:
    pal_fields = {}

  for j in np.arange(len(fields)):
    if fields[j] == 'Q':
      field = np.log10(data[fields[j]][:])
    else:
      field = data[fields[j]][:]

    if i == 0:
      pal_fields[fields[j]] = field
    else:
      field -= pal_fields[fields[j]]
      #reset colormaps and ranges for difference plots
      cmaps = ['RdBu_r','RdBu_r','RdBu_r','RdBu_r']
      cranges = [np.arange(-80,90,10),np.linspace(-3,3,13),np.arange(-40,41,5),np.arange(-80,81,10)]

    if i == 0:
      ax = fig.add_subplot(top_grid[j+4])
    else:
      ax = fig.add_subplot(bot_grid[4*i+j])

    c = ax.contourf(lats, p, np.squeeze(field), cranges[j], cmap=cmaps[j], rasterized=True) 
    for cc in c.collections:
      cc.set_edgecolor("face")
#    print(fields[j],np.min(field),np.max(field))
    ax.invert_yaxis()
    ax.set_yscale('log')

    cont = ax.contour(lats, p, np.squeeze(field), cranges[j][::2], colors='k',linewidths=1)
    ax.clabel(cont, inline =True,fontsize=6,fmt=clabel_fmt[j])
    if j == 0:
      ax.set(ylabel='Pressure (hPa)')
    else:
      ax.set_yticklabels([])

    ax.set_ylim((1e3,3.5))
    if j == 3:
      ax.text(1.03,0.5,o2_title[i],verticalalignment='center',rotation=270,transform=ax.transAxes,fontsize=8)

    ax.set_xticks([-90,-45,0,45,90])
    if i == len(simnames)-1:
      ax.set(xlabel='Latitude (deg)')
      ax.set_xticklabels(['  -90','-45','0','45','90  '],fontsize=8)
    else:
      ax.set_xticklabels([])

    if i == 0 or i == 1:
      if i == 0:
        cax = fig.add_subplot(top_grid[j])
      elif i == 1:
        cax = fig.add_subplot(bot_grid[j])
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cax.tick_params(axis='x',labelsize=6)
      if i == 0:
        cbar.set_label(cname[j],fontsize=8)


plt.savefig('figures/fig04.pdf')
plt.close()

