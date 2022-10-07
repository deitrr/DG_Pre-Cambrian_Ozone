import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pdb
import pathlib
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#simnames =["O3_during_goe_test","o3_after_goe_test","o3_proto2_test","o3_during_goe_hiCO2"]
simnames = ["o3_pal_avg_test","o3_proto2_test","o3_proto1_hiCO2","o3_after_goe_hiCO2","o3_during_goe_hiCO2","o3_pre_goe_hiCO2"]
#o2_title = ['Pre-GOE \n(O$_2$ $\\approx$ 5$\\times10^{-9}$)','Start of GOE\n(O$_2$ $\\approx$ 6$\\times10^{-8}$)','After GOE\n(O$_2$ $\\approx$ 7$\\times10^{-5}$)','PAL \n(O$_2$ $\\approx$ 0.21)']
#o2_title = ['Pre-GOE 2', 'Post-GOE 1','Post-GOE 3', 'Pre-GOE 2 (high CO$_2$)']
o2_title = ['Nominal O$_2=$ PAL','Nominal O$_2=10^{-2}$','Nominal O$_2=10^{-3}$','Nominal O$_2=10^{-4}$','Nominal O$_2=10^{-7}$','Nominal O$_2=10^{-9}$']
fields = ['LWCF','CLDHGH','SWCF','CLDLOW']
cname = ['LW cloud forcing\n(W m$^{-2}$)','High cloud fraction','SW cloud forcing\n(W m$^{-2}$)','Low cloud fraction']
cmaps = ['viridis','viridis','viridis_r','viridis']
cranges = [[-2,94],[0.04,0.81],[-154,-0.1],[0,0.92]]
#clabel_fmt = ['%d','%d','%d','%d']
#caxticks = [[150,190,230,270,310],[-9,-7,-5,-3,-1],
#fig, axes = plt.subplots(ncols=4,nrows=6,figsize=(7.5,9.0))
fig = plt.figure(figsize=(7.5,9))

outer_grid = gridspec.GridSpec(2,1,wspace=0.2,hspace=0.1,left=0.05,right=0.95,bottom=0.02,top=0.93,height_ratios=(1,5))
top_grid = gridspec.GridSpecFromSubplotSpec(2,4,subplot_spec=outer_grid[0],wspace=0.1,hspace=0.1,height_ratios=(1,15))
bot_grid = gridspec.GridSpecFromSubplotSpec(6,4,subplot_spec=outer_grid[1],wspace=0.1,hspace=0.1,height_ratios=(1,15,15,15,15,15))

for i in np.arange(len(simnames)):
  simpath = simnames[i] + "/merged_hist"
  file = simpath + "/" + simnames[i] + ".cam.h0.timmean_0031_0060.nc"

  data = nc.Dataset(file,'r')
  p = data.variables['lev'][:]
  lats = data.variables['lat'][:]
  lons = data.variables['lon'][:]
  lon2d, lat2d = np.meshgrid(lons,lats)
  if i == 0:
    pal_fields = {}

  #axes[0][i].set_title(o2_title[i],fontsize=8)

  for j in np.arange(len(fields)):
    if fields[j] == 'PT':
      field = data.variables['T'][:]*(1000/p[None,:,None,None])**0.286
    elif fields[j] == 'Q':
      field = np.log10(data.variables[fields[j]][:])
    elif fields[j] == 'OMEGA':
      field = data.variables[fields[j]][:]*100
    else:
      field = data.variables[fields[j]][:]

    if i == 0:
      pal_fields[fields[j]] = field
    else:
      field -= pal_fields[fields[j]]
      cmaps = ['RdBu_r','RdBu_r','RdBu_r','RdBu_r']
      cranges = [[-24,24],[-0.19,0.19],[-35,35],[-0.24,0.24]]

    if i == 0:
      ax = fig.add_subplot(top_grid[j+4])
    else:
      ax = fig.add_subplot(bot_grid[4*i+j])

    m = Basemap(lat_0=0,lon_0=0,ax=ax,fix_aspect=False,projection='eck4',resolution='c')
    m.drawcoastlines(linewidth=0.5)
    if j == 0:
      m.drawparallels([-45,0,45],labels = [True,False,False,False], fontsize=6)
    else:
      m.drawparallels([-45,0,45],labels = [False,False,False,False], fontsize=6)
    if i == 0 or i == 5:
      m.drawmeridians([-90,0,90],labels = [False,False,False,True], fontsize=6)
    else:
      m.drawmeridians([-90,0,90],labels = [False,False,False,False], fontsize=6)


    print(fields[j],np.min(field),np.max(field))
    #ax.invert_yaxis()
    #if j < 3:
    #ax.set_yscale('log')

    c = m.pcolormesh(lon2d,lat2d,np.squeeze(field),cmap=cmaps[j],latlon='True',rasterized=True,vmin=cranges[j][0],vmax=cranges[j][-1])

    #cont = ax.contour(lats, p, np.squeeze(field), cranges[j][::2], colors='k',linewidths=1)
    #ax.clabel(cont, inline =True,fontsize=6,fmt=clabel_fmt[j])
   # if j == 0:
    #  ax.set(ylabel='Pressure (hPa)')
    #else:
    #  ax.set_yticklabels([])

    if j == 3:
      ax.text(1.03,0.5,o2_title[i],verticalalignment='center',rotation=270,transform=ax.transAxes,fontsize=8)

    #if i == len(simnames)-1:
    #  ax.set(xlabel='Latitude (deg)')
    #else:
    #  ax.set_xticklabels([])

    if i == 0 or i == 1:
      if i == 0:
        cax = fig.add_subplot(top_grid[j])
      elif i == 1:
        cax = fig.add_subplot(bot_grid[j])
      #cax = make_axes_locatable(axes[i][j]).append_axes('top',size='5%',pad=0.05)
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cax.tick_params(axis='x',labelsize=6)
      if i == 0:
        cbar.set_label(cname[j],fontsize=8)


  #data.close()

#`plt.tight_layout(h_pad=0.1,w_pad=0.3)
plt.savefig('re_horiz_diff_cld.pdf')
plt.close()

