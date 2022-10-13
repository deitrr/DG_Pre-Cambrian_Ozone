import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec

#this just suppresses a deprecation warning from netCDF4
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

simnames = ["ref_const_O3",
            "const_CO2_O2_1e-2",
            "temp_cont_O2_1e-3",
            "temp_cont_O2_1e-4",
            "temp_cont_O2_1e-7",
            "temp_cont_O2_1e-9"]

o2_title = ['Nominal O$_2=$ PAL',
            'Nominal O$_2=10^{-2}$',
            'Nominal O$_2=10^{-3}$',
            'Nominal O$_2=10^{-4}$',
            'Nominal O$_2=10^{-7}$',
            'Nominal O$_2=10^{-9}$']

fields = ['LWCF','CLDHGH','SWCF','CLDLOW']
cname = ['LW cloud forcing\n(W m$^{-2}$)','High cloud fraction','SW cloud forcing\n(W m$^{-2}$)','Low cloud fraction']
cmaps = ['viridis','viridis','viridis_r','viridis']
cranges = [[-2,94],[0.04,0.81],[-154,-0.1],[0,0.92]]

fig = plt.figure(figsize=(7.5,9))
outer_grid = gridspec.GridSpec(2,1,wspace=0.2,hspace=0.1,left=0.05,right=0.95,bottom=0.02,top=0.93,height_ratios=(1,5))
top_grid = gridspec.GridSpecFromSubplotSpec(2,4,subplot_spec=outer_grid[0],wspace=0.1,hspace=0.1,height_ratios=(1,15))
bot_grid = gridspec.GridSpecFromSubplotSpec(6,4,subplot_spec=outer_grid[1],wspace=0.1,hspace=0.1,height_ratios=(1,15,15,15,15,15))

for i in np.arange(len(simnames)):
  simpath = "CAM_post_process/" + simnames[i]
  file = simpath + "/" + simnames[i] + ".cam.h0.timmean_0031_0060.nc"

  data = nc.Dataset(file,'r')
  p = data.variables['lev'][:]
  lats = data.variables['lat'][:]
  lons = data.variables['lon'][:]
  lon2d, lat2d = np.meshgrid(lons,lats)
  if i == 0:
    pal_fields = {}

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

    #Map using the Eckert IV projection which is equal area so poles are not overemphasized
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


    c = m.pcolormesh(lon2d,lat2d,np.squeeze(field),cmap=cmaps[j],latlon='True',rasterized=True,vmin=cranges[j][0],vmax=cranges[j][-1])

    if j == 3:
      ax.text(1.03,0.5,o2_title[i],verticalalignment='center',rotation=270,transform=ax.transAxes,fontsize=8)

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

plt.savefig('figures/fig08.pdf')
plt.close()

