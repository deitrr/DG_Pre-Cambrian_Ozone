import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pdb
import pathlib

#simnames =["O3_during_goe_test","o3_after_goe_test","o3_proto2_test","o3_during_goe_hiCO2"]
simnames = ["o3_pal_avg_test","o3_proto2_test","o3_proto1_hiCO2","o3_after_goe_hiCO2","o3_during_goe_hiCO2","o3_pre_goe_hiCO2"]
#o2_title = ['Pre-GOE \n(O$_2$ $\\approx$ 5$\\times10^{-9}$)','Start of GOE\n(O$_2$ $\\approx$ 6$\\times10^{-8}$)','After GOE\n(O$_2$ $\\approx$ 7$\\times10^{-5}$)','PAL \n(O$_2$ $\\approx$ 0.21)']
#o2_title = ['Pre-GOE 2', 'Post-GOE 1','Post-GOE 3', 'Pre-GOE 2 (high CO$_2$)']
o2_title = ['Nominal O$_2=$ PAL','Nominal O$_2=10^{-2}$','Nominal O$_2=10^{-3}$','Nominal O$_2=10^{-4}$','Nominal O$_2=10^{-7}$','Nominal O$_2=10^{-9}$']
fields = ['T','U']
cname = ['Zonal velocity (m s$^{-1}$), GCM','Zonal velocity (m s$^{-1}$), thermal wind approximation']
cmaps = ['RdBu_r','RdBu_r']
cranges = [np.arange(-60,65,5),np.arange(-60,65,5)]
clabel_fmt = ['%d','%d','%d','%d']
#caxticks = [[150,190,230,270,310],[-9,-7,-5,-3,-1],
#fig, axes = plt.subplots(ncols=4,nrows=6,figsize=(7.5,9.0))
fig = plt.figure(figsize=(7.5,9))

outer_grid = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.1,left=0.1,right=0.95,bottom=0.05,top=0.95,height_ratios=(1,))
#top_grid = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=outer_grid[0],wspace=0.1,hspace=0.1,height_ratios=(1,15))
bot_grid = gridspec.GridSpecFromSubplotSpec(7,2,subplot_spec=outer_grid[0],wspace=0.1,hspace=0.1,height_ratios=(1,15,15,15,15,15,15))

R = 287.0
Omega = 7.292e-5
a = 6.371e6

#cam_ic = nc.Dataset('cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc','r')
#p_int = cam_ic['ilev'][:].squeeze()*100.0
#dlnp = np.log(p_int)[1:] - np.log(p_int)[:-1]

for i in np.arange(len(simnames)):
  simpath = simnames[i] + "/merged_hist"
  file = simpath + "/" + simnames[i] + ".cam.h0.zonmean_0031_0060.nc"

  data = nc.Dataset(file,'r')
  p = data.variables['lev'][:]
  dlnp = -(np.log(p*100)[1:] - np.log(p*100)[:-1])
  lats = data.variables['lat'][:]
  f = 2*Omega*np.sin(lats*np.pi/180)
  dy = a * (lats[1]-lats[0])*np.pi/180
  y = a*lats*np.pi/180
  if i == 0:
    pal_fields = {}

  u_ref = data['U'][:].squeeze()
  u = u_ref.copy()
  T = data['T'][:].squeeze()

  Tav = np.zeros_like(lats)
  for ilev in np.arange(11,-1,-1):
    #pdb.set_trace()
    Tav += 0.5*(T[ilev,:]+T[ilev+1,:])/12
    dTdy = np.gradient(0.5*(T[ilev,:]+T[ilev+1,:]),y)
    u[ilev,:] = u[ilev+1,:] + R/f * dTdy * dlnp[ilev]
    #print(dTdy[13],f[13],u[ilev+1,13])
  dTavdy = np.gradient(Tav,y)
  print(u_ref[12,13],dTavdy[13],u_ref[12,13]+R/f[13]*dTavdy[13]*np.log(p[0]/p[12]))

  #pdb.set_trace()

  for j in np.arange(2):
    if j == 0:
      field = u_ref
    else:
      field = u

    #if i == 0:
    #  ax = fig.add_subplot(top_grid[j+2])
    #else:
    ax = fig.add_subplot(bot_grid[2*i+j+2])

    c = ax.contourf(lats, p, np.squeeze(field), cranges[j], cmap=cmaps[j], rasterized=True) 
    for cc in c.collections:
      cc.set_edgecolor("face")
#    print(fields[j],np.min(field),np.max(field))
    ax.invert_yaxis()
    #if j < 3:
    ax.set_yscale('log')

    cont = ax.contour(lats, p, np.squeeze(field), cranges[j][::2], colors='k',linewidths=1)
    ax.clabel(cont, inline =True,fontsize=6,fmt=clabel_fmt[j])
    if j == 0:
      ax.set(ylabel='Pressure (hPa)')
    else:
      ax.set_yticklabels([])

    if j == 1:
      ax.text(1.03,0.5,o2_title[i],verticalalignment='center',rotation=270,transform=ax.transAxes,fontsize=8)

    if i == len(simnames)-1:
      ax.set(xlabel='Latitude (deg)')
    else:
      ax.set_xticklabels([])
    
    if j == 1:
      ax.fill_between([-90,90],[p[12],p[12]],[p[29],p[29]],fc='0.7',ec='0.5',alpha=0.5,hatch='x')

    #if i == 0 or i == 1:
    if i == 0:
        #cax = fig.add_subplot(top_grid[j])
      #elif i == 1:
      cax = fig.add_subplot(bot_grid[j])
      #cax = make_axes_locatable(axes[i][j]).append_axes('top',size='5%',pad=0.05)
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cax.tick_params(axis='x',labelsize=6)
      #if i == 0:
      cbar.set_label(cname[j],fontsize=8)


  #data.close()

#`plt.tight_layout(h_pad=0.1,w_pad=0.3)
plt.savefig('u_thermal.pdf')
plt.close()

