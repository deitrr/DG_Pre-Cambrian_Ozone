import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import pdb
import pathlib
import scipy.integrate as intg
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#simnames =["O3_during_goe_test","o3_after_goe_test","o3_proto2_test","o3_during_goe_hiCO2"]
simnames = ["o3_pal_avg_test","o3_proto2_test","o3_proto1_hiCO2","o3_after_goe_hiCO2","o3_during_goe_hiCO2","o3_pre_goe_hiCO2"]
#o2_title = ['Pre-GOE \n(O$_2$ $\\approx$ 5$\\times10^{-9}$)','Start of GOE\n(O$_2$ $\\approx$ 6$\\times10^{-8}$)','After GOE\n(O$_2$ $\\approx$ 7$\\times10^{-5}$)','PAL \n(O$_2$ $\\approx$ 0.21)']
#o2_title = ['Pre-GOE 2', 'Post-GOE 1','Post-GOE 3', 'Pre-GOE 2 (high CO$_2$)']
#o2_title = ['PAL O$_2$','$\sim10^{-2}$ O$_2$','$\sim10^{-3}$ O$_2$','$\sim10^{-4}$ O$_2$','$\sim10^{-6}$ O$_2$','$\sim10^{-8}$ O$_2$']
o2_title = ['Nominal O$_2=$ PAL','Nominal O$_2=10^{-2}$','Nominal O$_2=10^{-3}$','Nominal O$_2=10^{-4}$','Nominal O$_2=10^{-7}$','Nominal O$_2=10^{-9}$']

#fields = ['OMEGA','CLOUD','IWC','LWC']
#cname = ['Vertical velocity\n(10$^{-2}$ Pa s$^{-1}$)','Cloud fraction','log(Ice cloud density\n[kg m$^{-3}$])','log(Liquid cloud density\n[kg m$^{-3}$])']
fields = ['streamf','CLOUD','IWC','LWC']
cname = ['Mass stream function\n(10$^{10}$ kg s$^{-1}$)','Cloud fraction','log(Ice cloud density\n[kg m$^{-3}$])','log(Liquid cloud density\n[kg m$^{-3}$])']
cmaps = ['RdBu_r','viridis_r','viridis_r','viridis_r']
#cranges = [np.arange(-5,5.5,0.5),np.linspace(0,0.55,12),np.arange(-15,-3,1),np.arange(-15,-3,1)]
cranges = [np.arange(-10,12,2),np.linspace(0,0.55,12),np.arange(-12,-3,1),np.arange(-12,-3,1)]
#cranges = [np.arange(-5,5.5,0.5),np.linspace(0,0.55,12),np.arange(-12,-1,1),np.arange(-12,0,1)]
clabel_fmt = ['%d','%#.2f','%d','%d']
caxticks = [[-4,-3,-2,-1,0,1,2,3,4],[-0.15,-0.1,-0.05,0,0.05,0.1,0.15]]
#fig, axes = plt.subplots(ncols=4,nrows=6,figsize=(7.5,9.0))
fig = plt.figure(figsize=(7.5,9))

outer_grid = gridspec.GridSpec(2,1,wspace=0.2,hspace=0.1,left=0.1,right=0.95,bottom=0.05,top=0.93,height_ratios=(1,5))
top_grid = gridspec.GridSpecFromSubplotSpec(2,4,subplot_spec=outer_grid[0],wspace=0.1,hspace=0.1,height_ratios=(1,15))
bot_grid = gridspec.GridSpecFromSubplotSpec(6,4,subplot_spec=outer_grid[1],wspace=0.1,hspace=0.1,height_ratios=(1,15,15,15,15,15))

R = 287.0
g = 9.80665
a = 6371000.

for i in np.arange(len(simnames)):
  simpath = simnames[i] + "/merged_hist"
  file = simpath + "/" + simnames[i] + ".cam.h0.zonmean_0031_0060.nc"

  data = nc.Dataset(file,'r')
  p = data.variables['lev'][:]
  lats = data.variables['lat'][:]
  if i == 0:
    pal_fields = {}

  file_tim = simpath + "/" + simnames[i] + ".cam.h0.timmean_0031_0060.nc"
  data_tim = nc.Dataset(file_tim,'r')

  p_int = data_tim.variables['ilev'][:].squeeze()*100.0
  T = data.variables['T'][:].squeeze()
  rho = p[:,None]*100 / (R * T)
  dz = np.diff(p_int)[:,None] / (rho*g)

  for j in np.arange(len(fields)):
    if fields[j] == 'PT':
      field = data.variables['T'][:]*(1000/p[None,:,None,None])**0.286
    elif fields[j] == 'Q':
      field = np.log10(data.variables[fields[j]][:])
    elif fields[j] == 'OMEGA':
      field = data.variables[fields[j]][:]*100
    elif fields[j] == 'streamf':
      v = data.variables['V'][:].squeeze()
      field = np.zeros_like(v)
      for ilev in np.arange(1,len(p)):
        psi = intg.simps(v[:ilev+1,:],x=p[:ilev+1]*100,axis=0)
        field[ilev,:] = psi * 2 * np.pi * a / g * np.cos(lats[None,:]*np.pi/180)/1e10
    elif fields[j] == 'IWC' or fields[j] == 'LWC':
      field = np.log10(data.variables[fields[j]][:])
      field[field<-12] = np.nan
      #field[field==0] = 1e-30
      #field[np.isnan(field)] = 1e-30
      #field[field<-15] = np.nan
      #field = data.variables[fields[j]][:]*1e6
      #pdb.set_trace()
      #field = np.log10(data.variables[fields[j]][:].squeeze()*dz)
    else:
      field = data.variables[fields[j]][:]


    if j <= 1:
      if i == 0:
        pal_fields[fields[j]] = field
      else:
        field -= pal_fields[fields[j]]
        cmaps = ['RdBu_r','RdBu_r','viridis_r','viridis_r']
        cranges = [np.arange(-3.5,4,0.5),np.linspace(-0.16,0.16,17),np.arange(-12,-3,1),np.arange(-12,-3,1)]
    print(fields[j],np.nanmin(field),np.nanmax(field))

    if i == 0:
      ax = fig.add_subplot(top_grid[j+4])
    else:
      ax = fig.add_subplot(bot_grid[4*i+j])

    c = ax.contourf(lats, p, np.squeeze(field), cranges[j], cmap=cmaps[j]) 
    for cc in c.collections:
      cc.set_edgecolor("face")
    ax.invert_yaxis()
    #if j < 3:
    ax.set_yscale('log')
    ax.set(ylim=(1000,50))

    if j > 0:
      cont = ax.contour(lats, p, np.squeeze(field), cranges[j][::2], colors='k',linewidths=1)
      ax.clabel(cont, inline =True,fontsize=6,fmt=clabel_fmt[j])
    else:
      cont = ax.contour(lats, p, np.squeeze(field), cranges[j][1::2], colors='k',linewidths=1)
      ax.clabel(cont, inline =True,fontsize=6,fmt=clabel_fmt[j])
    if j == 0:
      ax.set(ylabel='Pressure (hPa)')
      #axes[i][j].set(xlim=(-50,50))
    else:
      ax.set_yticklabels([])

    if j == 3:
      ax.text(1.03,0.5,o2_title[i],verticalalignment='center',rotation=270,transform=ax.transAxes,fontsize=8)

    ax.set_xticks([-90,-45,0,45,90])
    if i == len(simnames)-1:
      ax.set(xlabel='Latitude (deg)')
      ax.set_xticklabels(['  -90','-45','0','45','90  '])
    else:
      ax.set_xticklabels([])

    if i == 0:
      cax = fig.add_subplot(top_grid[j])
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cbar.set_label(cname[j],fontsize=8)
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cax.tick_params(axis='x',labelsize=5)

    elif i == 1:
      if j <= 1:
        cax = fig.add_subplot(bot_grid[j]) 
        #cax = make_axes_locatable(axes[i][j]).append_axes('top',size='5%',pad=0.05)
        cbar = plt.colorbar(c,cax=cax,orientation='horizontal',ticks=caxticks[j])
        cax.xaxis.set_ticks_position('top')
        cax.xaxis.set_label_position('top')
        cax.tick_params(axis='x',labelsize=5)

  data.close()

#plt.tight_layout(h_pad=0.1,w_pad=0.3)
plt.savefig('re_select_fields_zonal_cld_diff_comb_tempcont.pdf')
plt.close()

