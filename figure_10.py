import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
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

cname = ['Zonal velocity (m s$^{-1}$), GCM','Zonal velocity (m s$^{-1}$), thermal wind approximation']
cmaps = ['RdBu_r','RdBu_r']
cranges = [np.arange(-60,65,5),np.arange(-60,65,5)]
clabel_fmt = ['%d','%d','%d','%d']

fig = plt.figure(figsize=(7.5,9))

outer_grid = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.1,left=0.1,right=0.95,bottom=0.05,top=0.95,height_ratios=(1,))
bot_grid = gridspec.GridSpecFromSubplotSpec(7,2,subplot_spec=outer_grid[0],wspace=0.1,hspace=0.1,height_ratios=(1,15,15,15,15,15,15))

R = 287.0
Omega = 7.292e-5
a = 6.371e6

for i in np.arange(len(simnames)):
  simpath = "CAM_post_process/" + simnames[i]
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
    #integrate the thermal wind equation from the tropopause upward
    Tav += 0.5*(T[ilev,:]+T[ilev+1,:])/12
    dTdy = np.gradient(0.5*(T[ilev,:]+T[ilev+1,:]),y)
    u[ilev,:] = u[ilev+1,:] + R/f * dTdy * dlnp[ilev]


  for j in np.arange(2):
    if j == 0:
      field = u_ref
    else:
      field = u

    ax = fig.add_subplot(bot_grid[2*i+j+2])

    c = ax.contourf(lats, p, np.squeeze(field), cranges[j], cmap=cmaps[j])
    for cc in c.collections:
      cc.set_edgecolor("face")
    ax.invert_yaxis()
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

    if i == 0:
      cax = fig.add_subplot(bot_grid[j])
      cbar = plt.colorbar(c,cax=cax,orientation='horizontal')
      cax.xaxis.set_ticks_position('top')
      cax.xaxis.set_label_position('top')
      cax.tick_params(axis='x',labelsize=6)
      cbar.set_label(cname[j],fontsize=8)


plt.savefig('figures/fig10.pdf')
plt.close()

