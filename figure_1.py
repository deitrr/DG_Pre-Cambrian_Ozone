import numpy as np
import scipy.interpolate as sint
import matplotlib.pyplot as plt
import netCDF4 as nc

#this just suppresses a deprecation warning from netCDF4
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

path_to_archive = "DG_Pre-Cambrian_O3_archive/CAM_input/"

#These are the flux levels I selected from Daniel's photochem data
bgoe_flux = 3.2e11
dgoe_flux = 3.55e11
agoe_flux = 3.85e11
agoe_flux_1 = 4.2e11
agoe_flux_2 = 5.05e11

#Read in Daniel's photochem output (just for O2 really)
z, p, T, RH, f_O2, O2, O3 = np.loadtxt(path_to_archive +
                        "increase_O2_flux_fratio_0.3_O2_O3_mixrat_profiles.csv",
                        unpack=True,skiprows=1,delimiter=',',
                        usecols=(0,1,2,4,5,6,7))

z = np.unique(z)
T = np.unique(T)
f_O2 = np.unique(f_O2)

p = np.reshape(p,(len(f_O2),len(z)))
O2 = np.reshape(O2,(len(f_O2),len(z)))
O3 = np.reshape(O3,(len(f_O2),len(z)))

b_goe = np.where(f_O2==bgoe_flux)[0][0]
d_goe = np.where(f_O2==dgoe_flux)[0][0]
a_goe = np.where(f_O2==agoe_flux)[0][0]
pal = -81

a_goe1 = np.where(f_O2==agoe_flux_1)[0][0]
a_goe2 = np.where(f_O2==agoe_flux_2)[0][0]

#Read in CAM input files for a few necessary things
cam_o3 = nc.Dataset(path_to_archive+"ozone_1.9x2.5_L26_1850clim_c090420.nc")
cam_ic = nc.Dataset(path_to_archive+"cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc")

#convert hybrid levels to pressure for cam ic file
cam_ic_p = cam_ic['hyam'][:][:,None,None]*cam_ic['P0'][:] + \
            cam_ic['PS'][0][None,:,:,]*cam_ic['hybm'][:][:,None,None]

#interpolate to levels in o3 file
cam_o3_p = sint.interp1d(cam_ic['lev'][:].data,cam_ic_p.data,axis=0,
                         fill_value='extrapolate')(cam_o3['lev'][:].data)

#plot only the data in the pressure range used in CAM
plot_index = (p>=np.min(cam_o3_p/1e5))
fig, (ax1, ax2) = plt.subplots(ncols=2,nrows=1,figsize=(7.5,3))

cmap = plt.cm.seismic_r
colors = [cmap(0.8),cmap(0.6),cmap(0.4),cmap(0.2),cmap(0.),'k']
ax1.plot(O2[pal,plot_index[pal]],p[pal,plot_index[pal]],color=colors[5],label='Nom O$_2$ = PAL')
ax1.plot(O2[a_goe2,plot_index[a_goe2]],p[a_goe2,plot_index[a_goe2]],color=colors[4],label='Nom O$_2$ = $\sim10^{-2}$')
ax1.plot(O2[a_goe1,plot_index[a_goe1]],p[a_goe1,plot_index[a_goe2]],color=colors[3],label='Nom O$_2$ = $\sim10^{-3}$')
ax1.plot(O2[a_goe,plot_index[a_goe]],p[a_goe,plot_index[a_goe]],color=colors[2],label='Nom O$_2$ = $\sim10^{-4}$')
ax1.plot(O2[d_goe,plot_index[d_goe]],p[d_goe,plot_index[d_goe]],color=colors[1],label='Nom O$_2$ = $\sim10^{-7}$')
ax1.plot(O2[b_goe,plot_index[b_goe]],p[b_goe,plot_index[b_goe]],color=colors[0],label='Nom O$_2$ = $\sim10^{-9}$')
ax1.set(yscale='log',xscale='log',xlabel='O$_2$ mixing ratio (mol mol$^{-1}$)',ylabel='Pressure (bar)')
ax1.invert_yaxis()

ax2.set(yscale='log',xscale='log',xlabel='O$_3$ mixing ratio (mol mol$^{-1}$)')
ax2.invert_yaxis()

#average the earth varying profiles over the year
cam_avg_O3 = np.zeros_like(cam_o3['O3'][:])
for i in np.arange(12):
    cam_avg_O3[i] = np.mean(cam_o3['O3'],axis=0)

#interpolation for new O3 profiles
new_o3_bgoe = sint.interp1d(p[b_goe,:],O3[b_goe,:],
                            fill_value='extrapolate')(cam_o3_p/1e5)
new_o3_dgoe = sint.interp1d(p[d_goe,:],O3[d_goe,:],
                            fill_value='extrapolate')(cam_o3_p/1e5)
new_o3_agoe = sint.interp1d(p[a_goe,:],O3[a_goe,:],
                            fill_value='extrapolate')(cam_o3_p/1e5)
new_o3_pal = sint.interp1d(p[pal,:],O3[pal,:],
                            fill_value='extrapolate')(cam_o3_p/1e5)
new_o3_agoe1 = sint.interp1d(p[a_goe1,:],O3[a_goe1,:],
                            fill_value='extrapolate')(cam_o3_p/1e5)
new_o3_agoe2 = sint.interp1d(p[a_goe2,:],O3[a_goe2,:],
                            fill_value='extrapolate')(cam_o3_p/1e5)

ax2.plot(new_o3_pal[:,48,0],cam_o3_p[:,48,0]/1e5,c=colors[5],label='Nominal O$_2$ = PAL')
ax2.plot(new_o3_agoe2[:,48,0],cam_o3_p[:,48,0]/1e5,c=colors[4],label='Nominal O$_2$ = $10^{-2}$')
ax2.plot(new_o3_agoe1[:,48,0],cam_o3_p[:,48,0]/1e5,c=colors[3],label='Nominal O$_2$ = $10^{-3}$')
ax2.plot(new_o3_agoe[:,48,0],cam_o3_p[:,48,0]/1e5,c=colors[2],label='Nominal O$_2$ = $10^{-4}$')
ax2.plot(new_o3_dgoe[:,48,0],cam_o3_p[:,48,0]/1e5,c=colors[1],label='Nominal O$_2$ = $10^{-7}$')
ax2.plot(new_o3_bgoe[:,48,0],cam_o3_p[:,48,0]/1e5,c=colors[0],label='Nominal O$_2$ = $10^{-9}$')
ax2.legend(loc='upper left',fontsize=6)
xlims = ax2.get_xlim()
ax2.set_xlim(1e-22,xlims[1])

plt.tight_layout()
plt.savefig('figures/fig01.pdf')
plt.close()

cam_o3.close()
cam_ic.close()
