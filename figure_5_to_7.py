import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from climt import RRTMGShortwave, RRTMGLongwave, get_default_state, get_grid
from matplotlib import rcParams
plt.rcParams.update({'font.size':8})

path_to_archive = "DG_Pre-Cambrian_O3_archive/CAM_input/"
cesm_dir = 'CAM_post_process/cesm_profs'

cesm_names = ["const_CO2_O2_1e-9",
              "const_CO2_O2_1e-7",
              "const_CO2_O2_1e-4",
              "const_CO2_O2_1e-3",
              "const_CO2_O2_1e-2",
              "ref_const_O3"]

#fields we need from CAM
fields = ['lev','TS','T','Q','ozone']

#create library of profiles from cesm data, global average
global_profs = []
for i in np.arange(len(cesm_names)):
    data = nc.Dataset(cesm_dir+cesm_names[i]+'/'+cesm_names[i]+'.cam.h0.globmean_0031_0060.nc','r')
    dtmp = {}
    dtmp['name'] = cesm_names[i]
    for f in fields:
        dtmp[f] = data[f][:].data.squeeze()

    global_profs.append(dtmp)

data.close()


#get sigma-pressure coordinate coefficients
cam_ic = nc.Dataset(path_to_archive+"cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc")
a = cam_ic['hyai'][:] * cam_ic['P0'][:]
b = cam_ic['hybi'][:]
P0 = cam_ic['P0'][0]


#construct RRTMG model
rad_sw = RRTMGShortwave(mcica=False)
rad_lw = RRTMGLongwave(mcica=False)

grid = get_grid(nx=1, ny=1, nz=len(global_profs[0]['lev']))
state = get_default_state([rad_sw, rad_lw], grid_state = grid)


#this naming convention is a nightmare
#these ones don't change between simulations
state['atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels'][:] = a[::-1]
state['atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels'][:] = b[::-1]
state['surface_air_pressure'][0] = P0
state['air_pressure'][:] = cam_ic['lev'][:][::-1,None,None]*100
state['air_pressure_on_interface_levels'][:] = cam_ic['ilev'][:][::-1,None,None]*100
state['mole_fraction_of_carbon_dioxide_in_air'][:] = 284.7*1e-6
state['zenith_angle'][0] = 60.0*np.pi/180.0
state['surface_albedo_for_diffuse_shortwave'][0] = 0.15
state['surface_albedo_for_direct_shortwave'][0] = 0.15
state['surface_albedo_for_diffuse_near_infrared'][0] = 0.15
state['surface_albedo_for_direct_near_infrared'][0] = 0.15
p = state['air_pressure'][:].copy()/100  #for short hand
p_int = state['air_pressure_on_interface_levels'][:].copy()/100  #for short hand

#Cloud setup: "case study" from Goldblatt and Zahnle 2011
cldtop = np.array([25000.,50000,70000])
cldbot = np.array([30000.,60000,85000])
fcld = np.array([0.25,0.25,0.4])
wp = np.array([20.0,25,40])
type = ['ice','liq','liq']
#below, set radii to constant (will use wp to zero these out ice/liquid for some layers)
state['cloud_water_droplet_radius'][:] = 11.0
state['cloud_ice_particle_size'][:] = 75.0


#These deal with the cloud overlap for the 8 columns needed to model 3 cloud layers
abc = np.zeros(3) #temporarily store cloud fraction for current iteration
ABC = np.array([[1,1,1],
                [1,1,0],
                [1,0,1],
                [0,1,1],
                [1,0,0],
                [0,1,0],
                [0,0,1],
                [0,0,0]])


#set the water paths for each cloud layer
for itype in np.arange(len(cldtop)):
    ilayers = np.where((state['air_pressure'][:]>=cldtop[itype]) &
                       (state['air_pressure'][:]<cldbot[itype]))
    abc[itype] = fcld[itype]
    if type[itype] == 'liq':
        state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'][ilayers] = wp[itype]/1000
        state['mass_content_of_cloud_ice_in_atmosphere_layer'][ilayers] = 0.0
    else:
        state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'][ilayers] = 0.0
        state['mass_content_of_cloud_ice_in_atmosphere_layer'][ilayers] = wp[itype]/1000


#these arrays will store the fluxes for each calculation
fup_sw_all = np.zeros((3,len(cesm_names),len(p_int)))
fdn_sw_all = np.zeros_like(fup_sw_all)
fup_lw_all = np.zeros_like(fup_sw_all)
fdn_lw_all = np.zeros_like(fup_sw_all)


iref = 5 #index of pal case

#index of tropopause, roughly
ipr = (np.squeeze(p) < 190).data

#Main calculations occur here
for i in np.arange(3):
    for j in np.arange(len(cesm_names)):
        state['mole_fraction_of_ozone_in_air'][:] = global_profs[iref]['ozone'][::-1,None,None]
        state['air_temperature'][:] = global_profs[iref]['T'][::-1,None,None]
        state['surface_temperature'][0] = global_profs[iref]['TS']
        state['specific_humidity'][:] = global_profs[iref]['Q'][::-1,None,None]

        if i == 0:
            #change only stratospheric humidity (pressures below ~200 hPa)
            state['specific_humidity'][ipr] = global_profs[j]['Q'][::-1][ipr][:,None,None]

        elif i == 1:
            #change ozone for the entire column
            state['mole_fraction_of_ozone_in_air'][:] = global_profs[j]['ozone'][::-1][:][:,None,None]

        elif i == 2:
            #change ozone for the entire column and stratospheric humidity
            state['mole_fraction_of_ozone_in_air'][:] = global_profs[j]['ozone'][::-1][:][:,None,None]
            state['specific_humidity'][ipr] = global_profs[j]['Q'][::-1][ipr][:,None,None]

        #---------------------------------------------------------------------------
        #temporarily store the fluxes for the 8 columns
        fup_sw_col = []
        fdn_sw_col = []
        fup_lw_col = []
        fdn_lw_col = []

        #set clouds for the 8 columns and compute fluxes
        for iABC in np.arange(8):
            hlayers = np.where((state['air_pressure'][:]>cldtop[0]) &
                               (state['air_pressure'][:]<cldbot[0]))
            mlayers = np.where((state['air_pressure'][:]>cldtop[1]) &
                               (state['air_pressure'][:]<cldbot[1]))
            llayers = np.where((state['air_pressure'][:]>cldtop[2]) &
                               (state['air_pressure'][:]<cldbot[2]))
            state['cloud_area_fraction_in_atmosphere_layer'][hlayers] = ABC[iABC,0]
            state['cloud_area_fraction_in_atmosphere_layer'][mlayers] = ABC[iABC,1]
            state['cloud_area_fraction_in_atmosphere_layer'][llayers] = ABC[iABC,2]

            out_sw = rad_sw(state)

            fup_sw_col.append(out_sw[1]['upwelling_shortwave_flux_in_air'])
            fdn_sw_col.append(out_sw[1]['downwelling_shortwave_flux_in_air'])

            out_lw = rad_lw(state)

            fup_lw_col.append(out_lw[1]['upwelling_longwave_flux_in_air'])
            fdn_lw_col.append(out_lw[1]['downwelling_longwave_flux_in_air'])

        #do weighted average over 8 columns
        fup_sw = abc[0]*abc[1]*abc[2] * fup_sw_col[0] + \
                 abc[0]*abc[1]*(1-abc[2]) * fup_sw_col[1] + \
                 abc[0]*(1-abc[1])*abc[2] * fup_sw_col[2] + \
                 (1-abc[0])*abc[1]*abc[2] * fup_sw_col[3] + \
                 abc[0]*(1-abc[1])*(1-abc[2]) * fup_sw_col[4] + \
                 (1-abc[0])*abc[1]*(1-abc[2]) * fup_sw_col[5] + \
                 (1-abc[0])*(1-abc[1])*abc[2] * fup_sw_col[6] + \
                 (1-abc[0])*(1-abc[1])*(1-abc[2]) * fup_sw_col[7]

        fdn_sw = abc[0]*abc[1]*abc[2] * fdn_sw_col[0] + \
                 abc[0]*abc[1]*(1-abc[2]) * fdn_sw_col[1] + \
                 abc[0]*(1-abc[1])*abc[2] * fdn_sw_col[2] + \
                 (1-abc[0])*abc[1]*abc[2] * fdn_sw_col[3] + \
                 abc[0]*(1-abc[1])*(1-abc[2]) * fdn_sw_col[4] + \
                 (1-abc[0])*abc[1]*(1-abc[2]) * fdn_sw_col[5] + \
                 (1-abc[0])*(1-abc[1])*abc[2] * fdn_sw_col[6] + \
                 (1-abc[0])*(1-abc[1])*(1-abc[2]) * fdn_sw_col[7]

        fup_lw = abc[0]*abc[1]*abc[2] * fup_lw_col[0] + \
                 abc[0]*abc[1]*(1-abc[2]) * fup_lw_col[1] + \
                 abc[0]*(1-abc[1])*abc[2] * fup_lw_col[2] + \
                 (1-abc[0])*abc[1]*abc[2] * fup_lw_col[3] + \
                 abc[0]*(1-abc[1])*(1-abc[2]) * fup_lw_col[4] + \
                 (1-abc[0])*abc[1]*(1-abc[2]) * fup_lw_col[5] + \
                 (1-abc[0])*(1-abc[1])*abc[2] * fup_lw_col[6] + \
                 (1-abc[0])*(1-abc[1])*(1-abc[2]) * fup_lw_col[7]

        fdn_lw = abc[0]*abc[1]*abc[2] * fdn_lw_col[0] + \
                 abc[0]*abc[1]*(1-abc[2]) * fdn_lw_col[1] + \
                 abc[0]*(1-abc[1])*abc[2] * fdn_lw_col[2] + \
                 (1-abc[0])*abc[1]*abc[2] * fdn_lw_col[3] + \
                 abc[0]*(1-abc[1])*(1-abc[2]) * fdn_lw_col[4] + \
                 (1-abc[0])*abc[1]*(1-abc[2]) * fdn_lw_col[5] + \
                 (1-abc[0])*(1-abc[1])*abc[2] * fdn_lw_col[6] + \
                 (1-abc[0])*(1-abc[1])*(1-abc[2]) * fdn_lw_col[7]

        #store the averaged flux
        #for SW, factor of 1/2 comes from average zenith angle (cos(60 deg))
        fup_sw_all[i,j,:] = fup_sw.squeeze()*0.5
        fdn_sw_all[i,j,:] = fdn_sw.squeeze()*0.5
        fup_lw_all[i,j,:] = fup_lw.squeeze()
        fdn_lw_all[i,j,:] = fdn_lw.squeeze()


#net flux = down - up by convention
fnet_sw_all = fdn_sw_all - fup_sw_all
fnet_lw_all = fdn_lw_all - fup_lw_all
fnet = fnet_sw_all + fnet_lw_all

#shorthand for PAL reference case
fup_sw_ref = fup_sw_all[0,-1]
fdn_sw_ref = fdn_sw_all[0,-1]
fup_lw_ref = fup_lw_all[0,-1]
fdn_lw_ref = fdn_lw_all[0,-1]
fnet_sw_ref = fnet_sw_all[0,-1]
fnet_lw_ref = fnet_lw_all[0,-1]

#nominal O2 for plotting
o2_array = np.array([4.5545e-9,6.059e-8,7.255e-5,9.082e-4,9.941e-3,0.2112])


#pseudo tropopause forcing
#call tropopause where net flux in PAL case falls below 5 W m^-2
itrop = np.where(np.abs(fnet[0,-1,:])<5)[0][0] + 1

forcing_sw_tp = fnet_sw_all[:,:,itrop] - fnet_sw_all[:,-1,:][:,None,itrop]
forcing_lw_tp = fnet_lw_all[:,:,itrop] - fnet_lw_all[:,-1,:][:,None,itrop]
forcing_tot_tp = fnet[:,:,itrop] - fnet[:,-1,:][:,None,itrop]
forcing_sw_tp[:,-1] = 0.0
forcing_lw_tp[:,-1] = 0.0
forcing_tot_tp[:,-1] = 0.0

#repeat at toa
itrop = -1

forcing_sw_toa = fnet_sw_all[:,:,itrop] - fnet_sw_all[:,-1,:][:,None,itrop]
forcing_lw_toa = fnet_lw_all[:,:,itrop] - fnet_lw_all[:,-1,:][:,None,itrop]
forcing_tot_toa = fnet[:,:,itrop] - fnet[:,-1,:][:,None,itrop]
forcing_sw_toa[:,-1] = 0.0
forcing_lw_toa[:,-1] = 0.0
forcing_tot_toa[:,-1] = 0.0


#now we make the plots
cmap = plt.cm.seismic_r
colors = [cmap(0.8),cmap(0.6),cmap(0.4),cmap(0.2),cmap(0.),'k']
o2_level = ['Nominal O$_2=$ PAL','Nominal O$_2=10^{-2}$','Nominal O$_2=10^{-3}$','Nominal O$_2=10^{-4}$','Nominal O$_2=10^{-7}$','Nominal O$_2=10^{-9}$']

lwidths = [3,1.5,2,2,2,2]
filenames = ['fig06.pdf','fig05.pdf','fig07.pdf']
title = ['Forcing due to stratospheric humidity',
         'Forcing due to ozone',
         'Forcing due to ozone and stratospheric humidity']
dyn_range = [6.5,18,20]  #dynamic range for x axis of difference plots
lines = []
for ipert in np.arange(3):
    fig, axes = plt.subplots(ncols=3,nrows=5,figsize=(7.5,10))
    fig.suptitle(title[ipert],fontsize=12)

    line, = axes[0][0].plot(fup_sw_ref,p_int.squeeze(),'k',label=o2_level[-1],zorder=100)
    lines.append(line)
    axes[0][0].set(yscale='log')
    axes[0][0].set_xlabel('SW upward flux (W m$^{-2}$)',fontsize=10)
    axes[0][0].set_ylabel('Pressure (hPa)',fontsize=10)
    axes[0][0].invert_yaxis()
    inset00 = axes[0][0].inset_axes([0.2,0.45,0.5,0.5])
    inset00.semilogy(fup_sw_ref,p_int.squeeze(),'k',zorder=100)
    inset00.set_ylim((300,3.5))
    inset00.set_xlim((115,130))

    axes[0][1].plot(fdn_sw_ref,p_int.squeeze(),'k',zorder=100)
    axes[0][1].set(yscale='log')
    axes[0][1].set_xlabel('SW downward flux (W m$^{-2}$)',fontsize=10)
    axes[0][1].invert_yaxis()
    inset01 = axes[0][1].inset_axes([0.2,0.45,0.5,0.5])
    inset01.semilogy(fdn_sw_ref,p_int.squeeze(),'k',zorder=100)
    inset01.set_ylim((300,3.5))
    inset01.set_xlim((330,355))

    axes[0][2].plot(fnet_sw_ref,p_int.squeeze(),'k',zorder=100)
    axes[0][2].set(yscale='log')
    axes[0][2].set_xlabel('SW net flux (W m$^{-2}$)',fontsize=10)
    axes[0][2].invert_yaxis()
    inset02 = axes[0][2].inset_axes([0.2,0.45,0.5,0.5])
    inset02.semilogy(fnet_sw_ref,p_int.squeeze(),'k',zorder=100)
    inset02.set_ylim((300,3.5))
    inset02.set_xlim((212,230))

    axes[2][0].plot(fup_lw_ref,p_int.squeeze(),'k',zorder=100)
    axes[2][0].set(yscale='log')
    axes[2][0].set_xlabel('LW upward flux (W m$^{-2}$)',fontsize=10)
    axes[2][0].set_ylabel('Pressure (hPa)',fontsize=10)
    axes[2][0].invert_yaxis()
    inset20 = axes[2][0].inset_axes([0.4,0.45,0.5,0.5])
    inset20.semilogy(fup_lw_ref,p_int.squeeze(),'k',zorder=100)
    inset20.set_ylim((300,3.5))
    inset20.set_xlim((228,240))

    axes[2][1].plot(fdn_lw_ref,p_int.squeeze(),'k',zorder=100)
    axes[2][1].set(yscale='log')
    axes[2][1].set_xlabel('LW downward flux (W m$^{-2}$)',fontsize=10)
    axes[2][1].invert_yaxis()
    inset21 = axes[2][1].inset_axes([0.45,0.45,0.5,0.5])
    inset21.semilogy(fdn_lw_ref,p_int.squeeze(),'k',zorder=100)
    inset21.set_ylim((300,3.5))
    inset21.set_xlim((-2,30))

    axes[2][2].plot(fnet_lw_ref,p_int.squeeze(),'k',zorder=100)
    axes[2][2].set(yscale='log')
    axes[2][2].set_xlabel('LW net flux (W m$^{-2}$)',fontsize=10)
    axes[2][2].invert_yaxis()
    inset22 = axes[2][2].inset_axes([0.4,0.45,0.5,0.5])
    inset22.semilogy(fnet_lw_ref,p_int.squeeze(),'k',zorder=100)
    inset22.set_ylim((300,3.5))
    inset22.set_xlim((-235,-210))

    axes[1][0].plot(fup_sw_ref-fup_sw_ref,p_int.squeeze(),'k')
    axes[1][0].set(yscale='log')
    axes[1][0].set_xlabel('$\Delta$SW upward flux (W m$^{-2}$)',fontsize=10)
    axes[1][0].set_ylabel('Pressure (hPa)',fontsize=10)
    axes[1][0].invert_yaxis()
    xlims = axes[1][0].get_xlim()
    axes[1][0].set_xlim((0.5*(xlims[0]+xlims[1]-dyn_range[ipert]),0.5*(xlims[0]+xlims[1]+dyn_range[ipert])))

    axes[1][1].plot(fdn_sw_ref-fdn_sw_ref,p_int.squeeze(),'k')
    axes[1][1].set(yscale='log')
    axes[1][1].set_xlabel('$\Delta$SW downward flux (W m$^{-2}$)',fontsize=10)
    axes[1][1].invert_yaxis()
    xlims = axes[1][1].get_xlim()
    axes[1][1].set_xlim((0.5*(xlims[0]+xlims[1]-dyn_range[ipert]),0.5*(xlims[0]+xlims[1]+dyn_range[ipert])))

    axes[1][2].plot(fnet_sw_ref-fnet_sw_ref,p_int.squeeze(),'k')
    axes[1][2].set(yscale='log')
    axes[1][2].set_xlabel('$\Delta$SW net flux (W m$^{-2}$)',fontsize=10)
    axes[1][2].invert_yaxis()
    xlims = axes[1][2].get_xlim()
    axes[1][2].set_xlim((0.5*(xlims[0]+xlims[1]-dyn_range[ipert]),0.5*(xlims[0]+xlims[1]+dyn_range[ipert])))

    axes[3][0].plot(fup_lw_ref-fup_lw_ref,p_int.squeeze(),'k')
    axes[3][0].set(yscale='log')
    axes[3][0].set_xlabel('$\Delta$LW upward flux (W m$^{-2}$)',fontsize=10)
    axes[3][0].set_ylabel('Pressure (hPa)',fontsize=10)
    axes[3][0].invert_yaxis()
    xlims = axes[3][0].get_xlim()
    axes[3][0].set_xlim((0.5*(xlims[0]+xlims[1]-dyn_range[ipert]),0.5*(xlims[0]+xlims[1]+dyn_range[ipert])))

    axes[3][1].plot(fdn_lw_ref-fdn_lw_ref,p_int.squeeze(),'k')
    axes[3][1].set(yscale='log')
    axes[3][1].set_xlabel('$\Delta$LW downward flux (W m$^{-2}$)',fontsize=10)
    axes[3][1].invert_yaxis()
    xlims = axes[3][1].get_xlim()
    axes[3][1].set_xlim((0.5*(xlims[0]+xlims[1]-dyn_range[ipert]),0.5*(xlims[0]+xlims[1]+dyn_range[ipert])))

    axes[3][2].plot(fnet_lw_ref-fnet_lw_ref,p_int.squeeze(),'k')
    axes[3][2].set(yscale='log')
    axes[3][2].set_xlabel('$\Delta$LW net flux (W m$^{-2}$)',fontsize=10)
    axes[3][2].invert_yaxis()
    xlims = axes[3][2].get_xlim()
    axes[3][2].set_xlim((0.5*(xlims[0]+xlims[1]-dyn_range[ipert]),0.5*(xlims[0]+xlims[1]+dyn_range[ipert])))

    #insets
    for icase in np.arange(4,-1,-1):
        line, = axes[0][0].plot(fup_sw_all[ipert,icase],p_int.squeeze(),c=colors[icase],label=o2_level[icase],lw=lwidths[icase],zorder=icase)
        lines.append(line)
        inset00.plot(fup_sw_all[ipert,icase],p_int.squeeze(),c=colors[icase],label=o2_level[icase],lw=lwidths[icase],zorder=icase)

        axes[0][1].plot(fdn_sw_all[ipert,icase],p_int.squeeze(),c=colors[icase],label=o2_level[icase],lw=lwidths[icase],zorder=icase)
        inset01.plot(fdn_sw_all[ipert,icase],p_int.squeeze(),c=colors[icase],label=o2_level[icase],lw=lwidths[icase],zorder=icase)

        axes[0][2].plot(fnet_sw_all[ipert,icase],p_int.squeeze(),c=colors[icase],label=o2_level[icase],lw=lwidths[icase],zorder=icase)
        inset02.plot(fnet_sw_all[ipert,icase],p_int.squeeze(),c=colors[icase],label=o2_level[icase],lw=lwidths[icase],zorder=icase)

        axes[2][0].plot(fup_lw_all[ipert,icase],p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        inset20.plot(fup_lw_all[ipert,icase],p_int.squeeze(),c=colors[icase],label=o2_level[icase],lw=lwidths[icase],zorder=icase)

        axes[2][1].plot(fdn_lw_all[ipert,icase],p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        inset21.plot(fdn_lw_all[ipert,icase],p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)

        axes[2][2].plot(fnet_lw_all[ipert,icase],p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        inset22.plot(fnet_lw_all[ipert,icase],p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)

        axes[1][0].plot(fup_sw_all[ipert,icase]-fup_sw_ref,p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        axes[1][1].plot(fdn_sw_all[ipert,icase]-fdn_sw_ref,p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        axes[1][2].plot(fnet_sw_all[ipert,icase]-fnet_sw_ref,p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        axes[3][0].plot(fup_lw_all[ipert,icase]-fup_lw_ref,p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        axes[3][1].plot(fdn_lw_all[ipert,icase]-fdn_lw_ref,p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)
        axes[3][2].plot(fnet_lw_all[ipert,icase]-fnet_lw_ref,p_int.squeeze(),c=colors[icase],lw=lwidths[icase],zorder=icase)

    for irow in np.arange(4):
        for icol in np.arange(3):
            axes[irow][icol].hlines(200,axes[irow][icol].get_xlim()[0],axes[irow][icol].get_xlim()[1],colors='0.5',linestyles=':')

    axes[4][2].legend(lines,o2_level,loc='upper right',fontsize=7)
    axes[4][2].get_xaxis().set_visible(False)
    axes[4][2].get_yaxis().set_visible(False)
    axes[4][2].axis('off')

    ax = axes[4][0]
    ax.semilogx(o2_array[:],forcing_sw_toa[ipert][:],'^',color='darkviolet',ms=6,linestyle='-',label='SW')
    ax.plot(o2_array[:],forcing_lw_toa[ipert][:],'^',color='seagreen',ms=6,linestyle='-',label='LW')
    ax.plot(o2_array[:],forcing_tot_toa[ipert][:],'k^',ms=6,linestyle='-',label='SW+LW')
    ax.set(xlim=(1e-9,5e-1))
    ax.set_xlabel('Nominal O$_2$ level (mol mol$^{-1}$)',fontsize=10)
    ax.set_ylabel('TOA Forcing (W m$^{-2}$)',fontsize=10)
    ax.hlines(0,1e-9,1e-1,colors='0.5',linestyles=':')
    ax.legend(loc = 'best',fontsize=6)

    ax = axes[4][1]
    ax.semilogx(o2_array[:],forcing_sw_tp[ipert][:],'^',color='darkviolet',ms=6,linestyle='-',label='SW')
    ax.plot(o2_array[:],forcing_lw_tp[ipert][:],'^',color='seagreen',ms=6,linestyle='-',label='LW')
    ax.plot(o2_array[:],forcing_tot_tp[ipert][:],'k^',ms=6,linestyle='-',label='SW+LW')
    ax.set(xlim=(1e-9,5e-1))
    ax.set_xlabel('Nominal O$_2$ level (mol mol$^{-1}$)',fontsize=10)
    ax.set_ylabel('200 hPa Forcing (W m$^{-2}$)',fontsize=10)
    ax.hlines(0,1e-9,1e-1,colors='0.5',linestyles=':')

    plt.tight_layout()
    plt.savefig('figures/'+filenames[ipert])
    plt.close()
