import numpy as np
import netCDF4 as nc
import pdb
import pathlib
from cdo import *
cdo = Cdo()

#some Earth constants, etc
g = 9.80665
cp = 993.0
Rd = 287.058
Rv = 461.4
Lv = 2.45e6

def Temp_at_cloud_mmr(cmmr,T,P,cmmr0):
  #For a single column, find the highest level where
  #cloud mixing ratio reaches a critical value, cmmr0
  #Return the temperature and pressure of that level

  #first check if cmmr even reaches cmmr0
  if (cmmr<cmmr0).all():
    #return nans if column never gets that cloudy
    return np.nan, np.nan

  else:
    #iterate downwards from the top
    for i in np.arange(len(cmmr)-1):
      cmmr1 = cmmr[i]
      cmmr2 = cmmr[i+1]
      if cmmr1 < cmmr0 and cmmr2 > cmmr0:
        ilev = i
        break

    #interpolate to the midpoint between levels that span cmmr0
    Tnew = T[ilev] + (cmmr0 - cmmr[ilev]) * (T[ilev+1] - T[ilev]) / (cmmr[ilev+1] - cmmr[ilev])
    Pnew = P[ilev] + (cmmr0 - cmmr[ilev]) * (P[ilev+1] - P[ilev]) / (cmmr[ilev+1] - cmmr[ilev])
    return (Tnew, Pnew)


def add_cloud_top_TP(time_mean_file,mix_ratio_crit):
  #Calculates the cloud top temperature and pressure of all columns
  #adds these as new fields in time_mean_file

  #read in time mean data set 
  data = nc.Dataset(time_mean_file,'r+') #modifying nc files is dangerous but I accept the risks

  #fetch our grid info
  lons = data['lon'][:]
  lats = data['lat'][:]
  p = data['lev'][:]  #approximation to the pressure in hPa (not so good near surface)

  #ice cloud mixing ratio
  cmmr = data['CLDICE'][:].squeeze()

  #temperature
  T = data['T'][:].squeeze()

  #new data sets for cloud top temperature and pressure
  Tcloudt = np.zeros_like(data['TS'][:])
  Pcloudt = np.zeros_like(data['TS'][:])
  for ilat in np.arange(len(lats)):
      for ilon in np.arange(len(lons)):
        (Tcloudt[0,ilat,ilon], Pcloudt[0,ilat,ilon]) = Temp_at_cloud_mmr(cmmr[:,ilat,ilon],T[:,ilat,ilon],p,mix_ratio_crit)

  #create the necessary fields in the nc files
  if not 'ICLDT' in data.variables:
    data.createVariable('ICLDT','f8',('time','lat','lon'))
  if not 'ICLDP' in data.variables: 
    data.createVariable('ICLDP','f8',('time','lat','lon'))

  #set the data fields
  data['ICLDT'][:,:,:] = Tcloudt
  data['ICLDP'][:,:,:] = Pcloudt

  data.close()


def calc_lts_eis(merged_file,output_path,sim,overwrite):
  #Calculate LTS and EIS
  #The LCL is from Equation 24 of Lawrence 2005 (BAMS, 86, p255)
  #and EIS from Equations 4 and 5 of Wood & Bretherton 2006 (Journal of Climate, 19, p6425)

  print('\n  ****  LTS and EIS  ***')

  data = nc.Dataset(merged_file,'r')

  #get index and pressure of level closest to 700 hPa
  ilev = np.argmin(np.abs(data['lev'][:].squeeze()-700))
  plev = np.float(data['lev'][ilev])

  #get index of level closest to 850 hPa
  ilev850 = np.argmin(np.abs(data['lev'][:].squeeze()-850))

  #nothing else needed from this file
  data.close()

  lts_file = output_path / (sim + ".cam.h0.merged_0031_0060_lts.nc")
  eis_file = output_path / (sim + ".cam.h0.merged_0031_0060_eis.nc")

  #------LTS-------
  if not lts_file.exists() or overwrite:
    print('  * Calculating LTS... *')

    #potential temperature at 700 hPa
    pt_file = output_path / (sim + ".cam.h0.merged_0031_0060_pt700.tmp.nc")
    cdo.select("name=PT700",input="-sellevidx,%02d -expr,\'PT700=(T*(1000/%f)^0.286)\' "%(ilev+1,plev)+str(merged_file), output=str(pt_file))

    #potential temperature of surface
    pts_file = output_path / (sim + ".cam.h0.merged_0031_0060_ptsurf.tmp.nc")
    cdo.select("name=PTS",input=" -expr,\'PTS=(TS*(1e5/PS)^0.286)\' "+str(merged_file), output=str(pts_file))

    #merged the two PT files
    merge_pt = output_path / (sim + ".cam.h0.merged_0031_0060_ptmrg.tmp.nc")
    cdo.merge(input=str(pt_file)+" "+str(pts_file),output=str(merge_pt))

    #LTS calculation
    cdo.select("name=LTS",input="-expr,\'LTS=(PT700-PTS)\' "+str(merge_pt),output=str(lts_file))

  data_lts = nc.Dataset(str(lts_file),'r')

  if not eis_file.exists() or overwrite:
    #-----Z700-------
    print('  * Now z700...        *')

    #some file manipulation gymnastics to get CDO to do what I want 
    z700_file = output_path / (sim + ".cam.h0.merged_0031_0060_z700.tmp.nc")
    cdo.select("name=Z700",input="-sellevidx,%02d -expr,\'Z700=(Z3)\' "%(ilev+1)+str(merged_file), output=str(z700_file))
    data_z700 = nc.Dataset(str(z700_file),'r+')
    data_z700['lev'][:] = data_lts['lev'][:]
    data_z700.close()

    #-----LCL--------
    print('  * Now LCL...         *')

    #LCL from Lawrence 2005 and some file manipulation
    lcl_file = output_path / (sim + ".cam.h0.merged_0031_0060_lcl.tmp.nc")
    cdo.select("name=LCL",input="-sellevidx,30 -expr,\'LCL=((20+(TS-273.15)/5)*(100-RELHUM))\' "+str(merged_file), output=str(lcl_file))
    data_lcl = nc.Dataset(str(lcl_file),'r+')
    data_lcl['lev'][:] = data_lts['lev'][:]
    data_lcl.close()

    #-----Gamma850---
    print('  * Now Gamma850...    *')

    #saturation mixing ratio at 850 hPa
    smr850_file = output_path / (sim + ".cam.h0.merged_0031_0060_smr850.tmp.nc")
    cdo.select("name=SMR850",input="-sellevidx,%02d -expr,\'SMR850=(100*Q/RELHUM)\' "%(ilev850+1)+str(merged_file), output=str(smr850_file))

    #get surface temperature in usable format
    ts_file = output_path / (sim+".cam.h0.merged_0031_0060_ts.tmp.nc")
    cdo.select("name=TS",input=str(merged_file), output=str(ts_file))

    #same for temperature at 700 hPa
    t700_file = output_path / (sim+".cam.h0.merged_0031_0060_t700.tmp.nc")
    cdo.select("name=T700",input="-sellevidx,%02d -expr,\'T700=(T)\' "%(ilev+1)+str(merged_file), output=str(t700_file))

    #merge the temperature files
    merge_t = output_path / (sim + ".cam.h0.merged_0031_0060_tmrg.tmp.nc")
    cdo.merge(input=str(ts_file)+" "+str(t700_file),output=str(merge_t))

    #calculate temperature at 850 hPa
    t850_file = output_path / (sim + ".cam.h0.merged_0031_0060_t850.tmp.nc")
    cdo.select("name=T850",input="-expr,\'T850=(0.5*(TS+T700))\' "+str(merge_t), output=str(t850_file))

    #merge sat. mixing ratio and temperature files
    merge_850 = output_path / (sim + ".cam.h0.merged_0031_0060_850.tmp.nc")
    cdo.merge(input=str(smr850_file)+" "+str(t850_file),output=str(merge_850))

    #calculate moist lapse rate at 850 hPa and some file manipulation
    g850_file = output_path / (sim + ".cam.h0.merged_0031_0060_g850.tmp.nc")
    cdo.select("name=Gamma850",input="-expr,\'Gamma850=( (%f/%f) * (1-(1+%f*SMR850/(%f*T850)) / (1+%f^2*SMR850/(%f*%f*T850^2))))\' "
                %(g,cp,Lv,Rd,Lv,cp,Rv)+str(merge_850),output=str(g850_file))
    data_g850 = nc.Dataset(str(g850_file),'r+')
    data_g850['lev'][:] = data_lts['lev'][:]
    data_g850.close()

    #merge LTS, z700, LCL, lapse rate files
    fin_merge = output_path / (sim + ".cam.h0.merged_0031_0060_eismrg.tmp.nc")
    cdo.merge(input=str(lts_file)+" "+str(z700_file)+" "+str(lcl_file)+" "+str(g850_file),output=str(fin_merge))

    print('  * And finally EIS....*')

    #Now we can finally calculate EIS (Equation 4 of Wood & Bretherton 2006)
    cdo.select("name=EIS",input="-expr,\'EIS=(LTS-Gamma850*(Z700-LCL))\' "+str(fin_merge), output=str(eis_file))

    print('  **********************')
