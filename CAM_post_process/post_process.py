from cdo import *
cdo = Cdo()
import argparse
import pathlib
import os
import pdb
import numpy as np
import aux_post_proc as app

parser = argparse.ArgumentParser()
parser.add_argument('sims_path',nargs=1,help='Path to CAM_output folder')
parser.add_argument('-s','--sim_name', nargs=1, default=['all'], help='name of simulation to process')
parser.add_argument('-i','--init_year', nargs=1, default=[31], type=int, help='first year of averaging')
parser.add_argument('-l','--last_year', nargs=1, default=[60], type=int, help='last year of averaging')
parser.add_argument('-w','--overwrite', action='store_true', help = 'force overwrite of output files')

args = parser.parse_args()

#path to simulation merged files from FRDR repository
CAM_output = pathlib.Path(args.sims_path[0])

#names of the simulations
#sim_list = ['ref_const_O3']
if args.sim_name[0] == 'all':
  sim_list = os.listdir(str(CAM_output))
else:
  sim_list = args.sim_name

overwrite = args.overwrite

#all simulations are averaged over years 31 - 60
start_year = 31
end_year = 60

#define the fields we need time mean and global files
fields_select = '-select,name=TS,cb_ozone_c,RELHUM,Q,CLDLOW,CLDHGH,LWCF,SWCF,TGCLDIWP,TGCLDLWP,LANDFRAC,IWC,LWC,T,CLOUD,U,V,ilev,CLDICE '


for sim in sim_list:
   #set the path to the simulation merged files
   merged_file = CAM_output / sim / (sim + '.cam.h0.merged_%04d_%04d.nc'%(start_year,end_year))
   if not merged_file.exists():
     raise IOError('File %s not found!'%(str(merged_file)))

   print("Simulation '%s'"%str(sim))
   print('Processing data from %s...'%str(merged_file))
   #post-processed data will go here
   output_path = pathlib.Path(sim)
   if not output_path.exists():
     output_path.mkdir(parents=True,exist_ok=True)


   #----Process time mean + zonal, global, tropical, polar means using CDO--------------------------

   #set name of time mean file
   output_timmean = output_path / (sim + '.cam.h0.timmean_%04d_%04d.nc'%(start_year,end_year))
   if not output_timmean.exists() or overwrite:
     #do time mean
     print('  Time mean...')
     cdo.timmean(input=fields_select+str(merged_file),output=str(output_timmean))

   #name of zonal mean file
   output_zonmean = output_path / (sim + '.cam.h0.zonmean_%04d_%04d.nc'%(start_year,end_year))
   if not output_zonmean.exists() or overwrite:
     #do zonal mean
     print('  Zonal mean...')
     cdo.zonmean(input='-select,name=RELHUM,Q,IWC,LWC,T,CLOUD,U,V,ilev '+str(output_timmean),output=str(output_zonmean))

   #name of global mean file
   output_globmean = output_path / (sim + '.cam.h0.globmean_%04d_%04d.nc'%(start_year,end_year))
   if not output_globmean.exists() or overwrite:
     #calc global mean
     print('  Global mean...')
     cdo.fldmean(input=fields_select+str(output_timmean),output=str(output_globmean))

   #name of tropical mean file
   output_tropmean = output_path / (sim + '.cam.h0.tropmean_%04d_%04d.nc'%(start_year,end_year))
   if not output_tropmean.exists() or overwrite:
     #do tropical mean (latitudes -15 to 15)
     print('  Tropical mean...')
     cdo.fldmean(input='-select,name=TS,T -sellonlatbox,0,360,-15,15 '+str(output_timmean),output=str(output_tropmean))

   #name of south polar mean file
   output_spolmean = output_path / (sim + '.cam.h0.spolmean_%04d_%04d.nc'%(start_year,end_year))
   if not output_spolmean.exists() or overwrite:
     #south polar mean (latitudes -90 to -60)
     print('  South polar mean...')
     cdo.fldmean(input='-select,name=TS -sellonlatbox,0,360,-90,-60 '+str(output_timmean),output=str(output_spolmean))

   #name of north polar mean file
   output_npolmean = output_path / (sim + '.cam.h0.npolmean_%04d_%04d.nc'%(start_year,end_year))
   if not output_npolmean.exists() or overwrite:
     #north polar mean (latitudes 60 to 90)
     print('  North polar mean...')
     cdo.fldmean(input='-select,name=TS -sellonlatbox,0,360,60,90 '+str(output_timmean),output=str(output_npolmean))
   #------------------------------------------------------------------------------------------------


   #----Calculate cloud top temperature and pressure------------------------------------------------

   #Here I calculate cloud top temperature and pressure as the values of these fields at the highest
   #level where the ice cloud mixing ratio exceeds some critical value (10^-8, in this case)
   output_globmean_cldt = output_path / (sim + '.cam.h0.globmean_%04d_%04d_cldtp.nc'%(start_year,end_year)) 
   if not output_globmean_cldt.exists() or overwrite:
     #first I create a temporary nc file that has the fields I need
     tmp_cldtp_file = output_path / (sim + '.cam.h0.timmean_%04d_%04d_cldtp.tmp.nc'%(start_year,end_year))
     cdo.select('name=TS,T,CLDICE',input=str(output_timmean),output=str(tmp_cldtp_file))

     #this estimates ice cloud top temperature and pressure and adds it to above file
     print('  Calculating cloud tops...')
     app.add_cloud_top_TP(tmp_cldtp_file,1e-8)

     #now I need to create a mask for NaN values (this is tricky to get right...)
     output_timmean_cldt = output_path / (sim + '.cam.h0.timmean_%04d_%04d_cldtp2.tmp.nc'%(start_year,end_year))
     cdo.setmissval('nan',input='-select,name=ICLDT,ICLDP '+str(tmp_cldtp_file),output=str(output_timmean_cldt))

     #and finally create the global mean 
     print('  ...and their global means...')
     cdo.fldmean(input=str(output_timmean_cldt),output=str(output_globmean_cldt))
   #------------------------------------------------------------------------------------------------


   #----Calculate LTS and EIS-----------------------------------------------------------------------

   #Calculates LTS and EIS from formulae in Wood & Bretherton 2006 (Journal of Climate, 19, p6425)
   #Uses LCL definition from Lawrence 2005 (BAMS, 86, p255)
   #The process is gross so I put it in the aux_post_proc library
   app.calc_lts_eis(merged_file,output_path,sim,overwrite)
   #------------------------------------------------------------------------------------------------


   #clean up (remove .tmp.nc files)
   os.system('rm '+str(output_path)+'/*.tmp.nc')

   print(' '*26 + '...Done!')


