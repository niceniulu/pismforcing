#!/usr/bin/env python3
# coding: utf-8

import netCDF4  as  nc
import numpy  as  np
import os
import calendar
import argparse


def  get_times(ystart, nyear):
    times = np.zeros((nyear,12))
    timebnds = np.zeros((nyear,12,2))

    total = 0
    for i in range(nyear):
        theyear = ystart + i
        if  calendar.isleap(theyear):
            days_in_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        for j in range(12):
            times[i,j] = total
            timebnds[i,j,0] = total
            total = total + days_in_month[j]
            timebnds[i,j,1] = total

    times = np.reshape(times, [-1])
    timebnds = np.reshape(timebnds, [-1,2])

    return times, timebnds

#----------------------------------------------------------
#  generate PISM climate forcing method 2
#----------------------------------------------------------
def  gen_pism_atm_given_3(temp2, precip, lat, lon, x, y, fout='forcing.nc', \
    tempunit="K", precipunit="kg/m**2s", version=1.2, timefreq="monthly", ystart=1850 ):
    #--------------------------
    # write to new output file
    #--------------------------
    # prepare output file: fout
    fo = nc.Dataset(fout,'w',format='NETCDF4_CLASSIC')
    fo.createDimension('x',len(x))
    fo.createDimension('y',len(y))
    fo.createDimension('time',None)
    fo.createDimension('bnds',2)
    fo.description = 'PISM atmosphere given forcing for version' + str(version)
    # x
    xs = fo.createVariable('x',np.float64, ('x',))
    xs.standard_name = "projection_x_coordinate"
    xs.long_name = "X-coordinate in Cartesian system"
    xs.units = "m"
    xs.axis = "X"
    xs[:] = x[:]
    # y
    ys = fo.createVariable('y',np.float64, ('y',))
    ys.standard_name = "projection_y_coordinate"
    ys.long_name = "Y-coordinate in Cartesian system"
    ys.units = "m"
    ys.axis = "Y"
    ys[:] = y[:]
    # time
    times = fo.createVariable('time',np.float64, ('time',))
    times.standard_name = 'time'
    times.long_name = 'time'
    times.units = 'days since '+str(ystart)+'-01-01 00:00:00'
    times.calendar = 'gregorian'  #'365_day'
    times.bounds = 'time_bnds'
    timebnds = fo.createVariable('time_bnds',np.float64, ('time','bnds'))

    if timefreq == "monthly":
        nyear = int(temp2.shape[0]/12.)
        times[:], timebnds[:] = get_times(ystart, nyear)

        #ntime = 12
        #times[:] = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5,319, 349.5]
        #timebnds[:ntime,:]=[0, 31, 31, 59, 59, 90, 90, 120, 120, 151, 151, 181, 181, 212, 212, 243, 243, 273, 273, 304, 304, 334, 334, 365.]

    else:
        print('not defined, fix me.')
    # lon
    lons = fo.createVariable('lon',np.float32, ('y','x'))
    lons.units = "degrees"
    lons.long_name = "longitude"
    lons.standard_name = "longitude"
    lons[:,:] = lon[:,:]
    # lat
    lats = fo.createVariable('lat',np.float32, ('y','x'))
    lats.units = "degrees"
    lats.long_name = "latitude"
    lats.standard_name = "latitude"
    lats[:,:] = lat[:,:]
    # air_temp
    air_temp = fo.createVariable('air_temp',np.float32, ('time','y','x'))
    air_temp.standard_name="air_temperature"
    air_temp.long_name = "Surface air temperature"
    air_temp.coordinates = "lon lat"
    if tempunit in ['K','Kelvin']:
        air_temp[...] = temp2-273.15
        air_temp.units = "degC"
    else:
        print("unknown temperature unit")

    # precipitation
    precipitation = fo.createVariable('precipitation',np.float32, ('time','y','x'))
    precipitation.standard_name = "lwe_precipitation_rate"
    precipitation.long_name = "Total Precipitation"
    precipitation.coordinates = "lon lat"
    if precipunit in ['kg/m**2s','kg m-2 s-1','kg m-2 second-1']:
        if version < 1.0 :
            precipitation[...] = precip*1000.*24.*3600/910.  # precipitation; kg/m^2s  =>  ice equivalent  mm/day
            precipitation.units = "mm day-1"
        elif version >= 1.0:
            precipitation[...] = precip
            precipitation.units = 'kg m-2 s-1'
    else:
        print("unknown precipitation unit")
        #
    print('=====> done!!!')
    return




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Input variables')
    parser.add_argument("fecham",type=str, help='Input climate forcing ECHAM file')
    parser.add_argument("fgrid",type=str, help='Input PISM grid file')
    parser.add_argument("--year", type=int, dest='year', default=1850, help='Which real year')
    parser.add_argument("--fout", type=str, dest='fout', default='out.nc', help='output file')
    parser.add_argument("--versionPISM",  type=float, dest='vPISM', default=1.2, help='The PISM version ')
    parser.add_argument("--SMB",  type=str, dest='smb', default='PDD', help='The smb scheme')
    parser.add_argument("--remap",  type=str, dest='remap', default='bil', help='The remap method')
    #parser.add_argument("--co2",  action='store_true', help='if add co2? default:False ')
    arguments = parser.parse_args()
    print(arguments)

    remap = arguments.remap
    fecham = arguments.fecham  #'./T63/1x_Ext_tund_003_echam_2900-2999.nc'
    fgrid = arguments.fgrid   #'../grids/pismr_antarctica_16km.griddes'
    
    #--------------------
    # interpolation
    #--------------------
    os.system("cdo  -remap"+remap+","+fgrid+ "   -selvar,temp2,tsurf,aprl,aprc,evap,geosp  "+fecham+"    regrid.nc")

    fi = nc.Dataset("regrid.nc",'r')
    temp2 = fi.variables['temp2'][:,:,:] # time,lat,lon
    evap = fi.variables['evap'][:,:,:] # time,lat,lon
    aprl = fi.variables['aprl'][:,:,:] # time,lat,lon
    aprc = fi.variables['aprc'][:,:,:] # time,lat,lon
    geosp = fi.variables['geosp'][:,:,:] # time,lat,lon

    lat = fi.variables['lat'][:,:]
    lon = fi.variables['lon'][:,:]
    x = fi.variables['x'][:]
    y = fi.variables['y'][:]

    fi.close()

    temp2 = temp2 + evap
    precip = aprl + aprc

    gen_pism_atm_given_3(temp2, precip, lat, lon, x, y, fout=arguments.fout, \
        tempunit="K", precipunit="kg/m**2s", version=arguments.vPISM, timefreq="monthly", ystart=arguments.year )

    os.system("rm  regrid.nc ")

    exit()











fgrid = "../pism_Greenland_topg_thk_bheatflx_5km.nc"
ftas = "/Users/lniu/LU/working/3_project_Greenland/Data/historical/tas_Amon_FIO-ESM-2-0_historical_r1i1p1f1_gn_185001-201412.nc"
fprecip = "/Users/lniu/LU/working/3_project_Greenland/Data/historical/pr_Amon_FIO-ESM-2-0_historical_r1i1p1f1_gn_185001-201412.nc"
#gridinfo = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0+y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
fout = "pism_forcing_historical_fioesm.nc"
version = 1.2



#--------------------
# read files
#--------------------
fi = nc.Dataset("regrid_tas.nc",'r')
temp2 = fi.variables['tas'][:,:,:] # time,lat,lon
fi.close()

fi = nc.Dataset("regrid_pr.nc",'r')
precip = fi.variables['pr'][:,:,:] # time,lat,lon
fi.close()

# grid file
fg = nc.Dataset(fgrid,'r')


gen_pism_atm_given_3(temp2, precip, lat, lon, x, y, fout=fout)
