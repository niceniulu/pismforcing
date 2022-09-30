#!/usr/bin/env python3
# coding: utf-8

import sys
import numpy as np
import os, argparse

import netCDF4   as nc

from datetime import datetime, timedelta
from netCDF4 import num2date, date2num

import yaml

import pyfesom2 as pf
from pyfesom2 import load_mesh, get_data

#sys.path.append(os.environ['FUNCTION_PATH']+'/pyfesom2')
#from pyfesom2 import load_mesh, get_data

#sys.path.append("/pf/a/a270075/esm-master/esm-runscripts/functions/externals/pyfesom2/")





###########################################################
# prepare forcing file for ice from FESOM2.0
# Lu Niu, AWI, 2020.06.16
# ---------------------------------------------------------
# On fesom origional grids,
# add mesh infomation and make cdo readable.
# ---------------------------------------------
# Input file info:
#   fintimes=[2030, 2031]:    timename of input files e.g.  var.fesom.2000.nc;
#   ntimeperfile=1:           how many timestep in one file;
#   version_FESOM=2.0:              FESOM version, now only 2.0 valid;
# Mesh info:
#   meshpath = 'path'         path of FESOM mesh;
#   abg=[50, 15, 90],         eular angles of FESOM mesh
# Ouput file info:
#   findir='./'
#   fout='fesom_for_ice.nc'   output file name;
#   ntimemean=1         calculate mean value over every ntimemean time steps;
#   timestartout=[2030,12,31,0,0]       date of first time step;
#   dtout=1                             delta T of output time step information;
#   dtoutunit                           unit of delta T output time step;
#   depthmax=400.                       maximum depth that selected;
#   depthmean=False,                     if calculate mean, i.e. 2D or 3D ocean;
#   ifrafesom=True,
#   fnameonly=None, (fesom.nc)     # has to be defined if ifrawfesom=False
#   ifverbose=False                     if print more information;
# #########################################################
def function_FESOM2ice(
            fintimes=[2030, 2031], ntimeperfile=1,  version_FESOM=2.0,
            meshpath = 'path', abg=[50, 15, -90],
            findir='./',
            fout='FESOM_for_ice.nc',
            ntimemean=1, timestartout=[2030,12,31,0,0], dtout=1 ,dtoutunit='years',
            depthmax=500., depthmin=300., depthmean=False,
            ifrawfesom=True, fnameonly=None,
            ifverbose=True):
    if str(version_FESOM) != "2.0" :
        raise ValueError('!!!!! only support for Fesom 2.0 now, stop !!!!!')
    #
    # predefined var names and time info
    varnames = ['temp', 'salt' ]
    varnamesout = ['temp', 'salt' ]  # fixed names for output file , do not edit
    units = ['degC', 'psu']
    unitsout = ['Kelvin', 'g/kg']   # fixed names for output file , do not edit
    timeunits = 'days since 1950-01-01 0:0:0'
    calendar = 'standard'
    nvar = len(varnames)
    #
    # number of time dimension
    ntime = len(fintimes) * ntimeperfile
    #
    # read mesh information
    mesh = load_mesh(meshpath, abg=abg, usepickle=False )
    lons = mesh.x2
    lats = mesh.y2
    depths = mesh.zlev  # negative values
    nlevel = mesh.nlev
    npoints = mesh.n2d
    #
    # read all data into a new array
    data = np.zeros((ntime, npoints, nlevel, nvar ))
    if ifrawfesom:
        print("O2I:     processing FESOM raw data from FESOM output directly.  ")
        for j in range(nvar):
            count=0
            for i in fintimes:
                fname = findir + '/' + varnames[j] + '.fesom.' + str(i) + '.nc'
                print("     ",varnames[j] , '.fesom.' , str(i) , '.nc')
                if  ifverbose: print('--------------- reading file:', fname)
                f = nc.Dataset(fname,'r')
                # Note: real data levels is 1 level less than mesh levels (/time,points,levels,vars/)
                data[count:count+ntimeperfile,:,:nlevel-1,j] = f.variables[varnames[j]][:]
                f.close()
                count = count+ntimeperfile

    elif not ifrawfesom:
        print("        processing catted one fesom file from", fnameonly)
        if fnameonly is None:
            raise ValueError("!!!!!! define the input file name, (ifrawfesom is False), stop !!!!")

            #
        fname = findir + '/' + fnameonly
        f = nc.Dataset(fname,'r')
        for j in range(nvar):
            count=0
            data[count:count+ntimeperfile,:,:nlevel-1,j] = f.variables[varnames[j]][:]
            #count = count+ntimeperfile
            #
        f.close()

    #
    # compute mean over time and depth
    # choose depths
    indexs = np.where((np.abs(depths[:]) <= np.abs(depthmax)) & (np.abs(depths[:]) >= np.abs(depthmin))  )[0]
    nindex = len(indexs)
    #
    if ifverbose: print('--------------- ',nindex, 'levels have been chosen:',depths[indexs])
    data2 = data[:,:,indexs[:],:]
    if depthmean:
        data3 = np.mean(data2, axis=2,keepdims=True)
    else:
        data3 = data2
    #
    # average times
    if ntimemean > 1:
        if ifverbose: print('--------------- averge of ', ntimemean, 'continues time steps')
        tmp = data3.reshape(-1, ntimemean, *data3.shape[-3:] )
        data4 = np.mean(tmp, axis=1)
        print(data4.shape)
    else:
        data4 = data3
    #
    # changes axies : change from (time,nod2,level,vars) to (time,level,nod2,vars)
    datafinal = np.swapaxes(data4, 1, 2)
    for i in range(nvar):
        if units[i] == ('degC' or 'C'):
            datafinal[:,:,:,i] = datafinal[:,:,:,i] + 273.15 # degC to Kelvin
            #

    ###--------------------------------------
    ### write to output file
    ###--------------------------------------
    fo = nc.Dataset(fout,'w')
    #
    fo.createDimension('nod2', data3.shape[1])
    fo.createDimension('level', data3.shape[2])
    fo.createDimension('time', None)
    # time
    time_out = fo.createVariable('time', np.float32, ('time',))
    time_out.units = timeunits
    time_out.calendar = calendar
    ntimeout = datafinal.shape[0]
    def compute_times():
        dates = []
        if dtoutunit == 'years':
            for i in range(ntimeout):
                dates.append(datetime(timestartout[0]+dtout*i, timestartout[1], timestartout[2]) )
        else:
            print('not supported yet, to be added,stop')
            return
            #
        times = date2num(dates,units = timeunits, calendar=calendar)
        return times
        #
    time_out[:] = compute_times()
    # level
    lev_out = fo.createVariable('level', np.float32, ('level',))
    lev_out.axis = 'Z'
    lev_out.positive ='down'
    lev_out.long_name = 'depths'
    lev_out.units = 'm'
    lev_out.description = 'depth downward, negative values'
    if  depthmean:
        lev_out[:] = depthmax
    else:
        lev_out[:] = depths[indexs]
    # longitude
    lon_out = fo.createVariable('lon', np.float32, ('nod2',))
    lon_out.long_name = 'longitude'
    lon_out.units = 'degree_east'
    lon_out[:] = lons
    # latitude
    lat_out = fo.createVariable('lat', np.float32, ('nod2',))
    lat_out.long_name = 'latitude'
    lat_out.units = 'degree_north'
    lat_out[:] = lats
    # data
    for j in range(nvar):
        data_out = fo.createVariable(varnamesout[j], np.float64, ('time','level', 'nod2',))
        data_out.coordinates = 'lon lat'
        data_out.units = unitsout[j]
        data_out[:] = datafinal[:,:,:,j]
    #
    fo.close()
    return




################################################################
# prepare ocean forcing file for PISM
# Lu Niu, 2020.06.22 @AWI
#----------------------------------------------
# take output file from function_FESOM2ice, as input file
#################################################################
def  function_ocean2PISM(fin='FESOM_for_ice.nc',ficegrid='grid_10km.griddes', fout='ocean_for_PISM.nc',
         switch_PISM='th', version_PISM=1.2, fpicobasin=None, interpmethod='nn' ):
    if  not  os.path.isfile(fin):
        raise ValueError(fin, "file is not found, should be generated from function_FESOM2ice, exit")


    if float(version_PISM) < 1.0 :
        raise ValueError("!!!!! Not tested yet for PISM version: "+ str(version_PISM) + " stop !!!!!")
    #
    # regrid from FESOM grid to PISM grid: regrid.nc
    #interpmethod='dis'    # Distance-weighted average remapping from cdo
    print("cdo interpolation to ice grids:", interpmethod)
    os.system("cdo   remap"+interpmethod+","+ficegrid+"   "+fin+"    regrid.nc")

    ### the variables that need to be selected in the end, see later codes
    selectvars = []

    ###----------------------------------------------
    if  switch_PISM == 'th' or switch_PISM == 'pico' :
        print('O2I:     Making forcing file for option: -ocean  '+switch_PISM)

        # read regrid.nc.
        f = nc.Dataset('regrid.nc','r+')
        temp = f.variables['temp']
        salt = f.variables['salt']
        depths = f.variables['level']

        # get rid of depth dimensions (time, level, **)
        # simple average over levels, but not weighted
        tempnew = np.mean(temp,axis=1, keepdims=False)
        saltnew = np.mean(salt,axis=1, keepdims=False)
        newdim = (temp.dimensions[0],temp.dimensions[2], temp.dimensions[3] )

        # temp
        tempout = f.createVariable('theta_ocean', temp.datatype, newdim )
        # copy attributes from original file
        tempout.setncatts({i: temp.getncattr(i)  for i in temp.ncattrs()})
        tempout[:] = tempnew
        # salt
        saltout = f.createVariable('salinity_ocean', salt.datatype, newdim  )
        saltout.setncatts({i: salt.getncattr(i)  for i in salt.ncattrs()})
        saltout[:] = saltnew

        ######
        selectvars.append('theta_ocean')
        selectvars.append('salinity_ocean')

        # for PISM pico :
        if switch_PISM == 'pico' :
            # Reading basin file for PISM pico:
            if fpicobasin == None:
                raise ValueError("!!!!! The pico basin file (fpicobasin) has to be defined. Stop !!!!!")
            # regrid
            os.system("cdo   remapnn,"+ficegrid+"   "+fpicobasin+"    regridbasin.nc") # use the nearest point interpolation
            f2 = nc.Dataset('regridbasin.nc','r')
            basins = f2.variables['basins']
            if len(basins[:].shape)  != 2:
                raise ValueError("!!!!! The basins variable is not 2 dimensions (lonxlat), time dim might exists, Stop !!!!!")
            #
            newdimlonlat = ( temp.dimensions[2], temp.dimensions[3] )
            basinsout = f.createVariable('basins', basins.datatype, newdimlonlat )
            # copy attributes from original file
            basinsout.setncatts({i: basins.getncattr(i)  for i in basins.ncattrs()})
            basinsout[:] = basins[:]
            f2.close()
            os.system("rm  regridbasin.nc")
            #
            selectvars.append('basins')
            #
        #
        f.close()
    else:
        raise ValueError("PISM option -ocean",switch_PISM, "or PISM version",version_PISM, " not supported yet, stop")


    ###------------------------------
    varsneed = ','.join(selectvars)
    # extrac new variables to the new output file
    print("\n........select the needed variables:", varsneed)
    os.system("ncks -O -v  "+ varsneed +"  regrid.nc   "+fout)

    ## clean regrid.nc
    os.system("rm  regrid.nc")

    print("O2I:     Done")

    return 0


def function_ocean_landmask( \
    finout='ocean_for_PISM_ref1.nc', \
    fslm2d='/pf/a/a270075/workPalmod/pool/pism/climates/awiesm-2.1/slm_LGM_GLAC1D_Ctrl.nc', \
    fgriddes='pismr_nhem_20km.griddes', \
    cdoremap='bil', \
    name_mask='slm', name_temp='theta_ocean', name_salt='salinity_ocean' ,\
    nsmooth=3 ):
    ######################
    print("O2I:     mask land with default values.")
    print("O2I:     remap with cdo",cdoremap,"  interpolate to target grid:", fgriddes)
    os.system("cdo   remap"+cdoremap+","+fgriddes+"   "+fslm2d+"    mask.nc")
    os.system("cp   "+finout+"   tmp.nc")

    fm = nc.Dataset('mask.nc','r')
    mm = fm.variables[name_mask][:]
    mm2 = np.where(mm<0.5, 0., 1.)

    fio = nc.Dataset('tmp.nc','r+')
    temp = fio.variables[name_temp]
    salt = fio.variables[name_salt]
    #
    print("O2I:     apply land points with temp:275.15, salt: 0. ")
    temp[:] = np.where(mm2>0.5, 275.15, temp)
    salt[:] = np.where(mm2>0.5, 0., salt)
    fio.close()

    print("O2I:     smoothing grids ")
    for i in np.arange(0, nsmooth):
        print("..... smooth, number=", i)
        os.system(" cdo  smooth9   tmp.nc   tmp2.nc")
        os.system("  mv  tmp2.nc tmp.nc ")

    os.system("rm   mask.nc ")
    os.system("mv  tmp.nc  " + finout)
    print("O2I:     Done.")

    return 0


def read_data(finname, var):
    with nc.Dataset(finname,'r') as fin:
        data = fin.variables[var]
        # print(data.dimensions)
        # print("nz1" in data.dimensions)
        # print(len(data.dimensions) == 3)
        # print(data.dimensions != ('time', 'nod2', 'nz1'))
        ndim = len(data.dimensions)
        # check
        if ndim == 3 :
            print("Read var >>",var,"<< with the dims:")
            print(data.dimensions)
            if data.dimensions != ('time', 'nod2', 'nz1'):
                raise ValueError('dimension is not (time nod2 nz1)')

        elif ndim == 2:
            print("Read var:",var," with the dims:")
            print(data.dimensions)
            if data.dimensions != ('time', 'nod2'):
                raise ValueError('dimension is not (time nod2)')

        data1 = fin.variables[var][:]
    return data1, ndim


def function_FESOM2ice_regular_grid(
            version_FESOM=2.0,
            meshpath = 'path', abg=[50, 15, -90],
            findir='./',
            fout='FESOM_for_ice.nc',
            depthmax=500., depthmin=300., 
            fnameonly=None,):
    
    if str(version_FESOM) != "2.0" :
        raise ValueError('!!!!! only support for Fesom 2.0 now, stop !!!!!')
    
    #######################################
    # predefined var names and time info
    varnames = ['temp', 'salt' ]
    varnamesout = ['temp', 'salt' ]  # fixed names for output file , do not edit
    units = ['degC', 'psu']
    unitsout = ['degC', 'g/kg']   # fixed names for output file , do not edit
    missingvalue = [0,33.]
    nvar = len(varnames)
    #
    
    finname = findir + './' + fnameonly
    foutname = fout
    depth_range = [depthmin, depthmax]
    nlon = 1440
    nlat = 720
    #######################################
    
    mesh = pf.load_mesh(meshpath,abg=abg, usepickle=False)
    
    # get the depth information
    if depth_range is None:
        depths = mesh.zlev[:-1]
    else:
        idx = (np.abs(mesh.zlev) >= depth_range[0]) & (np.abs(mesh.zlev) <= depth_range[1])
        depths = mesh.zlev[idx]

    nlev = len(depths)

    # the new regular grids:
    lon = np.linspace(-180, 180, nlon)
    lat = np.linspace(-90, 90, nlat)
    lons, lats = np.meshgrid(lon,lat)


    # read the time dim:
    fin = nc.Dataset(finname,'r')
    timein = fin.variables['time']
    ntime = len(timein)
 
    
    # write data
    with nc.Dataset(foutname, 'w') as fout:
        fout.createDimension('lat', nlat)
        latout = fout.createVariable('lat',np.float32, ('lat',))
        latout.units = 'degree_north'
        latout.standard_name = 'latitude'
        latout[:] = lat

        fout.createDimension('lon', nlon)
        lonout = fout.createVariable('lon',np.float32,('lon',))
        lonout.units = 'degree_east'
        lonout.standard_name = 'longitude'
        lonout[:] = lon

        fout.createDimension('time', None)
        timeout = fout.createVariable('time',np.float32,('time',))
        # timeout.setncatts(timein.__dict__)        
        timeout.standard_name = "time" 
        timeout.units = "years since 01-01-01 00:00:00"
        timeout.calendar = "365_day"
        timeout.axis = "T"
        

        fout.createDimension('level', None)
        levout = fout.createVariable('level',np.float32,('level',))
        levout.units = 'm'


        fout.description =  'interpolated fesom2.0 data'
        
        for i in range(nvar):
            thevar = varnames[i]
            print(">>>>>>>>>> processing >>> ",thevar)
            # read data
            data, ndim = read_data(finname, thevar)
            print(data.shape)

            varout = fout.createVariable(thevar,np.float32,('time','level','lat','lon'))
            varout.units = unitsout[i]

            for it in range(0,ntime):
                print("> time step ", it)
                if ndim == 3:
                    levout[:] = depths
                    for iz in range(0,nlev):
                        ilev = pf.ind_for_depth(depths[iz], mesh)
                        print("> level:", depths[iz])
                        level_data = data[it, :, ilev]
                        level_data[level_data==0] = np.nan
                        #idist_fast = pf.fesom2regular(level_data, mesh, lons, lats, how='idist', k=20)
                        idist_fast = pf.fesom2regular(level_data, mesh, lons, lats )
                        new = np.where(np.isnan(idist_fast), missingvalue[i], idist_fast)
                        varout[it,iz,:,:] = new

                elif ndim == 2:
                    levout[:] = 0. # meter
                    print("> surface")
                    level_data = data[it, :]
                    idist_fast = pf.fesom2regular(level_data, mesh, lons, lats )
                    new = np.where(np.isnan(idist_fast), missingvalue[i], idist_fast)
                    varout[it,0,:,:] = new

    # close input data
    fin.close()
    return

