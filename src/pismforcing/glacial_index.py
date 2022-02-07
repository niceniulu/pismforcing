#!/usr/bin/env python3
# coding: utf-8


import sys , os
import logging

import numpy as np
import netCDF4   as nc
from pathlib import Path
from scipy import interpolate
import matplotlib.pyplot as plt



## import own packages
from . import select_vars_for_PISM, downscale_ECHAM2PISMgrids, run_dEBM, generate_atmos_forcing_for_PISM
from . import function_FESOM2ice, function_ocean2PISM, function_ocean_landmask


###############################################################################
# main functin for the Glacial index method for generating forcing files for both atmosphere and ocean
###############################################################################
'''
theyear: the chosen year
fpism: the current pism output file
'''
def main_index_forcing_file(theyear, fpism,  **settings):
    #####
    #--------------- get the index --------------------------------------
    #####
    print("\n\nStep 1 ===========================================")
    print("             get the index at year:", theyear)
    print("=======================================================\n\n" )
    theindex = get_index(theyear,settings['refyear0'],settings['refyear1'], settings['index_file'],
                dt_interp=settings['index_interp_dt'], varname=settings['index_varname'],errorrefyear=settings['index_errorrefyear'], ifplot=False )
    print("\n..... the index value is: ", theindex, "\n\n")

    #####
    #--------------- the atmosphere interface --------------------------------------
    #####
    print("\n\nStep 2 ===========================================")
    print("             Prepare atmosphere forcing file")
    print("=======================================================\n\n" )
    file_atmos = 'atmosphere_forcing_for_PISM_' + str(theyear)+'.nc'
    glacial_index_main_atmosphere_file_for_PISM(theyear=theyear, index=theindex, file_output=file_atmos, file_pism=fpism, **settings)

    # write to forcing.txt file
    with open("forcing.txt", "w") as finfo:
        finfo.write("fname_atmos="+file_atmos+"\n")



    #####
    #--------------- the ocean interface --------------------------------------
    #####
    if settings['ocean_couple']:
        print("\n\nStep 3 ===========================================")
        print("        Enable ocean couple, Prepare ocean forcing file")
        print("=======================================================\n\n" )
        if settings['FESOM_interp'] :
            # add new dict items, define file names of ocean forcing to PISM grids.
            settings['PISM_file_ocean_forcing_ref0'] = "ocean_for_PISM_ref0.nc"
            settings['PISM_file_ocean_forcing_ref1'] = "ocean_for_PISM_ref1.nc"
            #
            bool0 = os.path.isfile(settings['PISM_file_ocean_forcing_ref0'])
            bool1 = os.path.isfile(settings['PISM_file_ocean_forcing_ref1'])
            if bool0 and bool1:
                print("..... PISM_file_ocean_forcing_ref0:",settings['PISM_file_ocean_forcing_ref0'],  \
                                    "PISM_file_ocean_forcing_ref1:", settings['PISM_file_ocean_forcing_ref1']," exits")
            else:
                print("..... PISM_file_ocean_forcing_ref0:",settings['PISM_file_ocean_forcing_ref0'],  \
                                    "PISM_file_ocean_forcing_ref1:", settings['PISM_file_ocean_forcing_ref1']," NOT exits")
                print("         FESOM file has to be interpolate to PISM grids for the first step >>>>> ")
                initial_FESOM2PISMgrids(**settings)
                #
            ff_ref0 = settings['PISM_file_ocean_forcing_ref0']
            ff_ref1 = settings['PISM_file_ocean_forcing_ref1']
        else:
            ff_ref0 = settings['ocean_file_dir'] + settings['ocean_filename_ref0']
            ff_ref1 = settings['ocean_file_dir'] + settings['ocean_filename_ref1']
            #
        print("..... The ref0 file: ", ff_ref0 )
        print("..... The ref1 file: ", ff_ref1 )
        #
        file_ocean = 'ocean_forcing_for_PISM_'+ settings['PISM_ocean_switch']+'_'+str(theyear)+".nc"
        print("..... The output ocean file for PISM : ", file_ocean )
        get_index_data(theindex,ff_ref0, ff_ref1, fileout=file_ocean)

        # write to forcing.txt file
        with open("forcing.txt", "a") as finfo:
            finfo.write("fname_ocean="+file_ocean+"\n")
        #
    return  file_atmos, file_ocean
###############################################################################











def get_index(theyear,refyear0, refyear1, file, dt_interp=1, varname=None, errorrefyear=500, ifplot=False):
    '''
    # Get the index value based on a timeseries.
    # return: theindex
    # ------------------------------------------
    # theyear: choose the year.  #
    #refyear0: the reference year that the index is 0.
    #refyear1: the reference year that the index is 1.
    #file: the file of the timeseries, can be .nc, .dat, .txt.
    #dt_interp: interpolate the original timeseries to a smaller timeinterval with THIS certain timestep.
    #varname:  the variable name, required when open an netcdf file.
    #errorrefyear: refyear-errorrefyear to reyear+errorrefyear, for calcuating the mean over the refyear times.
    #ifplot: if plot for the timeseries.
    #------------------------------------------
    '''

    #--------------------------- read files -------------------------------------------
    if Path(file).suffix == '.nc':
        if varname is None:
            raise ValueError("!!!!! The varname (variable name of the time series) has to be defined! Stop !!!!!!")

        print('I:       reading the netcdf file:' + file )
        f = nc.Dataset(file, 'r')
        time = f.variables['time'][:]
        var = f.variables[varname][:]
        f.close()

    elif (Path(file).suffix == '.dat') or  (Path(file).suffix == '.txt'):
        f = open(file,'r')
        rr = f.read().split()
        data = []
        for dd in rr:
            data.append(float(dd))
        #
        data = np.reshape(data,(-1,2))
        time = data[:,0]
        var = data[:,1]
        f.close()
    else:
        raise ValueError("!!!!! unknow file format, stop !!!!!")
        #
    #---------------------------------------------------
    yearmin = np.min(time)
    yearmax = np.max(time)
    # check if the chosen years are in the time range
    if (refyear0 > yearmax or refyear0 < yearmin):
        raise ValueError("!!!!!  the refyear0:", refyear0 ," is out of range, stop !!!!!")
        #
    if (refyear1 > yearmax or refyear1 < yearmin):
        raise ValueError("!!!!!  the refyear1:", refyear1 ," is out of range, stop !!!!!")
        #
    #--------------- interpolate to smaller timestep (dt_interp) data --------------------------------------
    newyears = np.arange(yearmin, yearmax+dt_interp,dt_interp)
    # interpolation using scipy.interpolate , default is linear interpolation
    ff = interpolate.interp1d(time, var)
    newvars = ff(newyears)
    #
    if ifplot:
        plt.plot(time, var, 'o', newyears, newvars, '-' )
        plt.savefig('timeseries.png')
        plt.show()
    #
    #--------------- calculate the index based on the original timeseries ---------------------------------------
    # calculate the value of the timeseries at refyear0 (average over a certain time period)
    yearlow = refyear0-errorrefyear ; yearup = refyear0+errorrefyear
    varidx0 = (newyears>= yearlow) & (newyears<=  yearup)
    var0 = np.mean( newvars[varidx0] )
    print("I:       The value at refyear0" , refyear0,"( average over [" , yearlow, ",", yearup, "]) is :", var0 )
    # calculate the value of the timeseries at refyear1
    yearlow = refyear1-errorrefyear ; yearup = refyear1+errorrefyear
    varidx1 = (newyears>=yearlow ) & (newyears<=yearup )
    var1 = np.mean( newvars[varidx1] )
    print("I:       The value at refyear1" , refyear1," ( average over [" , yearlow, ",", yearup, "]) is :", var1 )
    # the value at the chosen year:
    theidx = np.where(newyears==theyear)
    thevar = newvars[theidx]
    print("I:       The value at the chosen year" , theyear , " is :", thevar)

    theindex = (thevar-var0)/ (var1-var0)
    print("I:       Done. The output (index value):" , theindex)

    return theindex





###############################################################################
# main function for making the atmospheric forcing files
###############################################################################

def glacial_index_main_atmosphere_file_for_PISM(theyear=None,
                                                index=None,
                                                file_output='atmosphere_forcing_for_PISM.nc',
                                                file_pism='20km_lrno_pismr_out_-021100--021001.nc',
                                                ECHAM6_file_ref0=None,
                                                ECHAM6_file_ref1=None,
                                                PISM_file_griddes=None,
                                                PISM_smb_scheme=None,
                                                atmos_downscale_interp='bil',
                                                atmos_lapserate_T='linear',
                                                atmos_lapserate_T_value=0,
                                                atmos_lapserate_P='desertification',
                                                binary_dEBM='./dEBMmain',
                                                time_frequency='12month',
                                                PISM_lakestore=False,
                                                kyr_dEBM=None,
                                                **opt):  ## other defines, but not used

    print()
    #------ the least required defines
    if ECHAM6_file_ref0 is None:    raise ValueError("ECHAM6_file_ref0 needs to be defined. ")
    if ECHAM6_file_ref1 is None:    raise ValueError("ECHAM6_file_ref1 needs to be defined. ")
    if PISM_file_griddes is None:    raise ValueError("PISM_file_griddes needs to be defined. ")
    if PISM_smb_scheme is None:    raise ValueError("PISM_smb_scheme needs to be defined. ")


    #----- make new ECHAM6 file based on the index
    print()
    print("I:       Generate dummy ECHAM6 file based on the index: ", index ,"\n")
    if kyr_dEBM is None:
        year_kyr = theyear/1000.
    else:
        year_kyr = kyr_dEBM

    print("         year in kyr:", year_kyr)
    print()

    obliquity = modify_data_ECHAM6(year_kyr, index, ECHAM6_file_ref0, ECHAM6_file_ref1, "dummy_ECHAM6.nc")

    print("I:        The obliquity used for dEBM is: ", obliquity )

    #---------------------------------------------------------------
    # select the necessary variables used for different SMB scheme.
    print("I:       Select needed variables from ECHAM6 format file (select_vars_for_PISM) |  output: select_ECHAM6.nc \n")
    select_vars_for_PISM(fin='dummy_ECHAM6.nc', schemeSMB=PISM_smb_scheme, lakestore=PISM_lakestore, fout='select_ECHAM6.nc')


    #-------------------------------------------------------------
    # downscale:  interpolat from ECHAM6 grids to PISM grids
    print("I:       Interpolate (downscale) ECHAM6 format file to PISM grids (downscale_ECHAM2PISMgrids) |  output: downscaled.nc ")
    downscale_ECHAM2PISMgrids(fecham='select_ECHAM6.nc', fpismgrid=PISM_file_griddes, fpismcurrent=file_pism,
                                interpmethod=atmos_downscale_interp ,
                                downscaleT=atmos_lapserate_T, factorT=atmos_lapserate_T_value,
                                downscaleP=atmos_lapserate_P, factorP=None,
                                fout='downscaled.nc' )


    #-------------------------------------------------------------
    # generate the forcing file for PISM

    if PISM_smb_scheme == 'dEBM':
        print("I:       Run binary dEBM for surface mass balance. | output: surface_mass_balance.nc ")
        run_dEBM(fin='downscaled.nc',obliquity=obliquity, exedEBM=binary_dEBM)
        #
        print("I:       Make atmosphere forcing file for PISM.  ")
        generate_atmos_forcing_for_PISM(fin='downscaled.nc',fout=file_output,scheme='dEBM',
                                        fileplus1='surface_mass_balance.nc',timefreq = time_frequency,lakestore=PISM_lakestore)

        print("        Done. need to set PISM switch: -surface given  ")
    elif PISM_smb_scheme == 'PDD':
        print("I:       Make atmosphere forcing file for PISM.  ")
        generate_atmos_forcing_for_PISM(fin='downscaled.nc',fout=file_output,scheme='PDD',
                                        fileplus1=None,timefreq = time_frequency)
        print("       Done. need to set PISM switch: -atmophere given  -surface pdd ")
    else:
        raise ValueError(" unknown surface mass balance scheme.  ")

    return


###############################################################################
#   make dummy ECHAM6 file at chosen year
###############################################################################
#
## Modified from Uta's script (pism_index_script.py)
## function for making the atmosphere forcing from ECHAM6
#  Note: (Difference to Function index_forcing)
#  - diff1: not using xarray because of some errors occurred (e.g. the datetime info). Fix it if possible.
#  - diff2: all other variables are scaled.
#
def  modify_data_ECHAM6(kyr, index, refdata0, refdata1, fileout ):
    # Import Third-Party Packages
    #import climlab
    from climlab.solar.insolation import daily_insolation
    from climlab.solar.orbital import OrbitalTable

    f0 = nc.Dataset(refdata0, 'r')
    f1 = nc.Dataset(refdata1,'r')

    #-------------------------------------------------------------------------
    # check if the following variabels is existed, which might be used by dEBM
    #-------------------------------------------------------------------------
    check_vars = ['lat','lon','temp2', 'aprc','aprl','geosp','srad0d','sradsu','srads','albedo', \
                    'trads','tradsu', 'aclcov', 'var54']
    for i in check_vars:
        try:
            vv = f0.variables[i]
            vv = f1.variables[i]
        except:
            if  i == 'var54':
                print(" The variable var54 is q2m, change q2m to var54 for echam2ice functions. (cdo chname,q2m,var54)")

            raise ValueError('\n!!!!!! Variable: '+ i +"  is required. Stop !!!!!" )
    #----------

    lat = f1.variables['lat']
    lon = f1.variables['lon']
    nx = lon.shape[0]
    ny = lat.shape[0]


    dtype = 1
    if dtype == 1:
        SOLm = np.zeros((365, ny))
        daya = np.arange(1.78, 366.8)  # (2.3,367.3)
    else:
        SOLm = np.zeros((360, ny))
        daya = np.arange(0, 360)

    srad0dI = np.zeros((12, ny, nx))
    sradtest = np.zeros((12, ny, nx))

    w_srad0d = f1.variables['srad0d'][:]
    w_sradsu = f1.variables['sradsu'][:]
    w_srads = f1.variables['srads'][:]
    w_alb = f1.variables['albedo'][:]
    w_sradsd = w_srads - w_sradsu  # downward solar radiation at the surface
    w_tao = np.where( w_srad0d>0., w_sradsd / w_srad0d, 0.)  # else surface/TOA downward solar radiation #LU: set w_tao to 0. if w_srad0d is 0)
    #
    c_srad0d = f0.variables['srad0d'][:]
    c_sradsu = f0.variables['sradsu'][:]
    c_srads = f0.variables['srads'][:]
    c_alb = f0.variables['albedo'][:]  # albedo
    c_sradsd = c_srads - c_sradsu  # downward solar radiation at the surface
    c_tao = np.where( c_srad0d>0.,  c_sradsd / c_srad0d ,0.) # else surface/TOA downward solar radiation #LU: set c_tao to 0. if c_srad0d is 0)


    # subset of orbital parameters for a specified year in kyr (LGM is -21, careful: not yet defined for kyear>0!)
    orb = OrbitalTable.interp(kyear=kyr)
    orb2 = OrbitalTable.interp(kyear=0)
    orb_PI = {"ecc": 0.016724, "long_peri": 282.157, "obliquity": 23.4468}
    orb_LGM = {"ecc": 0.018994, "long_peri": 294.42, "obliquity": 22.949}
    # insolation values for the specific year for all days of a year and all latitudes of the atmospheric grid
    latt = np.zeros(96) ; latt[:] = lat[:]
    SOL = daily_insolation(lat=latt, day=daya, orb=orb, S0=1365.0, day_type=dtype)

    print("I:       The orbital information:")
    print(orb)
    # write to orbital.txt file (for run dEBM)
    #forb = open("orbital.txt", "w")
    #forb.write("COBLD_echam="+str(orb['obliquity'].values))
    #forb.close()
    obliquity = orb['obliquity'].values

    # compute monthly mean insolation from daily values
    if dtype == 1:
        days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        firstday = 0
    else:
        days_per_month = 30
        firstday = 360 - 80
        # firstday=0

    for i in range(0, 12):
        if dtype == 1:
            days = np.arange(firstday, firstday + days_per_month[i])  # dtype =1
            firstday = firstday + days_per_month[i]
        else:
            if (firstday + days_per_month) < 360:
                days = np.arange(firstday, firstday + days_per_month)
                firstday = firstday + days_per_month
            else:
                first_day_of_next_month = firstday + days_per_month - 360
                days = np.concatenate((range(firstday, 360), range(0, first_day_of_next_month)))
                firstday = first_day_of_next_month


        SOLm = SOL[:, days].mean(axis=1)
        for j in range(0, ny):
            srad0dI[i, j, :] = SOLm[j]

    # combine "warm" and "cold" fields according to index
    taoI = index * w_tao + (1 - index) * c_tao  # scaled transmissivity taoI
    albI = index * w_alb + (1 - index) * c_alb  # scaled albedo
    sradsdI = (
        taoI * srad0dI
    )  # downward surf. radiation from taoI and TOA insolation from climlab
    sradsuI = -albI * sradsdI  # but we actually only need UPWARD surf. radiation
    # and net surface radiation... I hope we can use albedo here...
    sradsI = ( 1 - albI) * sradsdI

    ##################### scaled variables and write to new output file ######################################
    fout = nc.Dataset(fileout,'w')
    # copy dimesions from refdata0
    for name,dimensions in f0.dimensions.items():
        fout.createDimension(name, len(dimensions) if not dimensions.isunlimited() else None)
        #
    print()
    print("I:       calculate data based on the index:", index)
    print(" and refdata0:", refdata0)
    print(" and refdata1:", refdata1)
    print("I:       With equation: index * ( refdata1 - refdata0 ) + refdata0 " )
    print("I:       the output file information are copy from refdata0" )
    #
    for name,varin in f0.variables.items():
        try:
            var1 = f1.variables[name][:]
            print("         Scaling variable:", name)
            outvar = fout.createVariable(name, varin.datatype, varin.dimensions )
            outvar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            #
            var0 = f0.variables[name][:]
            varnew = (var1-var0) * index + var0
            outvar[:] = varnew
        except:
            print("     ", name, "is not existed in other file, pass")
            pass
        #
    print("         But calculate srads, sradsu, srad0d based on insolation at current time")
    fout.variables['srads'][:] = sradsI
    fout.variables['sradsu'][:] = sradsuI
    fout.variables['srad0d'][:] = srad0dI

    print("         Cut negative aprl, aprc to 0")
    pp = fout.variables['aprl'][:]  ;  ppnew = np.where(pp<0., 0., pp)
    fout.variables['aprl'][:] = ppnew
    #
    pp = fout.variables['aprc'][:]  ;  ppnew = np.where(pp<0., 0., pp)
    fout.variables['aprc'][:] = ppnew

    print("         Done." )

    f0.close()
    f1.close()
    fout.close()
    return  obliquity


###############################################################################
#   new data with chosen index
###############################################################################
def get_index_data(theindex, refdata0, refdata1, fileout='data.nc' ):
    print()
    print("I:       calculate data based on the index:", theindex)
    print(" and refdata0:", refdata0)
    print(" and refdata1:", refdata1)
    print("I:       With equation: index * ( refdata1 - refdata0 ) + refdata0 " )
    print("I:       the output file information are copy from refdata0" )
    #
    f0 = nc.Dataset(refdata0,'r')
    f1 = nc.Dataset(refdata1,'r')
    fout = nc.Dataset(fileout,'w')
    # copy dimesions from refdata0
    for name,dimensions in f0.dimensions.items():
        fout.createDimension(name, len(dimensions) if not dimensions.isunlimited() else None)
        #
    # copy and calcualte variables
    # copy variables attributes from refdata0
    for name,varin in f0.variables.items():
        print("I:       processing variable:", name)
        outvar = fout.createVariable(name, varin.datatype, varin.dimensions )
        outvar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        #
        var0 = f0.variables[name][:]
        var1 = f1.variables[name][:]
        varnew = (var1-var0) * theindex + var0
        outvar[:] = varnew
        #
    print("I:       Done." )
    return


###############################################################################
#   interpolate FESOM to PISM grids
###############################################################################

def initial_FESOM2PISMgrids(**opts):
    #output:   opts['PISM_file_ocean_forcing_ref0']
    #          opts['PISM_file_ocean_forcing_ref1']
    print()
    print("I:       For initial: Interpolate FESOM file to PISM grids. ")
    print("I:       The defined settings: ")
    for aa, bb in opts.items():
        print("         ",aa,"=====>",bb)

    # get the ref0 file
    print("I:       Processing the ref0 data ")
    function_FESOM2ice(fintimes=[0000,], ntimeperfile=opts['FESOM_ntimeperfile'], version_FESOM=opts['FESOM_version'],
                meshpath=opts['FESOM_file_ref0_mesh'] , abg=opts['FESOM_file_ref0_mesh_abg'],
                findir=opts['ocean_file_dir'], fout='ocean_file_for_ice_ref0.nc',
                ntimemean=opts['FESOM_ntimeperfile'], timestartout=[1950,12,31,0,0], dtout=1 ,dtoutunit='years',
                depthmax=opts['FESOM_depth_max'], depthmin=opts['FESOM_depth_min'], depthmean=False,
                ifrawfesom=False, fnameonly=opts['ocean_filename_ref0'] )
    function_ocean2PISM(fin='ocean_file_for_ice_ref0.nc',ficegrid=opts['PISM_file_griddes'], fout=opts['PISM_file_ocean_forcing_ref0'],
                switch_PISM=opts['PISM_ocean_switch'], version_PISM=opts['PISM_version'],fpicobasin=opts['PISM_file_picobasins'] )
    function_ocean_landmask(finout=opts['PISM_file_ocean_forcing_ref0'], fslm2d=opts['landseamask_file_ref0'], fgriddes=opts['PISM_file_griddes'])

    # get the ref1 file
    print("I:       Processing the ref1 data ")
    function_FESOM2ice(fintimes=[0000,], ntimeperfile=opts['FESOM_ntimeperfile'], version_FESOM=opts['FESOM_version'],
                meshpath=opts['FESOM_file_ref1_mesh'] , abg=opts['FESOM_file_ref1_mesh_abg'],
                findir=opts['ocean_file_dir'], fout='ocean_file_for_ice_ref1.nc',
                ntimemean=opts['FESOM_ntimeperfile'], timestartout=[1950,12,31,0,0], dtout=1 ,dtoutunit='years',
                depthmax=opts['FESOM_depth_max'], depthmin=opts['FESOM_depth_min'], depthmean=False,
                ifrawfesom=False, fnameonly=opts['ocean_filename_ref1'] )
    function_ocean2PISM(fin='ocean_file_for_ice_ref1.nc',ficegrid=opts['PISM_file_griddes'], fout=opts['PISM_file_ocean_forcing_ref1'],
                switch_PISM=opts['PISM_ocean_switch'], version_PISM=opts['PISM_version'],fpicobasin=opts['PISM_file_picobasins'] )
    function_ocean_landmask(finout=opts['PISM_file_ocean_forcing_ref1'], fslm2d=opts['landseamask_file_ref1'], fgriddes=opts['PISM_file_griddes'])

    #clean intermeidate files
    os.system("rm  ocean_file_for_ice_ref0.nc  ocean_file_for_ice_ref1.nc  ")
    print("I:       Done. ")
    return
