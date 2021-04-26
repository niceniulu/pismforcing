#!/usr/bin/env python3
# coding: utf-8


import sys , os
import netCDF4   as nc
import numpy as np


import matplotlib.pyplot as plt



def downscale_ECHAM2PISMgrids(fecham='select_ECHAM6.nc',
                              fpismgrid='pismr_nhem_20km.griddes',
                              fpismcurrent='20km_lrno_pismr_out_-021100--021001.nc',
                              interpmethod='bil' ,
                              downscaleT='linear', factorT=-5,
                              downscaleP='desertification', factorP=None,
                              fout='downscaled.nc' ):
    '''
    Interpolate from ECHAM grids to PISM grids.
    '''
    #------------------------------------------------
    # regrid from ECHAM6 grid to PISM grid:  fout
    #------------------------------------------------
    print("A2I:     regrid ECHAM6 to PISM grids using cdo:", interpmethod, "| output: ", fout )
    os.system("cdo   remap"+interpmethod+","+fpismgrid+"   "+fecham+"    " + fout )
    #

    #--------------------------
    # read data
    #--------------------------
    # read regrid (fout, echam6)
    ff = nc.Dataset(fout,'r')
    usurf_old = ff.variables['orog'][:]
    precip = ff.variables['aprt'][:]
    temp2 = ff.variables['temp2'][:]
    ff.close()

    # read pism current file
    fpc = nc.Dataset(fpismcurrent,'r')
    try:
        # pism output file with time dim. [time,,]
        thk = fpc.variables['thk'][0,:,:]
        topg = fpc.variables['topg'][0,:,:]
        sealevel = fpc.variables['effective_sea_level_elevation'][0,:,:]
    except:
        # pism inital file without time dim and sea level elevation.
        thk = fpc.variables['thk'][:,:]
        topg = fpc.variables['topg'][:,:]
        sealevel = 0.

    fpc.close()
    tmp = thk + topg - sealevel
    usurf_new = np.where(tmp>0., tmp, 0.)

    diff_usurf = usurf_new - usurf_old


    #------------------------------
    # Temperature downscale
    #------------------------------
    print("A2I:     Downscaling temperature with option: ", downscaleT, ", and number",factorT, "K/km")
    if  downscaleT == 'linear' :
        if factorT > 0:
            raise ValueError('!!!!! DownscaltT: '+str(factorT) +' should be negative (K/km), Stop')
        elif  factorT < 0 and factorT > -1 :
            raise ValueError('!!!!! DownscaltT: '+str(factorT) +' should be with unit K/km, Stop')
        #

        diff_Temp = diff_usurf * factorT/1000.
        T_downscaled = temp2 + diff_Temp
    # elif : to be implemented
    else :
        print("A2I:     no downscaling Temperature applied. ")
        T_downscaled = temp2

    #

    #------------------------------
    # Precipitation downscale
    #------------------------------
    print("A2I:     Downscaling precipitation with option: ", downscaleP, ", need factor:",factorP)
    if downscaleP == 'desertification':
        # downscaling precipitation following desertification scheme -> Ziemen et al, 2014
        #  P_downscaled = P * e ** (-gamma * max(0,(max(2000,Z_hi) - max(2000,Z_lo))))
        #
        thresh_elevation = 2000  # m
        gamma = -0.00069315
        usurf_hires = np.where(usurf_new > thresh_elevation, usurf_new, thresh_elevation )
        usurf_lores = np.where(usurf_old > thresh_elevation, usurf_old, thresh_elevation )
        tmp = usurf_hires - usurf_lores
        P_downscaled = precip * np.exp( gamma * np.where(tmp>0., tmp, 0. ))
    # elif : to be implemented
    else :
        print("A2I:     No downscaling Precipitation applied. ")
        P_downscaled = precip

    #------------------------------
    # Write to output
    #------------------------------
    print("A2I:     Write to output file: ", fout ," | new vars: air_temp, precipitation ")

    f = nc.Dataset(fout,'r+')
    temp2 = f.variables['temp2']
    aprt = f.variables['aprt']
    otype = temp2.datatype
    odims = temp2.dimensions

    # downscalsed temp2
    temp2out = f.createVariable('air_temp', otype, odims)
    # copy attributes from original file
    temp2out.setncatts({i: temp2.getncattr(i)  for i in temp2.ncattrs()})
    temp2out[:] = T_downscaled

    # downscalsed precip
    precipout = f.createVariable('precipitation', otype, odims)
    # copy attributes from original file
    precipout.setncatts({i: aprt.getncattr(i)  for i in aprt.ncattrs()})
    precipout[:] = P_downscaled

    f.close()

    print('A2I:     Done. ')

    return 0




def select_vars_for_PISM(fin='dummy_ECHAM6.nc',
                         schemeSMB='PDD',     # PDD | dEBM
                         lakestore=True,   # if run lakestore, then output variable: evaporation
                         fout='select_ECHAM6.nc'):
    '''
    Select necessary variables for generating PISM forcing.
    '''
    #
    g = 9.81  # m/s^2
    sigma = 5.6703744e-8

    with nc.Dataset(fin,'r') as src, nc.Dataset(fout,'w') as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name in src.dimensions.keys():
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name][:] = src[name][:]
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
        #

        print("A2I:     Select and modify vars for SMB scheme: ", schemeSMB)
        geosp = src.variables['geosp']
        aprl = src.variables['aprl']
        aprc = src.variables['aprc']
        temp2 = src.variables['temp2']
        types = temp2.datatype
        dims = temp2.dimensions
        # calculate
        orog = geosp[:]/g
        aprt = aprl[:] + aprc[:]
        #
        dst.createVariable('temp2', types, dims)
        dst['temp2'][:] = temp2[:]
        dst['temp2'].long_name = "2m temperature"
        dst['temp2'].units = "K"
        #
        dst.createVariable('aprt', types, dims)
        dst['aprt'][:] = aprt
        dst['aprt'].long_name = "total precipitation"
        dst['aprt'].units = "kg m-2 second-1"
        #
        dst.createVariable('orog', types, dims)
        dst['orog'][:] = orog
        dst['orog'].long_name = "orography"
        dst['orog'].units = "m"

        if schemeSMB == 'PDD':
            print("         Done.  ")
        elif schemeSMB == 'dEBM':
            srads = src.variables['srads']
            sradsu = src.variables['sradsu']
            aclcov = src.variables['aclcov']
            trads = src.variables['trads']
            tradsu = src.variables['tradsu']
            srad0d = src.variables['srad0d']
            var54 = src.variables['var54']
            # calculate
            swd = srads[:] - sradsu[:]
            cc = aclcov[:]
            emiss = (trads[:] - tradsu[:])/(sigma * temp2[:]**4)
            tau = np.where(srad0d[:]<5, 0.5, swd/srad0d[:] )
            TOAswd = srad0d[:]
            q2m = var54[:]
            #
            dst.createVariable('swd', types, dims)
            dst['swd'][:] = swd
            dst['swd'].long_name = "surface downward shortwave radiation "
            dst['swd'].units = "W/m2"
            #
            dst.createVariable('emiss', types, dims)
            dst['emiss'][:] = emiss
            dst['emiss'].long_name = "emissivity "
            dst['emiss'].units = "frac"
            #
            dst.createVariable('cc', types, dims)
            dst['cc'][:] = cc
            dst['cc'].long_name = "cloud cover "
            dst['cc'].units = "frac"
            #
            dst.createVariable('q2m', types, dims)
            dst['q2m'][:] = q2m
            dst['q2m'].long_name = "specific humidity at 2m "
            dst['q2m'].units = "frac"
            #
            dst.createVariable('tau', types, dims)
            dst['tau'][:] = tau
            dst['tau'].long_name = "atm. transmissivity  "
            dst['tau'].units = "frac"
            #
            dst.createVariable('TOAswd', types, dims)
            dst['TOAswd'][:] = TOAswd
            dst['TOAswd'].long_name = "TOA shortwave downward radiation  "
            dst['TOAswd'].units = "W/m2"
            #
            print("         Done. ")
        else:
            print("A2I:     unknown SMB scheme", schemeSMB, ", generate orog, temp2, aprt only")


        if  lakestore:
            print("A2I:     select evap for lakestore. ")
            dst.createVariable('evap', types, dims)
            dst['evap'][:] = src.variables['evap'][:]
            dst['evap'].long_name = "evaporation"
            dst['evap'].units = "kg m-2 s-1"
    #
    return



def run_dEBM(fin='downscaled.nc',
             obliquity=None,
             exedEBM=None):

    ff = nc.Dataset(fin,'r')

    print("A2I:     check required variables ")
    # -- checking variables -----
    vars = ['lat','lon','air_temp', 'precipitation','swd','emiss','cc','q2m','tau','TOAswd']
    for i in vars:
        try:
            vv = ff.variables[i]
        except:
            raise ValueError('!!!!!  Variable: '+ i +"  is missing. Stop !!!!!" )

    if obliquity is None:
        raise ValueError("!!!!! Need to specify the obliquity value for dEBM, Stop. ")
    #
    if exedEBM is None:
        raise ValueError("!!!!! Need to specify the path for exe dEBM, Stop. ")
    #
    # write to namelist.debm
    print("A2I:     generate namelist.debm ")
    fnml = open("namelist.debm", "w")
    fnml.write("&runctl\n")
    fnml.write("lresume=.false.\n")
    fnml.write("use_shortwave_radiation_TOA=.true.\n")
    fnml.write("use_mask=.false.\n")
    fnml.write("debug_switch=.false.\n")
    fnml.write("debug_lon=241\n")
    fnml.write("debug_lat=107\n")
    fnml.write("debug_mon=12\n")
    fnml.write("debug_year=1\n")
    fnml.write("/\n\n")
    fnml.write("&debm\n")
    fnml.write("filename_in='" + fin +"'\n")
    fnml.write("precipitation_varname='precipitation'\n")
    fnml.write("temperature_varname='air_temp'\n")
    fnml.write("shortwave_radiation_downward_varname='swd'\n")
    fnml.write("shortwave_radiation_TOA_varname='TOAswd'\n")
    fnml.write("cloud_cover_varname='cc'\n")
    fnml.write("emissivity_varname='emiss'\n")
    fnml.write("transmissivity_varname='tau'\n")
    fnml.write("mapping_varname='thk'\n")
    fnml.write("longitude_varname='lon'\n")
    fnml.write("latitude_varname='lat'\n")
    fnml.write("time_varname='time'\n")
    fnml.write("stddev=3.5\n")
    fnml.write("obliquity=" + str(obliquity) + "\n")
    fnml.write("/\n")
    fnml.close()

    # run dEBM
    print("A2I:     running dEBM, ./dEBMmain  namelist.debm ")
    os.system("cp   "+exedEBM+"    ./dEBMmain")
    os.system("./dEBMmain  namelist.debm  ")

    print("A2I:     Done. | output: surface_mass_balance.nc ")
    return



def generate_atmos_forcing_for_PISM(fin='downscaled.nc',
                                    fout='PISM.nc',
                                    scheme='PDD',
                                    fileplus1=None,
                                    timefreq = "12month",
                                    lakestore = False):

    print("A2I:     generate atmos. forcing for PISM with scheme: ", scheme)

    vars = 'air_temp,precipitation'
    # if lakestore, add evaporation
    if lakestore:
        vars = vars + ',evap'


    if scheme == 'PDD':
        print('A2I:     pdd')
        os.system(" ncks  -O -v  " + vars + "    "+fin+"  "+ fout )

        with nc.Dataset(fout,'r+' ) as ff:
            #
            air_temp = ff.variables['air_temp']
            air_temp.standard_name="air_temperature"
            air_temp.long_name = "Surface air temperature"
            #
            precipitation = ff.variables['precipitation']
            precipitation.standard_name = "lwe_precipitation_rate"
            precipitation.long_name = "Total Precipitation"
            #
    elif scheme == 'dEBM':
        # check if the surface mass balance file existed.
        if fileplus1 is None:
            raise ValueError("!!!!! Need to specify outputfile from dEBM (fileplus1), Stop. ")

        print("A2I:     debm ")

        os.system(" ncks  -O -v   " + vars + "      "+fin+"  "+ fout )

        with nc.Dataset(fout,'r+' ) as ff, nc.Dataset(fileplus1,'r') as fp:
            # fileplus1 (SMB)
            smb = fp.variables['SMB']

            # file output
            ff = nc.Dataset(fout,'r+' )
            air_temp = ff.variables['air_temp']
            precip = ff.variables['precipitation']
            # ice surface temperature
            tempice = np.where(air_temp[:]<273.15, air_temp[:], 273.15)
            # runoff (no ice mask)
            runoff = precip[:] - smb[:]

            # write to output file
            dtype = air_temp.datatype
            dims = air_temp.dimensions
            print(dtype, dims)
            ff.createVariable('ice_surface_temp', dtype, dims)
            ff['ice_surface_temp'][:] = tempice[:]
            ff['ice_surface_temp'].long_name = "ice surface temperature  "
            ff['ice_surface_temp'].units = "K"

            ff.createVariable('climatic_mass_balance', dtype, dims)
            ff['climatic_mass_balance'][:] = smb[:]
            ff['climatic_mass_balance'].long_name = "climatic surface mass balance  "
            ff['climatic_mass_balance'].units = "kg m-2 s-1"

            ff.createVariable('water_input_rate', dtype, dims)
            ff['water_input_rate'][:] = runoff
            ff['water_input_rate'].long_name = "runoff (precip-surface_mass_balance)  "
            ff['water_input_rate'].units = "kg m-2 s-1"

    # add time axis information
    print("A2I:     add time and time bounds dimensions")
    with nc.Dataset(fout,'r+' ) as ff:
        time = ff.variables['time']
        time.standard_name = 'time'
        time.long_name = 'time'
        time.units = 'days since 01-01-01 00:00:00'
        time.calendar = '365_day'
        time.bounds = 'time_bnds'
        #
        ff.createDimension('bnds',2)
        timebnds = ff.createVariable('time_bnds',np.float64, ('time','bnds'))
        if timefreq == "12month":
            ntime = 12
            time[:] = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5,319, 349.5]
            timebnds[:ntime,:]=[0, 31, 31, 59, 59, 90, 90, 120, 120, 151, 151, 181, 181, 212, 212, 243, 243, 273, 273, 304, 304, 334, 334, 365.]
        else:
            raise ValueError("..... not 1-year monthly data, stop.")

    print("A2I:     Done. ")
    return
