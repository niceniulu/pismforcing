#!/usr/bin/env python3
# coding: utf-8

import netCDF4  as  nc
import numpy  as  np
import os
import calendar
import argparse



def  lapserate_albedo_correction(fclimate,
                                 fusurf,
                                 fpismcurrent,
                                 if_albedo=True,
                                 interpmethod='bil' ,
                                 downscaleT='linear', factorT=-5,
                                 downscaleP='desertification', factorP=None,
                                 fout='out.nc' ):

    # read climate fields:
    ff = nc.Dataset(fclimate, 'r')
    temp2 = ff.variables['air_temp'][:]
    precip = ff.variables['precipitation'][:]
    ff.close()

    # read the surface height for lapse rate correction
    ff = nc.Dataset(fusurf, 'r')
    usurf_old = ff.variables['usurf'][:]
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

    #------------------------------
    # albedo correction: if there is ice, and temperature >0., there set T to 0.
    #------------------------------
    if  if_albedo :
        print("    apply albedo correction")
        condition = (thk>1.) & (T_downscaled > 0.) # where there is ice, and tempearute > 0.
        T_downscaled = np.where( condition , 0. , T_downscaled )

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
    os.system(" cp " +fclimate+"   "+fout)

    f = nc.Dataset(fout,'r+')
    f.variables['air_temp'][:] = T_downscaled
    f.variables['precipitation'][:] = P_downscaled

    f.close()






if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Input variables')
    parser.add_argument("fclimate",type=str, help='Input climate forcing file on PISM grids')
    parser.add_argument("fusurf",type=str, help='Input usurf file')
    parser.add_argument("fpismcurrent",type=str, help='Input PISM output file')
    parser.add_argument("--albedo", action='store_true' , help='if apply albedo correction')
    parser.add_argument("--interpmethod", type=str,  default='bil', help='interpolation method')
    parser.add_argument("--downscaleT", type=str,  default='linear', help='downscale Temperature method')
    parser.add_argument("--factorT", type=float, default=-5., help='Temperature downscale factor dT/dh, K/km')
    parser.add_argument("--downscaleP", type=str,  default='desertification', help='downscale Precipitation method')
    parser.add_argument("--factorP", type=str, default='None', help='Precipitation downscale factor')

    parser.add_argument("--fout",  type=str, default='out.nc', help='output file name')

    arguments = parser.parse_args()
    print(arguments)

    fclimate = arguments.fclimate # 'for_pism_1x_echam.nc'
    fpismcurrent = arguments.fpismcurrent  # 'out_Cre_1xCO2_100000.nc'
    fusurf = arguments.fusurf   # 'usurf_Antarctica_16km.nc'

    lapserate_albedo_correction(fclimate, fusurf, fpismcurrent, \
                                if_albedo=arguments.albedo, interpmethod=arguments.interpmethod , \
                                downscaleT=arguments.downscaleT, factorT=arguments.factorT, \
                                downscaleP=arguments.downscaleP, factorP=arguments.factorP, fout=arguments.fout )
